#!/usr/bin/env python3
"""
Check SDRscan results on a simulated dataset (tests/simulate_dataset.py)
against what was planted: an XY sex-determining region on chr_03:300-700 kb,
a Y deletion at 600-700 kb (lower male depth) and a Y-only insertion.
Prints PASS / FAIL / SKIP per check and exits non-zero if anything fails.
"""

import argparse
import glob
import os
import sys

import pandas as pd

SDR = ("chr_03", 300_000, 700_000)
LOW_DEPTH = ("chr_03", 600_000, 700_000)


def overlaps(chrom, start, end, region, slack=20_000):
    c, s, e = region
    return chrom == c and start < e + slack and end > s - slack


class Checks:
    def __init__(self):
        self.rows = []

    def add(self, name, ok, detail=""):
        status = "SKIP" if ok is None else ("PASS" if ok else "FAIL")
        self.rows.append((status, name, detail))
        print(f"{status:5} {name}" + (f"  ({detail})" if detail else ""))

    def failed(self):
        return any(r[0] == "FAIL" for r in self.rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--pilot", type=int, default=0)
    ap.add_argument("--males", type=int, required=True)
    ap.add_argument("--females", type=int, required=True)
    a = ap.parse_args()
    R = a.results
    c = Checks()
    ploidy = int(dict(l.split("\t") for l in open(f"{a.data}/expected.txt").read().split("\n") if l)["ploidy"])

    # Stage 1-3
    bams = glob.glob(f"{R}/01_align/bam/*.RG.bam")
    n_samples = len(pd.read_csv(f"{a.data}/samples.tsv", sep="\t"))
    c.add("stage 1: one BAM per sample", len(bams) == n_samples, f"{len(bams)}/{n_samples}")
    cw = f"{R}/02_vcf/called_with.tsv"
    called = pd.read_csv(cw, sep="	", header=None) if os.path.exists(cw) else pd.DataFrame()
    c.add("stage 2: joint VCF called at the declared ploidy",
          os.path.exists(f"{R}/02_vcf/joint.vcf.gz") and len(called) == n_samples and set(called[1]) == {ploidy},
          f"{len(called)} samples, ploidy {sorted(set(called[1])) if len(called) else '-'}")
    gt = f"{R}/03_merged/merged.genotypes.tsv"
    n_sites = sum(1 for _ in open(gt)) if os.path.exists(gt) else 0
    c.add("stage 3: merged genotype table", n_sites > 1000, f"{n_sites} SNPs")
    if n_sites:
        g = pd.read_csv(gt, sep="\t", header=None, nrows=5000, dtype=str).iloc[:, 2:]
        is_missing = g.apply(lambda col: col.str.startswith("."))
        called = g[~is_missing].stack()
        n_alleles = called.str.replace("|", "/", regex=False).str.split("/").str.len()
        c.add("stage 3: genotypes have the declared ploidy", bool((n_alleles == ploidy).all()), called.iloc[0])
        miss = float(is_missing.to_numpy().mean())
        c.add("stage 3: few missing genotypes (joint calling)", miss < 0.1, f"{miss:.1%} missing")

    if not a.pilot:
        # Stage 5
        sites = pd.read_csv(f"{R}/05_snp_windows/relaxed/sex_specific_sites.tsv", sep="\t")
        ms = sites[sites["CLASS"] == "male_specific"]
        in_sdr = ms.apply(lambda r: overlaps(r.CHROM, r.POS, r.POS, SDR, 0), axis=1).mean() if len(ms) else 0
        c.add("stage 5: male-specific SNPs found", len(ms) > 50, f"{len(ms)} sites")
        c.add("stage 5: male-specific SNPs fall in the SDR", in_sdr > 0.8, f"{in_sdr:.0%} in SDR")
        # With few samples, variants shared by chance look sex-specific in both
        # directions; in an XY system female-specific sites should stay at that
        # background, i.e. about the number of male-specific sites outside the SDR.
        fs = (sites["CLASS"] == "female_specific").sum()
        bg = int(round(len(ms) * (1 - in_sdr))) if len(ms) else 0
        c.add("stage 5: female-specific SNPs only at background level (XY system)", fs <= max(10, 1.5 * bg),
              f"{fs} female-specific vs {bg} male-specific outside the SDR")
        iv = pd.read_csv(f"{R}/05_snp_windows/relaxed/merged_candidate_intervals.tsv", sep="\t")
        top = iv[iv["SEX"] == "male_specific"].head(1)
        c.add("stage 5: top candidate interval is the SDR",
              len(top) == 1 and overlaps(top.CHROM.iloc[0], top.START.iloc[0], top.END.iloc[0], SDR),
              f"{top.CHROM.iloc[0]}:{top.START.iloc[0]}-{top.END.iloc[0]}" if len(top) else "none")
        dose = sorted(os.path.basename(p) for p in glob.glob(f"{R}/05_snp_windows/relaxed_dose*"))
        c.add("stage 5: dosage-class runs 1..ploidy-1", len(dose) == ploidy - 1, ", ".join(dose))
        d1 = pd.read_csv(f"{R}/05_snp_windows/relaxed_dose1/sex_specific_sites.tsv", sep="\t")
        n1 = (d1["CLASS"] == "male_specific").sum()
        c.add("stage 5: Y SNPs are single-dose (dose1 run finds them)", n1 > 0.5 * len(ms), f"{n1} sites")
        png = glob.glob(f"{R}/05_snp_windows/relaxed/sex_specific_snp_windows.png")
        c.add("stage 5: genome-wide plot", bool(png))

        # Stage 6
        cdir = f"{R}/06_candidates/chr_03_300000_700000"
        res = dict(l.split("\t") for l in open(f"{cdir}/permutation.result.txt").read().split("\n") if l)
        p, minp = float(res["p_value"]), float(res["min_possible_p"])
        if minp < 0.05:
            c.add("stage 6: permutation test significant", p < 0.05, f"p={p:.3g}, floor {minp:.3g}")
        else:
            c.add("stage 6: permutation p at the design's floor", abs(p - minp) < 1e-9, f"p={p:.3g}, floor {minp:.3g}")
        c.add("stage 6: strict-overlap table", os.path.exists(f"{R}/06_candidates/strict_overlap_counts.tsv"))
        pca = glob.glob(f"{cdir}/pca_presence.scores.tsv")
        if pca:
            sc = pd.read_csv(pca[0], sep="\t")
            m, f = sc[sc.label == "male"].PC1, sc[sc.label == "female"].PC1
            c.add("stage 6: PCA separates the sexes on PC1", m.max() < f.min() or m.min() > f.max(),
                  f"male {m.min():.2f}..{m.max():.2f}, female {f.min():.2f}..{f.max():.2f}")
        else:
            c.add("stage 6: PCA separates the sexes on PC1", None, "PCA skipped (too few sites)")
        if ploidy > 2:
            c.add("stage 6: dosage PCA for polyploids", bool(glob.glob(f"{cdir}/pca_dosage.scores.tsv")))
        ins = open(f"{cdir}/inspect.summary.txt").read()
        c.add("stage 6: genotype inspection", "Strong male-biased" in ins)

        # Stage 7
        loo = glob.glob(f"{R}/07_leave_one_out_*/summary.tsv")
        if a.males >= 3:
            t = pd.read_csv(loo[0], sep="\t") if loo else pd.DataFrame()
            c.add("stage 7: SDR survives every leave-one-out", len(t) == a.males and (t.region_found == "yes").all(),
                  f"{(t.region_found == 'yes').sum() if len(t) else 0}/{a.males}")
        else:
            c.add("stage 7: leave-one-out skipped (< 3 males)", not loo)

        # Stage 8-9
        for ref in ("F", "M"):
            f = f"{R}/08_kmer/MSK/windows/genome.MSK_density.{ref}.bed"
            if os.path.exists(f):
                w = pd.read_csv(f, sep="\t", header=None, names=["CHROM", "START", "END", "COUNT"])
                top = w.sort_values("COUNT", ascending=False).head(5)
                ok = all(overlaps(r.CHROM, r.START, r.END, ("chr_03", 300_000, 720_000)) for r in top.itertuples())
                c.add(f"stage 8: male-specific k-mers peak in the SDR ({ref} reference)", ok,
                      ", ".join(f"{r.CHROM}:{r.START}" for r in top.itertuples()))
            else:
                c.add(f"stage 8: male-specific k-mers peak in the SDR ({ref} reference)", False, "no output")
        c.add("stage 9: k-mer overlap summary", os.path.exists(f"{R}/09_kmer_overlap/overlap_summary.tsv"))
        c.add("stage 11: Manhattan plots", len(glob.glob(f"{R}/11_figures/*manhattan.png")) >= 2)
    else:
        for st in ("05_snp_windows", "06_candidates", "08_kmer"):
            c.add(f"pilot: stage {st[:2]} skipped", not os.path.exists(f"{R}/{st}"))

    # Stage 10-11 dotplots
    coords = glob.glob(f"{R}/10_dotplots/*/*.coords")
    c.add("stage 10: nucmer alignment", bool(coords))
    if not a.pilot:
        c.add("stage 11: genome dotplot + correspondence table",
              bool(glob.glob(f"{R}/11_figures/*_genome.png")) and bool(glob.glob(f"{R}/11_figures/*correspondence_best.tsv")))
        summ = [p for p in glob.glob(f"{R}/11_figures/*.summary.tsv")]
        rev = [pd.read_csv(p, sep="\t") for p in summ]
        n_rev = sum(int(t["reverse"].fillna(0).sum()) for t in rev)
        c.add("stage 11: Y inversion seen as reverse alignment in the chromosome dotplot", n_rev > 0,
              f"{n_rev} reverse alignments in {len(summ)} chromosome plots")
        corr = glob.glob(f"{R}/11_figures/*correspondence_best.tsv")
        if corr:
            t = pd.read_csv(corr[0], sep="\t")
            c.add("stage 11: chromosomes correspond 1:1", (t.chr_ref == t.chr_qry).all(),
                  ", ".join(f"{a_}->{b_}" for a_, b_ in zip(t.chr_ref, t.chr_qry)))

    # Stage 12
    cand = pd.read_csv(f"{R}/12_depth_het/depth_het.candidates.tsv", sep="\t")
    dc = cand[cand.evidence == "depth"]
    hit = dc[dc.apply(lambda r: overlaps(r.CHROM, r.START, r.END, LOW_DEPTH), axis=1)]
    c.add("stage 12: lower male depth over the Y deletion", len(hit) > 0 and (hit.mean_effect < 0).all(),
          f"{len(dc)} depth candidates; " + "; ".join(f"{r.CHROM}:{r.START}-{r.END} ({r.mean_effect:.2f})" for r in dc.head(3).itertuples()))
    off = dc[~dc.apply(lambda r: overlaps(r.CHROM, r.START, r.END, SDR), axis=1)]
    c.add("stage 12: no depth candidates outside the SDR", len(off) == 0, f"{len(off)} elsewhere")
    strong = cand[cand.support == "strong"]
    het_sdr = strong[(strong.evidence == "het") & strong.apply(lambda r: overlaps(r.CHROM, r.START, r.END, SDR), axis=1)]
    c.add("stage 12: higher male heterozygosity across the SDR", len(het_sdr) > 0 and (het_sdr.mean_effect > 0).all(),
          "; ".join(f"{r.CHROM}:{r.START}-{r.END} ({r.mean_effect:.2f})" for r in het_sdr.itertuples()) or "none")
    off_strong = strong[~strong.apply(lambda r: overlaps(r.CHROM, r.START, r.END, SDR), axis=1)]
    c.add("stage 12: no strong (multi-window) candidates outside the SDR", len(off_strong) == 0,
          f"{len(off_strong)} strong elsewhere, {len(cand) - len(strong)} weak single/double windows")

    # Stage 13
    sx = f"{R}/13_sexing/sex_assignment.tsv"
    if os.path.exists(sx):
        s = pd.read_csv(sx, sep="\t")
        truth = pd.read_csv(f"{a.data}/truth.tsv", sep="\t").set_index("sample")["true_sex"]
        s["truth"] = s["sample"].map(truth)
        unk = s[s.recorded_sex == "U"]
        if len(unk):
            c.add("stage 13: unknown samples sexed correctly", (unk.called_sex == unk.truth).all(),
                  ", ".join(f"{r.sample}: {r.called_sex} (true {r.truth})" for r in unk.itertuples()))
            c.add("stage 13: unknown samples sexed with medium or high confidence",
                  unk.confidence.isin(["medium", "high"]).all(),
                  ", ".join(f"{r.sample}: {r.confidence}" for r in unk.itertuples()))
        kn = s[s.recorded_sex != "U"]
        c.add("stage 13: recorded sexes confirmed", (kn.check == "agrees").all(),
              f"{(kn.check == 'agrees').sum()}/{len(kn)} agree")
    else:
        c.add("stage 13: sex assignment output", False, "missing")

    n_fail = sum(r[0] == "FAIL" for r in c.rows)
    print(f"\n{len(c.rows) - n_fail}/{len(c.rows)} checks passed")
    sys.exit(1 if c.failed() else 0)


if __name__ == "__main__":
    main()
