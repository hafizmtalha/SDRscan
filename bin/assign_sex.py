#!/usr/bin/env python3
"""
Assign sex to samples from genomic evidence at known sex-linked loci, and
check that samples of recorded sex agree with their genomes.

Evidence (use any combination; at least one is required):

  depth   mean normalised read depth over known sex-linked regions (--regions
          BED + --depth table from depth_het_scan.py). With a female (XX/ZW-W
          free) reference, X/Z-like regions show half depth in the heterogametic
          sex.
  snp     share of known male-specific SNPs (--markers, from
          classify_sex_specific_snps.py) where the sample carries ALT, minus the
          share of known female-specific SNPs where it does.
  kmer    share of male-specific k-mers the sample contains minus the share of
          female-specific k-mers (--kmer-fractions: sample, msk_fraction,
          fsk_fraction).

Each evidence becomes a score where positive = male-like and negative =
female-like, scaled so that +1 / -1 is the typical value of a male / female:
  - if the sample sheet has samples of known sex (M/F) with data, the scale is
    calibrated on them (midpoint between the two sex means = 0);
  - otherwise theoretical expectations are used (depth: heterogametic 0.5 vs
    homogametic 1.0, which needs --system; snp/kmer: 0 = no preference).

The combined score is the mean of the available evidence scores. Calls:
  M / F       |combined| >= --call-threshold and no evidence disagrees strongly
  uncertain   weak or conflicting evidence

Caveat: SNP and k-mer markers discovered from the known-sex samples will of
course separate those same samples, so for them only the depth evidence is an
independent check (column marker_evidence_circular).
"""

import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gtutils import GenotypeParser, load_ploidy, read_list


def depth_values(depth_file, regions_file):
    d = pd.read_csv(depth_file, sep="\t")
    reg = pd.read_csv(regions_file, sep="\t", header=None, comment="#").iloc[:, :3]
    reg.columns = ["CHROM", "START", "END"]
    mask = np.zeros(len(d), dtype=bool)
    for _, r in reg.iterrows():
        mask |= (d["CHROM"] == r["CHROM"]) & (d["END"] > r["START"]) & (d["START"] < r["END"])
    if not mask.any():
        sys.exit("ERROR: no depth windows overlap the --regions")
    samples = [c for c in d.columns if c not in ("CHROM", "START", "END")]
    print(f"depth: {int(mask.sum())} windows in sex-linked regions")
    return d.loc[mask, samples].mean(axis=0)


def snp_values(markers_file, genotypes, vcf_samples_file, ploidy_args):
    markers = pd.read_csv(markers_file, sep="\t")
    cls = {(c, int(p)): k for c, p, k in zip(markers["CHROM"], markers["POS"], markers["CLASS"])}
    vcf_samples = read_list(vcf_samples_file)
    parser = GenotypeParser(load_ploidy(vcf_samples, *ploidy_args))
    hits = {s: {"male_specific": [0, 0], "female_specific": [0, 0]} for s in vcf_samples}
    with open(genotypes) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            k = cls.get((parts[0], int(parts[1])))
            if k is None:
                continue
            for s, gt in zip(vcf_samples, parts[2:]):
                d = parser.dosage(gt, s)
                if d is None:
                    continue
                hits[s][k][0] += d > 0
                hits[s][k][1] += 1
    parser.report(sys.stderr)
    n_m = (markers["CLASS"] == "male_specific").sum()
    n_f = (markers["CLASS"] == "female_specific").sum()
    print(f"snp: {n_m} male-specific and {n_f} female-specific marker sites")

    def frac(a):
        return a[0] / a[1] if a[1] else np.nan
    vals = {}
    for s, h in hits.items():
        fm, ff = frac(h["male_specific"]), frac(h["female_specific"])
        vals[s] = np.nansum([fm, -ff]) if not (np.isnan(fm) and np.isnan(ff)) else np.nan
    return pd.Series(vals)


def kmer_values(path):
    k = pd.read_csv(path, sep="\t")
    return pd.Series((k["msk_fraction"].fillna(0) - k["fsk_fraction"].fillna(0)).values, index=k["sample"])


def to_score(values, known, default_mid, default_half, default_sign):
    """Scale raw values so known males ~ +1 and known females ~ -1."""
    m = values[[s for s, x in known.items() if x == "M" and s in values.index]].dropna()
    f = values[[s for s, x in known.items() if x == "F" and s in values.index]].dropna()
    if len(m) and len(f) and abs(m.mean() - f.mean()) > 1e-9:
        mid = (m.mean() + f.mean()) / 2
        half = abs(m.mean() - f.mean()) / 2
        sign = 1 if m.mean() > f.mean() else -1
        how = f"calibrated on {len(m)} M / {len(f)} F"
    elif default_mid is not None:
        mid, half, sign = default_mid, default_half, default_sign
        how = "theoretical expectation"
    else:
        return None, "not usable (no known-sex samples and no --system)"
    return ((values - mid) / half * sign).clip(-3, 3), how


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-sheet", required=True, help="TSV with columns sample, sex (M/F/U)")
    ap.add_argument("--system", choices=["XY", "ZW"], help="needed for depth evidence without known-sex samples")
    ap.add_argument("--regions", help="BED of known sex-linked regions (for depth)")
    ap.add_argument("--depth", help="per_sample_window_depth.tsv from depth_het_scan.py")
    ap.add_argument("--markers", help="sex_specific_sites.tsv (CHROM POS CLASS ...)")
    ap.add_argument("--genotypes", help="genotype table containing the samples to assign")
    ap.add_argument("--vcf-samples")
    ap.add_argument("--ploidy", type=int)
    ap.add_argument("--sample-ploidy")
    ap.add_argument("--kmer-fractions", help="TSV: sample, msk_fraction, fsk_fraction")
    ap.add_argument("--call-threshold", type=float, default=0.5)
    ap.add_argument("--outprefix", default="sex_assignment")
    args = ap.parse_args()

    sheet = pd.read_csv(args.sample_sheet, sep="\t", dtype=str)
    known = {s: x for s, x in zip(sheet["sample"], sheet["sex"]) if x in ("M", "F")}
    samples = list(sheet["sample"])

    evidence, notes = {}, []
    if args.depth and args.regions:
        sign = None if args.system is None else (-1 if args.system == "XY" else 1)
        s, how = to_score(depth_values(args.depth, args.regions), known,
                          0.75 if sign else None, 0.25, sign)
        notes.append(f"depth: {how}")
        if s is not None:
            evidence["depth"] = s
    if args.markers and args.genotypes:
        s, how = to_score(snp_values(args.markers, args.genotypes, args.vcf_samples,
                                     (args.ploidy, args.sample_ploidy)), known, 0.0, 0.5, 1)
        notes.append(f"snp: {how}")
        evidence["snp"] = s
    if args.kmer_fractions:
        s, how = to_score(kmer_values(args.kmer_fractions), known, 0.0, 0.5, 1)
        notes.append(f"kmer: {how}")
        evidence["kmer"] = s
    if not evidence:
        sys.exit("ERROR: no usable evidence. Give --depth + --regions (+ --system or known-sex samples), "
                 "--markers + --genotypes, and/or --kmer-fractions.")

    out = pd.DataFrame(index=samples)
    out.index.name = "sample"
    out["recorded_sex"] = [sheet.set_index("sample").loc[s, "sex"] for s in samples]
    for name, s in evidence.items():
        out[f"score_{name}"] = s.reindex(samples)
    score_cols = [c for c in out.columns if c.startswith("score_")]
    out["combined_score"] = out[score_cols].mean(axis=1)
    strong = out[score_cols].where(out[score_cols].abs() >= args.call_threshold)
    out["conflict"] = (strong.gt(0).any(axis=1) & strong.lt(0).any(axis=1))

    def call(r):
        c = r["combined_score"]
        if np.isnan(c) or r["conflict"] or abs(c) < args.call_threshold:
            return "uncertain"
        return "M" if c > 0 else "F"

    out["called_sex"] = out.apply(call, axis=1)
    # A typical male/female scores about +1/-1 on the calibrated scale.
    #   high   |score| >= 0.75 with >= 2 independent kinds of evidence
    #   medium |score| >= 0.75 from one kind of evidence
    #   low    call made, but the sample sits closer to the midpoint
    n_ev = out[score_cols].notna().sum(axis=1)
    typical = out["combined_score"].abs() >= 0.75
    out["confidence"] = np.where(out["called_sex"] == "uncertain", "-",
                                 np.where(typical & (n_ev >= 2), "high",
                                          np.where(typical, "medium", "low")))
    out["n_evidence"] = n_ev
    out["check"] = np.where(~out["recorded_sex"].isin(["M", "F"]), "assigned",
                            np.where(out["called_sex"] == "uncertain", "unclear",
                                     np.where(out["called_sex"] == out["recorded_sex"], "agrees", "MISMATCH")))
    out["marker_evidence_circular"] = out["recorded_sex"].isin(["M", "F"]) & \
        any(c in score_cols for c in ("score_snp", "score_kmer"))
    out.to_csv(f"{args.outprefix}.tsv", sep="\t", float_format="%.3f")

    with open(f"{args.outprefix}.notes.txt", "w") as f:
        f.write("\n".join(notes) + "\n")
    for n in notes:
        print(n)

    # Plot: scores per sample
    fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(samples)), 4))
    xs = np.arange(len(samples))
    for i, c in enumerate(score_cols):
        ax.scatter(xs + (i - 1) * 0.15, out[c], label=c.replace("score_", ""), s=30)
    ax.scatter(xs, out["combined_score"], marker="_", s=300, color="black", label="combined")
    ax.axhspan(-args.call_threshold, args.call_threshold, color="grey", alpha=0.15)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{s} ({r})" for s, r in zip(samples, out["recorded_sex"])], rotation=90)
    ax.set_ylabel("score (+ male-like, - female-like)")
    ax.set_title("Genomic sex assignment")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{args.outprefix}.png", dpi=200)

    print(out[["recorded_sex", "called_sex", "confidence", "combined_score", "check"]].to_string())
    mism = out[out["check"] == "MISMATCH"]
    if len(mism):
        print(f"\n{len(mism)} sample(s) disagree with their recorded sex: {', '.join(mism.index)}. "
              "Check sexing, labelling or contamination before using them.")
    print(f"Wrote: {args.outprefix}.tsv, .png, .notes.txt")


if __name__ == "__main__":
    main()
