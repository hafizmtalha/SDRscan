#!/usr/bin/env python3
"""
Read-depth and heterozygosity scan between the sexes, for any ploidy and any
number of samples (including a 1 male + 1 female pilot).

Depth: per sample, mean depth per window (samtools bedcov) divided by that
sample's genome-wide median window depth, so ~1 means "autosome-like". For
each window the mean normalised depth per sex and log2(male / female) are
reported. A region present in one copy in the heterogametic sex (X in XY
males, Z in ZW females) gives log2 ~ -1 (XY) or ~ +1 (ZW); sequence missing
from one sex gives a much stronger shift.

Heterozygosity: per sample, the share of called sites in the window whose
genotype carries both REF and ALT (0 < ALT copies < ploidy). Hemizygous
regions lose heterozygosity in the heterogametic sex; a diverged but still
paired region (young Y/W) gains it.

"Separated" windows are those where every male differs from every female in
the same direction (all males below all females, or above). With one sample
per sex every window is trivially "separated", so in that case candidates
rest on effect size alone and outputs are labelled PILOT.

Outputs (<outprefix>.*):
  windows.tsv                per-window depth/het per sex, log2 ratio, flags
  per_sample_window_depth.tsv  normalised depth per sample (input for sexing)
  candidates.tsv             merged runs of windows passing the thresholds
  per_chrom.tsv              per-chromosome medians (whole-chromosome signals)
  genome.png/.pdf            genome-wide log2 depth ratio and het difference
"""

import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gtutils import GenotypeParser, load_ploidy, read_list

EPS = 1e-3


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bedcov", required=True,
                    help="samtools bedcov output: CHROM START END sum_depth(sample1) sum_depth(sample2) ...")
    ap.add_argument("--bam-samples", required=True, help="sample names in the order of the bedcov columns")
    ap.add_argument("--sample-sheet", required=True, help="TSV with columns sample, sex (M/F/U)")
    ap.add_argument("--genotypes", help="genotype table for heterozygosity (optional)")
    ap.add_argument("--vcf-samples")
    ap.add_argument("--sample-ploidy")
    ap.add_argument("--log2-threshold", type=float, default=0.6,
                    help="min |log2(M/F)| for a depth candidate window (0.6 ~ 1.5-fold)")
    ap.add_argument("--het-threshold", type=float, default=0.1,
                    help="min |het(M) - het(F)| for a heterozygosity candidate window")
    ap.add_argument("--min-depth", type=float, default=0.05,
                    help="ignore windows whose median normalised depth across samples is below this")
    ap.add_argument("--strong-windows", type=int, default=3,
                    help="candidates spanning at least this many adjacent windows are marked 'strong', "
                         "shorter ones 'weak' (column support)")
    ap.add_argument("--outprefix", default="depth_het")
    args = ap.parse_args()

    sheet = pd.read_csv(args.sample_sheet, sep="\t", dtype=str)
    sex_of = dict(zip(sheet["sample"], sheet["sex"]))
    bam_samples = read_list(args.bam_samples)

    cov = pd.read_csv(args.bedcov, sep="\t", header=None)
    cov.columns = ["CHROM", "START", "END"] + bam_samples
    length = (cov["END"] - cov["START"]).clip(lower=1)
    depth = cov[bam_samples].div(length, axis=0)
    norm = depth / depth.median(axis=0).replace(0, np.nan)

    males = [s for s in bam_samples if sex_of.get(s) == "M"]
    females = [s for s in bam_samples if sex_of.get(s) == "F"]
    if not males or not females:
        sys.exit(f"ERROR: need at least one male and one female (have {len(males)} M, {len(females)} F)")
    pilot = min(len(males), len(females)) < 2
    tag = f"{len(males)} M vs {len(females)} F" + (" - PILOT, no replication" if pilot else "")
    print(f"Depth/heterozygosity scan: {tag}")

    win = cov[["CHROM", "START", "END"]].copy()
    win["depth_M"] = norm[males].mean(axis=1)
    win["depth_F"] = norm[females].mean(axis=1)
    win["log2_MF"] = np.log2((win["depth_M"] + EPS) / (win["depth_F"] + EPS))
    win["depth_separated"] = ((norm[males].max(axis=1) < norm[females].min(axis=1)) |
                              (norm[males].min(axis=1) > norm[females].max(axis=1)))
    usable = norm.median(axis=1) >= args.min_depth

    # Heterozygosity per sample per window
    have_het = bool(args.genotypes)
    if have_het:
        vcf_samples = read_list(args.vcf_samples)
        parser = GenotypeParser(load_ploidy(vcf_samples, None, args.sample_ploidy))
        ploidy = parser.ploidy
        idx = {(c, s): i for i, (c, s) in enumerate(zip(win["CHROM"], win["START"]))}
        size = int(length.median())
        het = np.zeros((len(win), len(vcf_samples)))
        called = np.zeros((len(win), len(vcf_samples)))
        with open(args.genotypes) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                pos = int(parts[1])
                i = idx.get((parts[0], (pos - 1) // size * size))
                if i is None:
                    continue
                for j, (s, gt) in enumerate(zip(vcf_samples, parts[2:])):
                    d = parser.dosage(gt, s)
                    if d is None:
                        continue
                    called[i, j] += 1
                    het[i, j] += 0 < d < ploidy[s]
        parser.report(sys.stderr)
        rate = pd.DataFrame(np.where(called > 0, het / np.maximum(called, 1), np.nan), columns=vcf_samples)
        hm = [s for s in males if s in vcf_samples]
        hf = [s for s in females if s in vcf_samples]
        win["het_M"] = rate[hm].mean(axis=1).values
        win["het_F"] = rate[hf].mean(axis=1).values
        win["het_diff_MF"] = win["het_M"] - win["het_F"]
        win["het_separated"] = ((rate[hm].max(axis=1) < rate[hf].min(axis=1)) |
                                (rate[hm].min(axis=1) > rate[hf].max(axis=1))).values

    # Candidate windows
    depth_hit = usable & (win["log2_MF"].abs() >= args.log2_threshold)
    het_hit = (win["het_diff_MF"].abs() >= args.het_threshold) if have_het else pd.Series(False, index=win.index)
    if not pilot:
        depth_hit &= win["depth_separated"]
        if have_het:
            het_hit &= win["het_separated"]
    win["depth_candidate"] = depth_hit
    win["het_candidate"] = het_hit
    win.to_csv(f"{args.outprefix}.windows.tsv", sep="\t", index=False, float_format="%.4g")

    out_norm = pd.concat([cov[["CHROM", "START", "END"]], norm], axis=1)
    out_norm.to_csv(f"{args.outprefix}.per_sample_window_depth.tsv", sep="\t", index=False, float_format="%.4g")

    # Merge consecutive candidate windows with the same direction
    rows = []
    for kind, col, hit in (("depth", "log2_MF", depth_hit), ("het", "het_diff_MF", het_hit)):
        if col not in win:
            continue
        sub = win[hit].copy()
        sub["dir"] = np.sign(sub[col])
        for chrom, g in sub.groupby("CHROM", sort=False):
            g = g.sort_values("START")
            run = []
            for _, r in g.iterrows():
                if run and (r["START"] > run[-1]["END"] or r["dir"] != run[-1]["dir"]):
                    rows.append(_summarise(kind, col, run))
                    run = []
                run.append(r)
            if run:
                rows.append(_summarise(kind, col, run))
    cand = pd.DataFrame(rows, columns=["evidence", "CHROM", "START", "END", "n_windows", "mean_effect", "interpretation"])
    cand = cand.sort_values(["n_windows", "mean_effect"], ascending=[False, False], key=lambda s: s.abs())
    # Single windows cross the thresholds by chance, especially in pilots;
    # a run of several adjacent windows is the pattern a sex-linked region gives
    cand["support"] = np.where(cand["n_windows"] >= args.strong_windows, "strong", "weak")
    cand.to_csv(f"{args.outprefix}.candidates.tsv", sep="\t", index=False, float_format="%.4g")

    agg = {"log2_MF": "median"}
    if have_het:
        agg["het_diff_MF"] = "median"
    per_chrom = win[usable].groupby("CHROM", sort=False).agg(agg)
    per_chrom["n_depth_candidates"] = win[depth_hit].groupby("CHROM").size()
    per_chrom = per_chrom.fillna({"n_depth_candidates": 0})
    per_chrom.to_csv(f"{args.outprefix}.per_chrom.tsv", sep="\t", float_format="%.4g")

    _plot(win, usable, have_het, f"{args.outprefix}.genome", tag)
    print(f"Candidate depth windows: {int(depth_hit.sum())}, het windows: {int(het_hit.sum())}")
    print(f"Merged candidates: {len(cand)} -> {args.outprefix}.candidates.tsv")
    if pilot:
        print("PILOT: with one sample per sex these are differences between two individuals; "
              "confirm any candidate with more samples.")


def _summarise(kind, col, run):
    eff = float(np.mean([r[col] for r in run]))
    if kind == "depth":
        meaning = ("lower depth in males (X-like if XY, or male-missing sequence)" if eff < 0
                   else "lower depth in females (Z-like if ZW, or female-missing sequence)")
    else:
        meaning = "higher heterozygosity in males" if eff > 0 else "higher heterozygosity in females"
    return [kind, run[0]["CHROM"], int(run[0]["START"]), int(run[-1]["END"]), len(run), eff, meaning]


def _plot(win, usable, have_het, prefix, tag):
    w = win[usable].copy()
    chroms = list(dict.fromkeys(w["CHROM"]))
    lengths = w.groupby("CHROM", sort=False)["END"].max()
    offsets = (lengths.cumsum() - lengths).to_dict()
    x = w["START"] + w["CHROM"].map(offsets)
    colours = np.where(w["CHROM"].map({c: i % 2 for i, c in enumerate(chroms)}) == 0, "#3b6ea5", "#9aa5b1")
    panels = [("log2_MF", "log2 depth (M/F)")] + ([("het_diff_MF", "het(M) - het(F)")] if have_het else [])
    fig, axes = plt.subplots(len(panels), 1, figsize=(18, 4 * len(panels)), sharex=True, squeeze=False)
    for ax, (col, label) in zip(axes[:, 0], panels):
        ax.scatter(x, w[col], s=4, c=colours)
        ax.axhline(0, color="black", linewidth=0.6)
        ax.set_ylabel(label)
    axes[0, 0].set_title(f"Depth and heterozygosity between sexes ({tag})")
    axes[-1, 0].set_xticks([offsets[c] + lengths[c] / 2 for c in chroms])
    axes[-1, 0].set_xticklabels(chroms, rotation=90, fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{prefix}.png", dpi=200)
    plt.savefig(f"{prefix}.pdf")


if __name__ == "__main__":
    main()
