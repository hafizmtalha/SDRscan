#!/usr/bin/env python3
"""
Two-panel genome-wide plot of male- and female-specific SNP counts per window,
plus a table of the top-N windows for each sex.

Inputs are bedtools intersect -c outputs (CHROM START END COUNT).
"""

import argparse
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def read_windows(path):
    return pd.read_csv(path, sep="\t", header=None, names=["CHROM", "START", "END", "COUNT"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--male", default="male_specific_windows.bed")
    ap.add_argument("--female", default="female_specific_windows.bed")
    ap.add_argument("--outprefix", default="sex_specific_snp_windows")
    ap.add_argument("--top-windows", default="top_candidate_windows.tsv")
    ap.add_argument("--top-n", type=int, default=50)
    ap.add_argument("--chrom-regex", default=None,
                    help="Keep only chromosomes matching this regex (drops unplaced scaffolds)")
    ap.add_argument("--label-regex", default=r"(\d+)$",
                    help="Regex whose first group becomes the axis label (default: trailing number)")
    ap.add_argument("--title", default="Genome-wide density of sex-specific SNPs")
    args = ap.parse_args()

    male = read_windows(args.male)
    female = read_windows(args.female)
    if args.chrom_regex:
        male = male[male["CHROM"].str.contains(args.chrom_regex, regex=True)]
        female = female[female["CHROM"].str.contains(args.chrom_regex, regex=True)]

    chrom_order = list(dict.fromkeys(male["CHROM"].tolist()))
    offsets, ticks, labels, lengths = {}, [], [], {}
    offset = 0
    for chrom in chrom_order:
        chrom_len = max(male.loc[male["CHROM"] == chrom, "END"].max(),
                        female.loc[female["CHROM"] == chrom, "END"].max())
        offsets[chrom] = offset
        lengths[chrom] = chrom_len
        ticks.append(offset + chrom_len / 2)
        m = re.search(args.label_regex, chrom)
        labels.append(f"Chr{int(m.group(1))}" if m and m.group(1).isdigit() else chrom)
        offset += chrom_len

    for df in (male, female):
        df["GENOME_POS"] = df["START"] + df["CHROM"].map(offsets)

    male_top = male.sort_values("COUNT", ascending=False).head(args.top_n).assign(SEX="male_specific")
    female_top = female.sort_values("COUNT", ascending=False).head(args.top_n).assign(SEX="female_specific")
    top = pd.concat([male_top, female_top], ignore_index=True)[["SEX", "CHROM", "START", "END", "COUNT"]]
    top.to_csv(args.top_windows, sep="\t", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)
    for i, chrom in enumerate(chrom_order):
        if i % 2 == 0:
            for ax in axes:
                ax.axvspan(offsets[chrom], offsets[chrom] + lengths[chrom], alpha=0.08)

    for ax, df, ylabel in ((axes[0], male, "Male-specific\nSNP count"),
                           (axes[1], female, "Female-specific\nSNP count")):
        ax.scatter(df["GENOME_POS"], df["COUNT"], s=8)
        ax.plot(df["GENOME_POS"], df["COUNT"], linewidth=0.5)
        ax.set_ylabel(ylabel)
    axes[0].set_title(args.title)
    axes[1].set_xlabel("Chromosome")
    axes[1].set_xticks(ticks)
    axes[1].set_xticklabels(labels, rotation=90)

    plt.tight_layout()
    plt.savefig(f"{args.outprefix}.png", dpi=300)
    plt.savefig(f"{args.outprefix}.pdf")
    print(f"Wrote: {args.outprefix}.png / .pdf")
    print(f"Wrote: {args.top_windows}")


if __name__ == "__main__":
    main()
