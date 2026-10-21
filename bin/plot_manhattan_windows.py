#!/usr/bin/env python3
"""
Manhattan-style plot of window counts for several tracks on one panel, e.g.
male- and female-specific SNPs on two references, or sex-specific k-mer
density mapped to four haplotypes.

Each --track is LABEL=PATH[:COLOR[:MARKER]], where PATH is a window-count BED
(CHROM START END COUNT). Chromosomes are laid end to end in the order they
first appear in the first track; tracks on different references are aligned
by chromosome number (taken from --label-regex), so chr7 of every reference
falls in the same slot.

Example:
  plot_manhattan_windows.py \\
    --track "REF1 male-specific SNPs=fh1/male_windows.bed:green:o" \\
    --track "REF1 female-specific SNPs=fh1/female_windows.bed:purple:o" \\
    --track "REF2 male-specific SNPs=mh1/male_windows.bed:blue:^" \\
    --track "REF2 female-specific SNPs=mh1/female_windows.bed:orange:^" \\
    --chrom-regex 'chr_?\\d+$' --title "Genome-wide sex-specific SNP density" \\
    --outprefix snp_manhattan
"""

import argparse
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def parse_track(spec):
    label, rest = spec.split("=", 1)
    parts = rest.split(":")
    return label, parts[0], (parts[1] if len(parts) > 1 else None), (parts[2] if len(parts) > 2 else "o")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--track", action="append", required=True)
    ap.add_argument("--chrom-regex", default=None, help="Keep only chromosomes matching this regex")
    ap.add_argument("--label-regex", default=r"(\d+)$", help="Regex capturing the chromosome number")
    ap.add_argument("--min-count", type=int, default=1, help="Hide windows below this count")
    ap.add_argument("--title", default="Genome-wide sex-specific density")
    ap.add_argument("--ylabel", default="Count per window")
    ap.add_argument("--outprefix", default="manhattan")
    args = ap.parse_args()

    tracks = []
    for spec in args.track:
        label, path, color, marker = parse_track(spec)
        df = pd.read_csv(path, sep="\t", header=None, names=["CHROM", "START", "END", "COUNT"])
        if args.chrom_regex:
            df = df[df["CHROM"].str.contains(args.chrom_regex, regex=True)]
        df["CHRNUM"] = df["CHROM"].str.extract(args.label_regex)[0].astype(float)
        df = df.dropna(subset=["CHRNUM"])
        df["CHRNUM"] = df["CHRNUM"].astype(int)
        tracks.append((label, df, color, marker))

    lengths = pd.concat([t[1] for t in tracks]).groupby("CHRNUM")["END"].max().sort_index()
    offsets = lengths.cumsum() - lengths

    fig, ax = plt.subplots(figsize=(18, 6))
    for i, (num, length) in enumerate(lengths.items()):
        if i % 2 == 0:
            ax.axvspan(offsets[num], offsets[num] + length, color="grey", alpha=0.08)
    for label, df, color, marker in tracks:
        df = df[df["COUNT"] >= args.min_count]
        x = df["START"] + df["CHRNUM"].map(offsets)
        ax.scatter(x, df["COUNT"], s=18, c=color, marker=marker, label=label, alpha=0.8)

    ax.set_xticks(offsets + lengths / 2)
    ax.set_xticklabels([f"Chr{n}" for n in lengths.index])
    ax.set_xlim(0, (offsets + lengths).max())
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(args.ylabel)
    ax.set_title(args.title)
    ax.legend(loc="upper right", fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{args.outprefix}.png", dpi=300)
    plt.savefig(f"{args.outprefix}.pdf")
    print(f"Wrote: {args.outprefix}.png / .pdf")


if __name__ == "__main__":
    main()
