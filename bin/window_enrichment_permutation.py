#!/usr/bin/env python3
"""
Is one chromosome enriched for sex-specific signal compared with the rest of
the genome?

Takes a window-count BED (CHROM START END COUNT, from SNP or k-mer density)
and compares the mean count on --chrom with the means of --nperm random sets
of the same number of windows drawn from the whole genome. Also writes a
per-chromosome summary and a histogram of the null distribution.
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--windows", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--nperm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--chrom-regex", default=None, help="Keep only chromosomes matching this regex")
    ap.add_argument("--outprefix", default="window_enrichment")
    args = ap.parse_args()

    df = pd.read_csv(args.windows, sep="\t", header=None, names=["CHROM", "START", "END", "COUNT"])
    if args.chrom_regex:
        df = df[df["CHROM"].str.contains(args.chrom_regex, regex=True)]
    target = df[df["CHROM"] == args.chrom]["COUNT"].to_numpy()
    if target.size == 0:
        raise SystemExit(f"No windows for {args.chrom}")
    allc = df["COUNT"].to_numpy()

    obs = target.mean()
    rng = np.random.default_rng(args.seed)
    null = np.array([rng.choice(allc, target.size, replace=False).mean() for _ in range(args.nperm)])
    p = (np.sum(null >= obs) + 1) / (args.nperm + 1)

    per_chrom = (df.groupby("CHROM")["COUNT"]
                 .agg(n_windows="size", total="sum", mean="mean", max="max")
                 .sort_values("mean", ascending=False))
    per_chrom.to_csv(f"{args.outprefix}.per_chrom.tsv", sep="\t")
    with open(f"{args.outprefix}.result.txt", "w") as f:
        f.write(f"chrom\t{args.chrom}\nwindows\t{target.size}\nobserved_mean\t{obs:.6g}\n"
                f"genome_mean\t{allc.mean():.6g}\nnperm\t{args.nperm}\np_value\t{p:.6g}\n")

    plt.figure(figsize=(7, 4))
    plt.hist(null, bins=50)
    plt.axvline(obs, color="red", linestyle="--", label=f"{args.chrom} (p={p:.3g})")
    plt.xlabel("Mean count per window (random window sets)")
    plt.ylabel("Permutations")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{args.outprefix}.null.png", dpi=300)

    print(f"{args.chrom}: mean {obs:.4g} vs genome {allc.mean():.4g}, p = {p:.4g}")
    print(f"Wrote: {args.outprefix}.result.txt, .per_chrom.tsv, .null.png")


if __name__ == "__main__":
    main()
