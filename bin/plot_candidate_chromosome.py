#!/usr/bin/env python3
"""
Zoomed two-panel plot of male- and female-specific counts per window on one
chromosome (optionally restricted to a region and with a highlighted core).
"""

import argparse

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
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, default=None)
    ap.add_argument("--end", type=int, default=None)
    ap.add_argument("--highlight", type=int, nargs=2, metavar=("START", "END"), default=None,
                    help="Shade a core interval")
    ap.add_argument("--outprefix", default=None)
    args = ap.parse_args()

    outprefix = args.outprefix or f"{args.chrom}.sex_specific_zoom"
    panels = []
    for path in (args.male, args.female):
        df = read_windows(path)
        df = df[df["CHROM"] == args.chrom].copy()
        if args.start is not None:
            df = df[df["END"] >= args.start]
        if args.end is not None:
            df = df[df["START"] <= args.end]
        df["MID_MB"] = (df["START"] + df["END"]) / 2 / 1e6
        panels.append(df)

    fig, axes = plt.subplots(2, 1, figsize=(14, 6), sharex=True)
    for ax, df, label in zip(axes, panels, ("Male-specific\ncount", "Female-specific\ncount")):
        ax.scatter(df["MID_MB"], df["COUNT"], s=12)
        ax.plot(df["MID_MB"], df["COUNT"], linewidth=0.7)
        ax.set_ylabel(label)
        if args.highlight:
            ax.axvspan(args.highlight[0] / 1e6, args.highlight[1] / 1e6, alpha=0.15)
    axes[0].set_title(f"{args.chrom}: sex-specific SNP density")
    axes[1].set_xlabel("Position (Mb)")

    plt.tight_layout()
    plt.savefig(f"{outprefix}.png", dpi=300)
    plt.savefig(f"{outprefix}.pdf")
    print(f"Wrote: {outprefix}.png / .pdf")


if __name__ == "__main__":
    main()
