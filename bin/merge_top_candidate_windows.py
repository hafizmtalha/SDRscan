#!/usr/bin/env python3
"""
Merge overlapping / nearby top windows into candidate intervals.

Input:  top_candidate_windows.tsv (SEX CHROM START END COUNT)
Output: merged_candidate_intervals.tsv
        (SEX CHROM START END N_WINDOWS MAX_COUNT MEAN_COUNT), strongest first.
"""

import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="top_candidate_windows.tsv")
    ap.add_argument("--output", default="merged_candidate_intervals.tsv")
    ap.add_argument("--merge-gap", type=int, default=50000,
                    help="Merge windows that overlap or lie within this many bp (default: one window step)")
    ap.add_argument("--min-count", type=int, default=2,
                    help="Ignore windows with fewer sex-specific SNPs than this")
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    df = df[df["COUNT"] >= args.min_count].sort_values(["SEX", "CHROM", "START", "END"])

    merged = []
    for (sex, chrom), group in df.groupby(["SEX", "CHROM"], sort=False):
        cur_start = cur_end = None
        counts = []
        for _, r in group.iterrows():
            start, end, count = int(r["START"]), int(r["END"]), int(r["COUNT"])
            if cur_start is not None and start <= cur_end + args.merge_gap:
                cur_end = max(cur_end, end)
                counts.append(count)
                continue
            if cur_start is not None:
                merged.append([sex, chrom, cur_start, cur_end, len(counts), max(counts), sum(counts) / len(counts)])
            cur_start, cur_end, counts = start, end, [count]
        if cur_start is not None:
            merged.append([sex, chrom, cur_start, cur_end, len(counts), max(counts), sum(counts) / len(counts)])

    out = pd.DataFrame(merged, columns=["SEX", "CHROM", "START", "END", "N_WINDOWS", "MAX_COUNT", "MEAN_COUNT"])
    out = out.sort_values(["SEX", "MAX_COUNT", "N_WINDOWS", "MEAN_COUNT"],
                          ascending=[True, False, False, False]).reset_index(drop=True)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote: {args.output}")
    print(out.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
