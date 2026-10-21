#!/usr/bin/env python3
"""
Validation layer: count how many strict (high-confidence) sex-specific sites
fall inside the candidate intervals found with relaxed thresholds.

Requires bedtools on PATH.
"""

import argparse
import os
import subprocess
import tempfile

import pandas as pd


def count_overlaps(bed_a, bed_b):
    out = subprocess.check_output(["bedtools", "intersect", "-a", bed_a, "-b", bed_b, "-c"], text=True)
    return [int(line.split("\t")[-1]) for line in out.strip().splitlines()]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--intervals", default="merged_candidate_intervals.tsv")
    ap.add_argument("--strict-male", default="male_specific_sites.strict.bed")
    ap.add_argument("--strict-female", default="female_specific_sites.strict.bed")
    ap.add_argument("--out", default="strict_overlap_counts.tsv")
    ap.add_argument("--top-n", type=int, default=20)
    ap.add_argument("--pad", type=int, default=0)
    args = ap.parse_args()

    iv = (pd.read_csv(args.intervals, sep="\t")
          .sort_values(["MAX_COUNT", "N_WINDOWS", "MEAN_COUNT"], ascending=False)
          .head(args.top_n))
    bed = pd.DataFrame({
        "CHROM": iv["CHROM"],
        "START": (iv["START"] - args.pad).clip(lower=0).astype(int),
        "END": (iv["END"] + args.pad).astype(int),
        "SEX": iv["SEX"],
        "N_WINDOWS": iv["N_WINDOWS"],
        "MAX_COUNT": iv["MAX_COUNT"],
    })

    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tmp:
        bed.to_csv(tmp.name, sep="\t", header=False, index=False)
    try:
        bed["STRICT_MALE_SITES"] = count_overlaps(tmp.name, args.strict_male)
        bed["STRICT_FEMALE_SITES"] = count_overlaps(tmp.name, args.strict_female)
    finally:
        os.remove(tmp.name)

    bed.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote: {args.out}")
    print(bed.sort_values(["STRICT_MALE_SITES", "STRICT_FEMALE_SITES"], ascending=False)
          .head(15).to_string(index=False))


if __name__ == "__main__":
    main()
