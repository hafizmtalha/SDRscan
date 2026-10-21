#!/usr/bin/env python3
"""
Sex-label permutation test for a candidate region.

Genotypes stay fixed; sample sex labels are shuffled (keeping the number of
males and females) and the number of "strong male-biased" sites in the region
is recomputed each time. The empirical p-value is the fraction of shuffles
with at least as many strong sites as the real labelling.
"""

import argparse
import re

import numpy as np
import pandas as pd

MISSING = {".", "./.", ".|.", "./././.", ".|.|.|.", "././././.", ".|.|.|.|."}


def read_list(path):
    with open(path) as f:
        return [x.strip() for x in f if x.strip()]


def alt_presence(gt):
    if gt is None or gt.strip() in MISSING:
        return -1
    return int(any(a not in {"", ".", "0"} for a in re.split(r"[\/|]", gt.strip())))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--vcf-samples", default="vcf_samples.txt")
    ap.add_argument("--males", default="males.txt")
    ap.add_argument("--females", default="females.txt")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--nperm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min-called-male", type=int, default=4)
    ap.add_argument("--min-called-female", type=int, default=3)
    ap.add_argument("--min-male-alt", type=int, default=4)
    ap.add_argument("--max-female-alt", type=int, default=1)
    ap.add_argument("--outprefix", default="perm")
    args = ap.parse_args()

    samples = read_list(args.vcf_samples)
    males = set(read_list(args.males))
    females = set(read_list(args.females))

    rows = []
    with open(args.genotypes) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != args.chrom:
                continue
            if args.start <= int(parts[1]) <= args.end:
                rows.append([alt_presence(gt) for gt in parts[2:]])
    A = np.array(rows, dtype=int)
    if A.size == 0:
        raise SystemExit("No sites in the requested region.")

    def count_strong(male_idx, female_idx):
        m, fem = A[:, male_idx], A[:, female_idx]
        ok = ((m >= 0).sum(1) >= args.min_called_male) & ((fem >= 0).sum(1) >= args.min_called_female)
        strong = ((m == 1).sum(1) >= args.min_male_alt) & ((fem == 1).sum(1) <= args.max_female_alt)
        return int((ok & strong).sum())

    labelled = [i for i, s in enumerate(samples) if s in males or s in females]
    true_m = [i for i, s in enumerate(samples) if s in males]
    true_f = [i for i, s in enumerate(samples) if s in females]
    obs = count_strong(true_m, true_f)

    rng = np.random.default_rng(args.seed)
    null = np.empty(args.nperm, dtype=int)
    for k in range(args.nperm):
        perm = rng.permutation(labelled)
        null[k] = count_strong(perm[:len(true_m)], perm[len(true_m):len(true_m) + len(true_f)])
    p = (np.sum(null >= obs) + 1) / (args.nperm + 1)

    pd.DataFrame({"null_strong_male_sites": null}).to_csv(f"{args.outprefix}.null.tsv", sep="\t", index=False)
    with open(f"{args.outprefix}.result.txt", "w") as f:
        f.write(f"region\t{args.chrom}:{args.start}-{args.end}\nsites\t{A.shape[0]}\n"
                f"observed\t{obs}\nnperm\t{args.nperm}\np_value\t{p:.6g}\n")
    print(f"Sites in region: {A.shape[0]}")
    print(f"Observed strong male-biased sites: {obs}")
    print(f"Empirical p-value (null >= observed): {p:.6g}")


if __name__ == "__main__":
    main()
