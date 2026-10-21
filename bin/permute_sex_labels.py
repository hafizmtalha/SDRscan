#!/usr/bin/env python3
"""
Sex-label permutation test for a candidate region, for any ploidy.

Genotypes stay fixed; sample sex labels are shuffled (keeping the number of
each sex) and the number of "strong" sites biased towards --focal-sex is
recomputed each time. The empirical p-value is the fraction of shuffles with
at least as many strong sites as the real labelling.

A strong site has ALT in >= --min-alt of the focal sex and in <= --max-opposite
of the other sex, with >= --min-called genotypes in each. Thresholds are counts
or fractions of the respective sex.
"""

import argparse
import sys

import numpy as np
import pandas as pd

from gtutils import GenotypeParser, load_ploidy, read_list, threshold


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--vcf-samples", default="vcf_samples.txt")
    ap.add_argument("--males", default="males.txt")
    ap.add_argument("--females", default="females.txt")
    ap.add_argument("--focal-sex", choices=["M", "F"], default="M",
                    help="M tests male-biased sites (XY), F female-biased (ZW)")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    pl = ap.add_mutually_exclusive_group(required=True)
    pl.add_argument("--ploidy", type=int)
    pl.add_argument("--sample-ploidy")
    ap.add_argument("--nperm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min-called", default="0.75")
    ap.add_argument("--min-alt", default="0.75")
    ap.add_argument("--max-opposite", default="1")
    ap.add_argument("--outprefix", default="perm")
    args = ap.parse_args()

    samples = read_list(args.vcf_samples)
    males = set(read_list(args.males))
    females = set(read_list(args.females))
    focal, other = (males, females) if args.focal_sex == "M" else (females, males)
    parser = GenotypeParser(load_ploidy(samples, args.ploidy, args.sample_ploidy))

    rows = []
    with open(args.genotypes) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != args.chrom:
                continue
            if args.start <= int(parts[1]) <= args.end:
                ds = [parser.dosage(gt, s) for s, gt in zip(samples, parts[2:])]
                rows.append([-1 if d is None else int(d > 0) for d in ds])
    A = np.array(rows, dtype=int)
    if A.size == 0:
        raise SystemExit("No sites in the requested region.")

    nf, no = len(focal), len(other)
    min_called_f, min_called_o = threshold(args.min_called, nf), threshold(args.min_called, no)
    min_alt = threshold(args.min_alt, nf)
    max_opp = threshold(args.max_opposite, no)

    def count_strong(f_idx, o_idx):
        fo, ot = A[:, f_idx], A[:, o_idx]
        ok = ((fo >= 0).sum(1) >= min_called_f) & ((ot >= 0).sum(1) >= min_called_o)
        strong = ((fo == 1).sum(1) >= min_alt) & ((ot == 1).sum(1) <= max_opp)
        return int((ok & strong).sum())

    labelled = [i for i, s in enumerate(samples) if s in focal or s in other]
    true_f = [i for i, s in enumerate(samples) if s in focal]
    true_o = [i for i, s in enumerate(samples) if s in other]
    obs = count_strong(true_f, true_o)

    rng = np.random.default_rng(args.seed)
    null = np.empty(args.nperm, dtype=int)
    for k in range(args.nperm):
        perm = rng.permutation(labelled)
        null[k] = count_strong(perm[:nf], perm[nf:nf + no])
    p = (np.sum(null >= obs) + 1) / (args.nperm + 1)

    pd.DataFrame({"null_strong_sites": null}).to_csv(f"{args.outprefix}.null.tsv", sep="\t", index=False)
    with open(f"{args.outprefix}.result.txt", "w") as f:
        f.write(f"region\t{args.chrom}:{args.start}-{args.end}\nfocal_sex\t{args.focal_sex}\n"
                f"sites\t{A.shape[0]}\nrule\tALT in >= {min_alt}/{nf} focal, <= {max_opp}/{no} other\n"
                f"observed\t{obs}\nnperm\t{args.nperm}\np_value\t{p:.6g}\n")
    parser.report(sys.stderr)
    print(f"Sites in region: {A.shape[0]}")
    print(f"Strong {args.focal_sex}-biased sites (ALT in >= {min_alt}/{nf}, other sex <= {max_opp}/{no}): {obs}")
    print(f"Empirical p-value (null >= observed): {p:.6g}")


if __name__ == "__main__":
    main()
