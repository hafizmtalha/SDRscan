#!/usr/bin/env python3
"""
Summarise sex-biased ALT presence per site in a candidate region, to check
whether a window peak is driven by many repeated sex-biased SNPs rather than
a few noisy sites.

Input genotypes: bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n' candidate.vcf.gz

Outputs:
  <prefix>.site_summary.tsv  per-site counts
  <prefix>.top_sites.tsv     top male-biased and female-biased sites
  <prefix>.summary.txt       short text summary
"""

import argparse
import re

import pandas as pd

MISSING = {".", "./.", ".|.", "./././.", ".|.|.|.", "././././.", ".|.|.|.|."}


def read_list(path):
    with open(path) as f:
        return [x.strip() for x in f if x.strip()]


def gt_called(gt):
    return gt is not None and gt.strip() not in MISSING


def gt_has_alt(gt):
    if not gt_called(gt):
        return False
    return any(a not in {"", ".", "0"} for a in re.split(r"[\/|]", gt.strip()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf-samples", default="vcf_samples.txt")
    ap.add_argument("--males", default="males.txt")
    ap.add_argument("--females", default="females.txt")
    ap.add_argument("--genotypes", default="candidate_region.genotypes.tsv")
    ap.add_argument("--prefix", default="candidate_region")
    ap.add_argument("--min-called-per-sex", type=int, default=3)
    ap.add_argument("--top-n", type=int, default=50)
    ap.add_argument("--strong-male-alt", type=int, default=4,
                    help="'Strong male-biased' = male_alt >= this and female_alt <= --strong-max-opposite")
    ap.add_argument("--strong-female-alt", type=int, default=3)
    ap.add_argument("--strong-max-opposite", type=int, default=1)
    args = ap.parse_args()

    vcf_samples = read_list(args.vcf_samples)
    males = set(read_list(args.males))
    females = set(read_list(args.females))
    if males & females:
        raise SystemExit(f"ERROR: samples in both sex lists: {sorted(males & females)}")
    missing = [s for s in males | females if s not in vcf_samples]
    if missing:
        raise SystemExit(f"ERROR: samples not in {args.vcf_samples}: {missing}")

    rows = []
    total = 0
    with open(args.genotypes) as f:
        for line in f:
            total += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            gts = parts[2:]
            if len(gts) != len(vcf_samples):
                raise SystemExit(f"ERROR: {len(gts)} genotype columns != {len(vcf_samples)} samples at {parts[0]}:{pos}")
            mc = fc = ma = fa = 0
            for s, gt in zip(vcf_samples, gts):
                called, alt = gt_called(gt), gt_has_alt(gt)
                if s in males:
                    mc += called
                    ma += alt
                elif s in females:
                    fc += called
                    fa += alt
            if mc < args.min_called_per_sex or fc < args.min_called_per_sex:
                continue
            rows.append([parts[0], pos, mc, fc, ma, fa, ma - fa])

    df = pd.DataFrame(rows, columns=["CHROM", "POS", "male_called", "female_called",
                                     "male_alt", "female_alt", "alt_diff"])
    df.to_csv(f"{args.prefix}.site_summary.tsv", sep="\t", index=False)

    top_m = df.sort_values("alt_diff", ascending=False).head(args.top_n).assign(BIAS="male_biased")
    top_f = df.sort_values("alt_diff", ascending=True).head(args.top_n).assign(BIAS="female_biased")
    top = pd.concat([top_m, top_f], ignore_index=True)
    top = top[["BIAS", "CHROM", "POS", "male_called", "female_called", "male_alt", "female_alt", "alt_diff"]]
    top.to_csv(f"{args.prefix}.top_sites.tsv", sep="\t", index=False)

    strong_m = ((df["male_alt"] >= args.strong_male_alt) & (df["female_alt"] <= args.strong_max_opposite)).sum()
    strong_f = ((df["female_alt"] >= args.strong_female_alt) & (df["male_alt"] <= args.strong_max_opposite)).sum()

    summary = (
        f"Input genotypes: {args.genotypes}\n"
        f"Total lines read: {total}\n"
        f"Sites kept (min_called_per_sex={args.min_called_per_sex}): {len(df)}\n"
        f"Strong male-biased   (male_alt>={args.strong_male_alt} & female_alt<={args.strong_max_opposite}): {strong_m}\n"
        f"Strong female-biased (female_alt>={args.strong_female_alt} & male_alt<={args.strong_max_opposite}): {strong_f}\n"
    )
    with open(f"{args.prefix}.summary.txt", "w") as f:
        f.write(summary)
    print(summary, end="")
    print(f"Wrote: {args.prefix}.site_summary.tsv, {args.prefix}.top_sites.tsv, {args.prefix}.summary.txt")


if __name__ == "__main__":
    main()
