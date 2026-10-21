#!/usr/bin/env python3
"""
Summarise sex-biased ALT presence per site in a candidate region, to check
whether a window peak is driven by many repeated sex-biased SNPs rather than
a few noisy sites. Works for any ploidy (declare it with --ploidy or
--sample-ploidy).

Input genotypes: bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n' candidate.vcf.gz

Outputs:
  <prefix>.site_summary.tsv  per-site counts (and mean ALT dosage per sex)
  <prefix>.top_sites.tsv     top male-biased and female-biased sites
  <prefix>.summary.txt       short text summary

"Strong" thresholds are counts or fractions of each sex (0.75 = 75%).
"""

import argparse
import sys

import pandas as pd

from gtutils import GenotypeParser, load_ploidy, read_list, threshold


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf-samples", default="vcf_samples.txt")
    ap.add_argument("--males", default="males.txt")
    ap.add_argument("--females", default="females.txt")
    ap.add_argument("--genotypes", default="candidate_region.genotypes.tsv")
    ap.add_argument("--prefix", default="candidate_region")
    pl = ap.add_mutually_exclusive_group(required=True)
    pl.add_argument("--ploidy", type=int)
    pl.add_argument("--sample-ploidy")
    ap.add_argument("--min-called-per-sex", default="0.5", help="count or fraction of each sex")
    ap.add_argument("--top-n", type=int, default=50)
    ap.add_argument("--strong-alt", default="0.75",
                    help="'strong' site: ALT in at least this many of the biased sex (count or fraction)")
    ap.add_argument("--strong-max-opposite", default="1",
                    help="... and at most this many of the other sex (count or fraction)")
    args = ap.parse_args()

    vcf_samples = read_list(args.vcf_samples)
    males = set(read_list(args.males))
    females = set(read_list(args.females))
    if males & females:
        raise SystemExit(f"ERROR: samples in both sex lists: {sorted(males & females)}")
    missing = [s for s in males | females if s not in vcf_samples]
    if missing:
        raise SystemExit(f"ERROR: samples not in {args.vcf_samples}: {missing}")
    parser = GenotypeParser(load_ploidy(vcf_samples, args.ploidy, args.sample_ploidy))

    nm, nf = len(males), len(females)
    min_called_m = threshold(args.min_called_per_sex, nm)
    min_called_f = threshold(args.min_called_per_sex, nf)

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
            md = fd = 0
            for s, gt in zip(vcf_samples, gts):
                if s not in males and s not in females:
                    continue
                d = parser.dosage(gt, s)
                if d is None:
                    continue
                if s in males:
                    mc += 1; ma += d > 0; md += d
                else:
                    fc += 1; fa += d > 0; fd += d
            if mc < min_called_m or fc < min_called_f:
                continue
            rows.append([parts[0], pos, mc, fc, ma, fa, ma - fa,
                         round(md / mc, 3), round(fd / fc, 3)])

    df = pd.DataFrame(rows, columns=["CHROM", "POS", "male_called", "female_called", "male_alt",
                                     "female_alt", "alt_diff", "male_mean_dosage", "female_mean_dosage"])
    df.to_csv(f"{args.prefix}.site_summary.tsv", sep="\t", index=False)

    top_m = df.sort_values("alt_diff", ascending=False).head(args.top_n).assign(BIAS="male_biased")
    top_f = df.sort_values("alt_diff", ascending=True).head(args.top_n).assign(BIAS="female_biased")
    top = pd.concat([top_m, top_f], ignore_index=True)
    top = top[["BIAS"] + [c for c in df.columns]]
    top.to_csv(f"{args.prefix}.top_sites.tsv", sep="\t", index=False)

    sm, sf = threshold(args.strong_alt, nm), threshold(args.strong_alt, nf)
    om, of = threshold(args.strong_max_opposite, nm), threshold(args.strong_max_opposite, nf)
    strong_m = ((df["male_alt"] >= sm) & (df["female_alt"] <= of)).sum()
    strong_f = ((df["female_alt"] >= sf) & (df["male_alt"] <= om)).sum()

    summary = (
        f"Input genotypes: {args.genotypes}\n"
        f"Total lines read: {total}\n"
        f"Sites kept (called >= {min_called_m} males and >= {min_called_f} females): {len(df)}\n"
        f"Strong male-biased   (male_alt>={sm} & female_alt<={of}): {strong_m}\n"
        f"Strong female-biased (female_alt>={sf} & male_alt<={om}): {strong_f}\n"
    )
    with open(f"{args.prefix}.summary.txt", "w") as f:
        f.write(summary)
    parser.report(sys.stderr)
    print(summary, end="")
    print(f"Wrote: {args.prefix}.site_summary.tsv, {args.prefix}.top_sites.tsv, {args.prefix}.summary.txt")


if __name__ == "__main__":
    main()
