#!/usr/bin/env python3
"""
Batch version of inspect_candidate_region.py: take the top-N merged candidate
intervals, pad them, subset the VCF with bcftools, export genotypes and
summarise sex-biased sites for each interval.

Requires bcftools and tabix on PATH. Outputs go to --outdir, with
candidate_interval_summary.tsv as the overview.
"""

import argparse
import os
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def run(cmd, **kw):
    subprocess.run(cmd, check=True, **kw)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", default="merged.biallelic.snps.vcf.gz")
    ap.add_argument("--vcf-samples", default="vcf_samples.txt")
    ap.add_argument("--males", default="males.txt")
    ap.add_argument("--females", default="females.txt")
    ap.add_argument("--intervals", default="merged_candidate_intervals.tsv")
    ap.add_argument("--top-n", type=int, default=6)
    ap.add_argument("--pad", type=int, default=250000, help="bp added to each side of an interval")
    ap.add_argument("--min-called-per-sex", type=int, default=3)
    ap.add_argument("--outdir", default="candidate_inspection")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    intervals = (pd.read_csv(args.intervals, sep="\t")
                 .sort_values(["MAX_COUNT", "N_WINDOWS", "MEAN_COUNT"], ascending=False)
                 .head(args.top_n))

    summary = []
    for rank, (_, r) in enumerate(intervals.iterrows(), start=1):
        start = max(0, int(r["START"]) - args.pad)
        end = int(r["END"]) + args.pad
        region = f"{r['CHROM']}:{start}-{end}"
        tag = os.path.join(args.outdir, f"{rank}_{r['SEX']}_{r['CHROM']}_{start}_{end}")
        print(f"==> {region}")

        run(["bcftools", "view", "-r", region, args.vcf, "-Oz", "-o", f"{tag}.vcf.gz"])
        run(["tabix", "-f", "-p", "vcf", f"{tag}.vcf.gz"])
        with open(f"{tag}.genotypes.tsv", "w") as out:
            run(["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n", f"{tag}.vcf.gz"], stdout=out)

        run([sys.executable, os.path.join(HERE, "inspect_candidate_region.py"),
             "--vcf-samples", args.vcf_samples, "--males", args.males, "--females", args.females,
             "--genotypes", f"{tag}.genotypes.tsv", "--prefix", tag,
             "--min-called-per-sex", str(args.min_called_per_sex)])

        with open(f"{tag}.summary.txt") as f:
            lines = dict(l.split(":", 1) for l in f.read().splitlines() if ":" in l)
        summary.append({
            "interval_class": r["SEX"], "CHROM": r["CHROM"], "START_PAD": start, "END_PAD": end,
            **{k.strip(): v.strip() for k, v in lines.items() if k.startswith(("Sites kept", "Strong"))},
            "site_summary_file": f"{tag}.site_summary.tsv",
        })

    path = os.path.join(args.outdir, "candidate_interval_summary.tsv")
    pd.DataFrame(summary).to_csv(path, sep="\t", index=False)
    print(f"Summary written to: {path}")


if __name__ == "__main__":
    main()
