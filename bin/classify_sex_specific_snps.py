#!/usr/bin/env python3
"""
Classify sex-specific SNPs from a multi-sample genotype table.

Each sample's genotype is reduced to ALT-present / ALT-absent / missing, which
works for any ploidy (0/1, 0/0/0/1, 0|1|1|1, ...). A site is

  male_specific   if >= MIN_ALT_MALES males carry ALT and
                     <= MAX_ALT_IN_OPPOSITE_SEX females carry ALT
  female_specific if the reverse holds

and both sexes have at least MIN_CALLED_* called genotypes.

By default any non-zero allele counts as ALT-present. --alt-dosage restricts
this to genotypes carrying exactly the listed number of ALT copies, e.g.
--alt-dosage 1 counts only 0/0/0/1, --alt-dosage 2 only 0/0/1/1.

Inputs:
  --genotypes    bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n' merged.vcf.gz
  --vcf-samples  bcftools query -l merged.vcf.gz   (same order as genotype columns)
  --males / --females  one sample name per line

Output: TSV with CHROM POS CLASS male_called female_called male_alt female_alt
"""

import argparse
import re
import sys

import pandas as pd

MISSING = {".", "./.", ".|.", "./././.", ".|.|.|.", "././././.", ".|.|.|.|."}


def read_list(path):
    with open(path) as f:
        return [x.strip() for x in f if x.strip()]


def gt_called(gt):
    if gt is None:
        return False
    return gt.strip() not in MISSING


def alt_dosage(gt):
    """Number of non-reference alleles in a called genotype."""
    alleles = re.split(r"[\/|]", gt.strip())
    return sum(1 for a in alleles if a not in {"", ".", "0"})


def gt_has_alt(gt, dosages=None):
    if not gt_called(gt):
        return False
    d = alt_dosage(gt)
    if dosages is None:
        return d > 0
    return d in dosages


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", default="merged.genotypes.tsv")
    ap.add_argument("--vcf-samples", default="vcf_samples.txt")
    ap.add_argument("--males", default="males.txt")
    ap.add_argument("--females", default="females.txt")
    ap.add_argument("--out", default="sex_specific_sites.tsv")
    ap.add_argument("--min-called-males", type=int, default=4)
    ap.add_argument("--min-called-females", type=int, default=4)
    ap.add_argument("--min-alt-males", type=int, default=4)
    ap.add_argument("--min-alt-females", type=int, default=4)
    ap.add_argument("--max-alt-opposite", type=int, default=0,
                    help="Max ALT-present samples allowed in the opposite sex")
    ap.add_argument("--alt-dosage", type=int, nargs="+", default=None,
                    help="Only count genotypes with exactly these ALT copy numbers as ALT-present")
    args = ap.parse_args()

    vcf_samples = read_list(args.vcf_samples)
    male_samples = set(read_list(args.males))
    female_samples = set(read_list(args.females))
    dosages = set(args.alt_dosage) if args.alt_dosage else None

    if not vcf_samples:
        sys.exit("ERROR: vcf samples file is empty")
    if not male_samples or not female_samples:
        sys.exit("ERROR: males or females list is empty")
    overlap = male_samples & female_samples
    if overlap:
        sys.exit(f"ERROR: samples in both sex lists: {sorted(overlap)}")
    missing = [s for s in male_samples | female_samples if s not in vcf_samples]
    if missing:
        sys.exit(f"ERROR: samples not found in VCF sample list: {missing}")

    sample_idx = dict(enumerate(vcf_samples))
    rows = []
    total = n_male = n_female = 0

    with open(args.genotypes) as f:
        for line in f:
            total += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            gts = parts[2:]
            if len(gts) != len(vcf_samples):
                sys.exit(f"ERROR at {chrom}:{pos}: {len(gts)} genotype columns, "
                         f"{len(vcf_samples)} samples")

            male_called = female_called = male_alt = female_alt = 0
            for i, gt in enumerate(gts):
                s = sample_idx[i]
                called = gt_called(gt)
                alt = gt_has_alt(gt, dosages)
                if s in male_samples:
                    male_called += called
                    male_alt += alt
                elif s in female_samples:
                    female_called += called
                    female_alt += alt

            if male_called < args.min_called_males or female_called < args.min_called_females:
                continue

            if male_alt >= args.min_alt_males and female_alt <= args.max_alt_opposite:
                cls = "male_specific"
                n_male += 1
            elif female_alt >= args.min_alt_females and male_alt <= args.max_alt_opposite:
                cls = "female_specific"
                n_female += 1
            else:
                continue
            rows.append([chrom, pos, cls, male_called, female_called, male_alt, female_alt])

    df = pd.DataFrame(rows, columns=["CHROM", "POS", "CLASS", "male_called", "female_called",
                                     "male_alt", "female_alt"])
    df.to_csv(args.out, sep="\t", index=False)

    print(f"Total sites scanned:   {total}")
    print(f"Male-specific sites:   {n_male}")
    print(f"Female-specific sites: {n_female}")
    print(f"Output written to:     {args.out}")


if __name__ == "__main__":
    main()
