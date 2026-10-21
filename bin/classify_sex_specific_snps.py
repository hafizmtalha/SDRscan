#!/usr/bin/env python3
"""
Classify sex-specific SNPs from a multi-sample genotype table, for any ploidy.

Each sample's genotype is reduced to ALT-present / ALT-absent / missing. A site is

  male_specific   if >= MIN_ALT_MALES males carry ALT and
                     <= MAX_ALT_OPPOSITE females carry ALT
  female_specific if the reverse holds

and both sexes have at least MIN_CALLED_* called genotypes.

Thresholds are counts (7) or fractions of that sex's sample size (0.8 = 80%,
rounded up), so the same settings carry over between studies of different
sizes. --max-alt-opposite is a count or a fraction of the opposite sex.

Ploidy must be declared with --ploidy N or --sample-ploidy FILE
(sample<TAB>ploidy, for mixed-ploidy populations). Genotypes with a different
number of alleles than declared are treated as missing and reported.

By default any genotype with at least one ALT allele counts as ALT-present.
--alt-dosage restricts this to the listed ALT copy numbers (1..ploidy), e.g.
for a tetraploid --alt-dosage 1 counts only 0/0/0/1, --alt-dosage 2 only 0/0/1/1.
--min-alt-dosage N counts genotypes with at least N ALT copies.

Inputs:
  --genotypes    bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n' merged.vcf.gz
  --vcf-samples  bcftools query -l merged.vcf.gz   (same order as genotype columns)
  --males / --females  one sample name per line

Output: TSV with CHROM POS CLASS male_called female_called male_alt female_alt
"""

import argparse
import sys

import pandas as pd

from gtutils import GenotypeParser, check_dosages, load_ploidy, read_list, threshold


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", default="merged.genotypes.tsv")
    ap.add_argument("--vcf-samples", default="vcf_samples.txt")
    ap.add_argument("--males", default="males.txt")
    ap.add_argument("--females", default="females.txt")
    ap.add_argument("--out", default="sex_specific_sites.tsv")
    pl = ap.add_mutually_exclusive_group(required=True)
    pl.add_argument("--ploidy", type=int, help="Ploidy of every sample")
    pl.add_argument("--sample-ploidy", help="TSV sample<TAB>ploidy")
    ap.add_argument("--min-called-males", default="0.75", help="count or fraction of males")
    ap.add_argument("--min-called-females", default="0.75", help="count or fraction of females")
    ap.add_argument("--min-alt-males", default="0.75", help="count or fraction of males")
    ap.add_argument("--min-alt-females", default="0.75", help="count or fraction of females")
    ap.add_argument("--max-alt-opposite", default="0",
                    help="max ALT-present samples in the opposite sex (count or fraction)")
    dose = ap.add_mutually_exclusive_group()
    dose.add_argument("--alt-dosage", type=int, nargs="+", default=None,
                      help="only these ALT copy numbers count as ALT-present")
    dose.add_argument("--min-alt-dosage", type=int, default=1,
                      help="at least this many ALT copies count as ALT-present (default 1)")
    args = ap.parse_args()

    vcf_samples = read_list(args.vcf_samples)
    males = read_list(args.males)
    females = read_list(args.females)
    male_set, female_set = set(males), set(females)

    if not vcf_samples:
        sys.exit("ERROR: vcf samples file is empty")
    if not males or not females:
        sys.exit("ERROR: males or females list is empty")
    if male_set & female_set:
        sys.exit(f"ERROR: samples in both sex lists: {sorted(male_set & female_set)}")
    missing = [s for s in male_set | female_set if s not in vcf_samples]
    if missing:
        sys.exit(f"ERROR: samples not found in VCF sample list: {missing}")

    ploidy = load_ploidy(vcf_samples, args.ploidy, args.sample_ploidy)
    check_dosages(args.alt_dosage, ploidy)
    dosages = set(args.alt_dosage) if args.alt_dosage else None
    parser = GenotypeParser(ploidy)

    nm, nf = len(males), len(females)
    min_called_m = threshold(args.min_called_males, nm)
    min_called_f = threshold(args.min_called_females, nf)
    min_alt_m = threshold(args.min_alt_males, nm)
    min_alt_f = threshold(args.min_alt_females, nf)
    max_opp_in_f = threshold(args.max_alt_opposite, nf)   # females allowed at male-specific sites
    max_opp_in_m = threshold(args.max_alt_opposite, nm)   # males allowed at female-specific sites
    ploidies = sorted(set(ploidy.values()))
    alt_rule = (f"dosage in {sorted(dosages)}" if dosages else f"dosage >= {args.min_alt_dosage}")
    print(f"Samples: {nm} males, {nf} females; ploidy {ploidies}; ALT-present = {alt_rule}")
    print(f"male_specific:   called M>={min_called_m}, F>={min_called_f}; ALT in M>={min_alt_m}, F<={max_opp_in_f}")
    print(f"female_specific: called M>={min_called_m}, F>={min_called_f}; ALT in F>={min_alt_f}, M<={max_opp_in_m}")

    sex_of = [("M" if s in male_set else "F" if s in female_set else None) for s in vcf_samples]
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
                sys.exit(f"ERROR at {chrom}:{pos}: {len(gts)} genotype columns, {len(vcf_samples)} samples")

            mc = fc = ma = fa = 0
            for s, sex, gt in zip(vcf_samples, sex_of, gts):
                if sex is None:
                    continue
                d = parser.dosage(gt, s)
                if d is None:
                    continue
                alt = (d in dosages) if dosages else (d >= args.min_alt_dosage)
                if sex == "M":
                    mc += 1
                    ma += alt
                else:
                    fc += 1
                    fa += alt

            if mc < min_called_m or fc < min_called_f:
                continue
            if ma >= min_alt_m and fa <= max_opp_in_f:
                cls = "male_specific"
                n_male += 1
            elif fa >= min_alt_f and ma <= max_opp_in_m:
                cls = "female_specific"
                n_female += 1
            else:
                continue
            rows.append([chrom, pos, cls, mc, fc, ma, fa])

    pd.DataFrame(rows, columns=["CHROM", "POS", "CLASS", "male_called", "female_called",
                                "male_alt", "female_alt"]).to_csv(args.out, sep="\t", index=False)
    parser.report(sys.stderr)
    print(f"Total sites scanned:   {total}")
    print(f"Male-specific sites:   {n_male}")
    print(f"Female-specific sites: {n_female}")
    print(f"Output written to:     {args.out}")


if __name__ == "__main__":
    main()
