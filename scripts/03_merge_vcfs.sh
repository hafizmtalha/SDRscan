#!/usr/bin/env bash
# Stage 3: normalise and filter the joint VCF, then export the genotype table
# used by every SNP-based analysis.
#
#   bcftools norm (split multiallelics) -> QUAL >= MIN_QUAL && DP >= MIN_DP
#   -> SNPs only -> biallelic only -> genotype table + sample order + sex lists
source "$(dirname "$0")/common.sh"

IN="${OUTDIR}/02_vcf"
OUT="${OUTDIR}/03_merged"
mkdir -p "$OUT"
[ -s "$IN/joint.vcf.gz" ] || die "$IN/joint.vcf.gz missing; run stage 2"

# The genotypes must have been called at the ploidy currently configured
write_sample_ploidy "$OUT/sample_ploidy.tsv"
cmp -s "$OUT/sample_ploidy.tsv" "$IN/called_with.tsv" || \
  die "samples or ploidy differ from those used for calling ($IN/called_with.tsv); re-run stage 2"
log "Ploidy per sample: $OUT/sample_ploidy.tsv"

log "Normalising and filtering (QUAL>=${MIN_QUAL}, DP>=${MIN_DP}, biallelic SNPs)"
bcftools norm -m -both -f "$REF" "$IN/joint.vcf.gz" -Ou \
  | bcftools filter -i "QUAL>=${MIN_QUAL} && INFO/DP>=${MIN_DP}" -Ou \
  | bcftools view -v snps -m2 -M2 -Oz -o "$OUT/merged.filtered.snps.vcf.gz"
tabix -f -p vcf "$OUT/merged.filtered.snps.vcf.gz"

bcftools query -l "$OUT/merged.filtered.snps.vcf.gz" > "$OUT/vcf_samples.txt"
bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$OUT/merged.filtered.snps.vcf.gz" > "$OUT/merged.genotypes.tsv"
samples_of M > "$OUT/males.txt"
samples_of F > "$OUT/females.txt"

missing=$(comm -23 <(samples_of | sort) <(sort "$OUT/vcf_samples.txt"))
[ -z "$missing" ] || die "samples missing from the VCF: $missing (read-group sample names must equal sample IDs)"

log "SNPs kept: $(wc -l < "$OUT/merged.genotypes.tsv")"
log "Done. Outputs in $OUT"
