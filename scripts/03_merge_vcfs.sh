#!/usr/bin/env bash
# Stage 3: normalise, rename, merge and filter per-sample VCFs, then export
# the genotype table used by every SNP-based analysis.
#
#   bcftools norm -> reheader (sample name = sample ID) -> merge
#   -> QUAL >= MIN_QUAL && DP >= MIN_DP -> SNPs only -> biallelic only
#   -> genotype table + sample order + male/female lists
source "$(dirname "$0")/common.sh"

IN="${OUTDIR}/02_vcf"
OUT="${OUTDIR}/03_merged"
mkdir -p "$OUT/norm"

mapfile -t SAMPLE_LIST < <(samples_of)

# Every VCF must have been called at the ploidy currently declared for its sample
for S in "${SAMPLE_LIST[@]}"; do
  called=$(cat "$IN/${S}.ploidy" 2>/dev/null || echo "unknown")
  want=$(sample_ploidy "$S")
  [ "$called" = "$want" ] || die "$S was called at ploidy '$called' but the config says $want; re-run stage 2"
done
write_sample_ploidy "$OUT/sample_ploidy.tsv"
log "Ploidy per sample written to $OUT/sample_ploidy.tsv"

NORMED=()
for S in "${SAMPLE_LIST[@]}"; do
  out="$OUT/norm/${S}.norm.vcf.gz"
  if [ ! -s "$out" ] || [ "$IN/${S}.vcf.gz" -nt "$out" ]; then
    log "$S: normalise + rename"
    # freebayes names the sample after the read group; force it to the sample ID
    # so bcftools merge does not see duplicate or "unknown" sample names.
    bcftools norm -m -both -f "$REF" "$IN/${S}.vcf.gz" -Oz -o "$OUT/norm/${S}.tmp.vcf.gz"
    echo "$S" > "$OUT/norm/${S}.name"
    bcftools reheader -s "$OUT/norm/${S}.name" -o "$out" "$OUT/norm/${S}.tmp.vcf.gz"
    rm -f "$OUT/norm/${S}.tmp.vcf.gz" "$OUT/norm/${S}.name"
    tabix -f -p vcf "$out"
  fi
  NORMED+=("$out")
done

log "Merging ${#NORMED[@]} samples"
bcftools merge --threads "$THREADS" "${NORMED[@]}" -Oz -o "$OUT/merged.vcf.gz"
tabix -f -p vcf "$OUT/merged.vcf.gz"

log "Filtering (QUAL>=${MIN_QUAL}, DP>=${MIN_DP}, biallelic SNPs)"
bcftools filter -i "QUAL>=${MIN_QUAL} && INFO/DP>=${MIN_DP}" "$OUT/merged.vcf.gz" -Ou \
  | bcftools view -v snps -m2 -M2 -Oz -o "$OUT/merged.filtered.snps.vcf.gz"
tabix -f -p vcf "$OUT/merged.filtered.snps.vcf.gz"

bcftools query -l "$OUT/merged.filtered.snps.vcf.gz" > "$OUT/vcf_samples.txt"
bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$OUT/merged.filtered.snps.vcf.gz" > "$OUT/merged.genotypes.tsv"
samples_of M > "$OUT/males.txt"
samples_of F > "$OUT/females.txt"

log "Sites kept: $(wc -l < "$OUT/merged.genotypes.tsv")"
log "Done. Outputs in $OUT"
