#!/usr/bin/env bash
# Stage 2: per-sample variant calling with freebayes.
#
# Samples are called one by one (MAX_PARALLEL_SAMPLES at a time) rather than
# with freebayes-parallel over all BAMs, which overloaded a shared server.
# Set PLOIDY in the config (2 = diploid, 4 = tetraploid).
#
# Usage: scripts/02_call_variants.sh [sample ...]   (default: all samples)
source "$(dirname "$0")/common.sh"

BAMDIR="${OUTDIR}/01_align/bam"
OUT="${OUTDIR}/02_vcf"
mkdir -p "$OUT/logs"

SAMPLE_LIST=("$@")
[ ${#SAMPLE_LIST[@]} -gt 0 ] || mapfile -t SAMPLE_LIST < <(samples_of)

call_one() {
  local S="$1" BAM
  BAM="$BAMDIR/${S}.RG.dedup.bam"
  [ -f "$BAM" ] || BAM="$BAMDIR/${S}.RG.bam"
  if [ -s "$OUT/${S}.vcf.gz" ]; then
    echo "$S: already called, skipping"
    return
  fi
  echo "$S: calling (ploidy ${PLOIDY})"
  freebayes -f "$REF" -p "$PLOIDY" -C "$MIN_ALT_COUNT" \
    --min-mapping-quality "$MIN_MAPQ" --min-base-quality "$MIN_BASEQ" \
    --genotype-qualities --use-mapping-quality \
    "$BAM" > "$OUT/${S}.vcf" 2> "$OUT/logs/${S}.freebayes.log"
  # A 0-byte VCF means freebayes failed (we hit this with a truncated BAM).
  if [ ! -s "$OUT/${S}.vcf" ]; then
    echo "ERROR: empty VCF for $S, see $OUT/logs/${S}.freebayes.log" >&2
    return 1
  fi
  bgzip -f "$OUT/${S}.vcf"
  tabix -f -p vcf "$OUT/${S}.vcf.gz"
}
export -f call_one
export OUT BAMDIR REF PLOIDY MIN_ALT_COUNT MIN_MAPQ MIN_BASEQ

printf '%s\n' "${SAMPLE_LIST[@]}" | xargs -P "$MAX_PARALLEL_SAMPLES" -I{} bash -c 'call_one "$@"' _ {}
log "Done. Per-sample VCFs in $OUT"
