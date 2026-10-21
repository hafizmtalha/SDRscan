#!/usr/bin/env bash
# Stage 2: per-sample variant calling with freebayes at the declared ploidy.
#
# Ploidy comes from PLOIDY in the config, or the sample sheet's "ploidy"
# column for mixed-ploidy populations. The ploidy used for each sample is
# written next to its VCF (<sample>.ploidy); if PLOIDY changes later the sample
# is re-called instead of silently reusing calls made at another ploidy.
#
# Samples are called one by one (MAX_PARALLEL_SAMPLES at a time) rather than
# with freebayes-parallel over all BAMs, which overloaded a shared server.
#
# Usage: scripts/02_call_variants.sh [sample ...]   (default: all samples)
source "$(dirname "$0")/common.sh"

BAMDIR="${OUTDIR}/01_align/bam"
OUT="${OUTDIR}/02_vcf"
mkdir -p "$OUT/logs"

SAMPLE_LIST=("$@")
[ ${#SAMPLE_LIST[@]} -gt 0 ] || mapfile -t SAMPLE_LIST < <(samples_of)

JOBS="$OUT/jobs.tsv"
: > "$JOBS"
for S in "${SAMPLE_LIST[@]}"; do
  P=$(sample_ploidy "$S")
  if [ -s "$OUT/${S}.vcf.gz" ] && [ "$(cat "$OUT/${S}.ploidy" 2>/dev/null)" = "$P" ]; then
    log "$S: already called at ploidy $P, skipping"
    continue
  fi
  printf '%s\t%s\n' "$S" "$P" >> "$JOBS"
done
log "Calling $(wc -l < "$JOBS") samples (ploidy: $(cut -f2 "$JOBS" | sort -u | tr '\n' ' '))"

call_one() {
  local S="$1" P="$2" BAM
  BAM="$BAMDIR/${S}.RG.dedup.bam"
  [ -f "$BAM" ] || BAM="$BAMDIR/${S}.RG.bam"
  echo "$S: calling at ploidy $P"
  rm -f "$OUT/${S}.vcf.gz" "$OUT/${S}.ploidy"
  freebayes -f "$REF" -p "$P" -C "$MIN_ALT_COUNT" \
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
  echo "$P" > "$OUT/${S}.ploidy"
}
export -f call_one
export OUT BAMDIR REF MIN_ALT_COUNT MIN_MAPQ MIN_BASEQ

[ -s "$JOBS" ] && xargs -P "$MAX_PARALLEL_SAMPLES" -L 1 bash -c 'call_one "$@"' _ < "$JOBS"
log "Done. Per-sample VCFs in $OUT"
