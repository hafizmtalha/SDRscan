#!/usr/bin/env bash
# Stage 2: joint variant calling with freebayes at the declared ploidy.
#
# All samples are called together, so every sample gets a genotype at every
# variable site. (Calling samples one by one and merging the VCFs afterwards
# leaves samples that match the reference as *missing* instead of 0/0, which
# hides exactly the sex-specific sites this pipeline looks for.)
#
# To keep the load on shared servers bounded, the genome is split into chunks
# of CALL_CHUNK_BP and MAX_PARALLEL_JOBS single-threaded freebayes jobs run at
# a time (a single freebayes-parallel run with 60 jobs once made a 48-core
# server unresponsive).
#
# Ploidy comes from PLOIDY, or the sample sheet's "ploidy" column for mixed
# populations (passed to freebayes with --cnv-map). The samples and ploidies
# used are recorded in called_with.tsv; if they change, all chunks are re-called.
#
# Usage: scripts/02_call_variants.sh
source "$(dirname "$0")/common.sh"

BAMDIR="${OUTDIR}/01_align/bam"
OUT="${OUTDIR}/02_vcf"
CHUNKS="$OUT/chunks"
mkdir -p "$CHUNKS" "$OUT/logs"
CHUNK_BP="${CALL_CHUNK_BP:-1000000}"
JOBS="${MAX_PARALLEL_JOBS:-${MAX_PARALLEL_SAMPLES:-2}}"

# Samples, BAMs and ploidy for this run
BAMLIST="$OUT/bams.txt"; WANT="$OUT/called_with.new"
: > "$BAMLIST"; : > "$WANT"
for S in $(samples_of); do
  B="$BAMDIR/${S}.RG.dedup.bam"; [ -f "$B" ] || B="$BAMDIR/${S}.RG.bam"
  [ -f "$B" ] || die "no BAM for $S in $BAMDIR (run stage 1)"
  echo "$B" >> "$BAMLIST"
  printf '%s\t%s\n' "$S" "$(sample_ploidy "$S")" >> "$WANT"
done

if [ -f "$OUT/called_with.tsv" ] && ! cmp -s "$WANT" "$OUT/called_with.tsv"; then
  log "Samples or ploidy changed since the last call: re-calling all chunks"
  rm -f "$CHUNKS"/*.vcf.gz "$OUT/joint.vcf.gz" "$OUT/joint.vcf.gz.tbi"
fi
mv "$WANT" "$OUT/called_with.tsv"

if [ -s "$OUT/joint.vcf.gz" ]; then
  log "Joint VCF exists: $OUT/joint.vcf.gz"
  exit 0
fi

# Mixed ploidy -> per-sample copy-number map for freebayes
CNV_MAP=""
if [ "$(cut -f2 "$OUT/called_with.tsv" | sort -u | wc -l)" -gt 1 ]; then
  CNV_MAP="$OUT/cnv_map.tsv"
  cp "$OUT/called_with.tsv" "$CNV_MAP"
fi
log "Joint calling $(wc -l < "$BAMLIST") samples (ploidy: $(cut -f2 "$OUT/called_with.tsv" | sort -un | tr '\n' ' '))"

[ -f "${REF}.fai" ] || samtools faidx "$REF"
bedtools makewindows -g "${REF}.fai" -w "$CHUNK_BP" \
  | awk '{printf "%s:%d-%d\t%s_%d\n", $1, $2, $3, $1, $2}' > "$OUT/regions.txt"
log "$(wc -l < "$OUT/regions.txt") chunks of ${CHUNK_BP} bp, ${JOBS} at a time"

call_chunk() {
  local region="$1" name="$2" out cnv=()
  out="$CHUNKS/${name}.vcf.gz"
  [ -s "$out" ] && return 0
  [ -n "$CNV_MAP" ] && cnv=(--cnv-map "$CNV_MAP")
  freebayes -f "$REF" -p "$PLOIDY" "${cnv[@]}" -C "$MIN_ALT_COUNT" \
      --min-mapping-quality "$MIN_MAPQ" --min-base-quality "$MIN_BASEQ" \
      --genotype-qualities --use-mapping-quality \
      -r "$region" -L "$BAMLIST" 2> "$OUT/logs/${name}.log" \
    | bgzip > "$out.tmp" && mv "$out.tmp" "$out"
}
export -f call_chunk
export REF CHUNKS OUT BAMLIST MIN_ALT_COUNT MIN_MAPQ MIN_BASEQ PLOIDY CNV_MAP
xargs -P "$JOBS" -L 1 bash -c 'call_chunk "$@"' _ < "$OUT/regions.txt"

n_done=$(ls "$CHUNKS"/*.vcf.gz 2>/dev/null | wc -l)
[ "$n_done" -eq "$(wc -l < "$OUT/regions.txt")" ] || die "only $n_done of $(wc -l < "$OUT/regions.txt") chunks finished; see $OUT/logs"

# Concatenate in region order; drop records duplicated at chunk boundaries
cut -f2 "$OUT/regions.txt" | sed "s#^#$CHUNKS/#; s#\$#.vcf.gz#" > "$OUT/chunk_files.txt"
bcftools concat -f "$OUT/chunk_files.txt" -Ou \
  | bcftools sort -Ou -T "$OUT/sort_tmp" \
  | bcftools norm -d exact -Oz -o "$OUT/joint.vcf.gz"
tabix -f -p vcf "$OUT/joint.vcf.gz"
log "Done. $(bcftools index -n "$OUT/joint.vcf.gz") records in $OUT/joint.vcf.gz"
