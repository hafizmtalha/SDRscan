#!/usr/bin/env bash
# Stage 12: read-depth and heterozygosity scan between the sexes.
#
# Works for any ploidy and any number of samples. It is the main analysis
# for a 1 male + 1 female pilot (ALLOW_SINGLE_PER_SEX=1), and a ploidy-aware
# complement to findZX for larger designs. Its per-sample depth table is also
# an input for sex assignment (stage 13).
#
#   samtools bedcov over DEPTH_WINDOW windows (MAPQ >= MIN_MAPQ, main chromosomes)
#   -> depth normalised per sample -> log2(male/female) per window
#   -> heterozygosity per sample/window from the merged genotypes (if stage 3 ran)
#   -> candidate windows/intervals, per-chromosome summary, genome plot
source "$(dirname "$0")/common.sh"
require_both_sexes pilot_ok

OUT="${OUTDIR}/12_depth_het"
BAMDIR="${OUTDIR}/01_align/bam"
mkdir -p "$OUT"

# Every sample with a BAM, including unknown-sex ones (needed for stage 13)
: > "$OUT/bam_samples.txt"
BAMS=()
for S in $(samples_of); do
  B="$BAMDIR/${S}.RG.dedup.bam"; [ -f "$B" ] || B="$BAMDIR/${S}.RG.bam"
  [ -f "$B" ] || { log "no BAM for $S, skipping"; continue; }
  BAMS+=("$B"); echo "$S" >> "$OUT/bam_samples.txt"
done
[ ${#BAMS[@]} -gt 0 ] || die "no BAMs found in $BAMDIR (run stage 1)"

if [ ! -s "$OUT/bedcov.tsv" ] || [ "$OUT/bam_samples.txt" -nt "$OUT/bedcov.tsv" ]; then
  bedtools makewindows -g "${REF}.fai" -w "$DEPTH_WINDOW" \
    | awk -v re="$CHROM_REGEX" '$1 ~ re' > "$OUT/windows.bed"
  log "Depth over $(wc -l < "$OUT/windows.bed") windows for ${#BAMS[@]} samples"
  samtools bedcov -Q "$MIN_MAPQ" "$OUT/windows.bed" "${BAMS[@]}" > "$OUT/bedcov.tsv"
fi

# Expected depth ratio for a region missing one of P copies is (P-1)/P; "auto"
# uses half of that log2 drop as the threshold.
THR="${DEPTH_LOG2_THRESHOLD:-auto}"
if [ "$THR" = "auto" ]; then
  THR=$(awk -v p="$(max_ploidy)" 'BEGIN{printf "%.3f", -0.5*log((p-1)/p)/log(2)}')
  log "Depth threshold (auto, ploidy $(max_ploidy)): |log2(M/F)| >= $THR"
fi

HET_ARGS=()
if [ -s "${OUTDIR}/03_merged/merged.genotypes.tsv" ] && [ -s "$PLOIDY_TABLE" ]; then
  HET_ARGS=(--genotypes "${OUTDIR}/03_merged/merged.genotypes.tsv"
            --vcf-samples "${OUTDIR}/03_merged/vcf_samples.txt" --sample-ploidy "$PLOIDY_TABLE")
else
  log "No merged genotypes (stage 3): depth only"
fi

python "$BIN/depth_het_scan.py" --bedcov "$OUT/bedcov.tsv" --bam-samples "$OUT/bam_samples.txt" \
  --sample-sheet "$SAMPLES" "${HET_ARGS[@]}" \
  --log2-threshold "$THR" --het-threshold "$HET_DIFF_THRESHOLD" \
  --outprefix "$OUT/depth_het" | tee "$OUT/depth_het.log"

log "Done. Outputs in $OUT"
