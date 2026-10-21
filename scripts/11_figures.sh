#!/usr/bin/env bash
# Stage 11: summary figures.
#   - one Manhattan-style plot per SNP run (male- and female-specific windows)
#   - one Manhattan-style plot per k-mer mode with every reference as a track
#   - genome-wide and single-chromosome dotplots + chromosome correspondence
#
# To put SNP runs from two references (two OUTDIRs) on one figure, call
# bin/plot_manhattan_windows.py directly with one --track per file.
source "$(dirname "$0")/common.sh"

OUT="${OUTDIR}/11_figures"
mkdir -p "$OUT"
TAG="${WINDOW}_${STEP}"
COLOURS=(purple green orange blue red brown)
MARKERS=(o '^' s D v P)

# SNP runs
for entry in "${SNP_RUNS[@]}"; do
  run="${entry%%|*}"
  d="${OUTDIR}/05_snp_windows/$run"
  [ -f "$d/male_specific_windows_${TAG}.bed" ] || continue
  python "$BIN/plot_manhattan_windows.py" \
    --track "Male-specific SNPs=$d/male_specific_windows_${TAG}.bed:green:o" \
    --track "Female-specific SNPs=$d/female_specific_windows_${TAG}.bed:purple:^" \
    --chrom-regex "$CHROM_REGEX" --title "Genome-wide sex-specific SNP density ($run)" \
    --ylabel "Sex-specific SNPs per window" --outprefix "$OUT/snp_${run}_manhattan"
done

# k-mer runs: all references on one plot per mode
for entry in "${KMER_RUNS[@]}"; do
  mode="${entry%%|*}"
  tracks=()
  i=0
  for spec in "${KMER_REFS[@]}"; do
    tag="${spec%%=*}"
    f="${OUTDIR}/08_kmer/${mode}/windows/genome.${mode}_density.${tag}.bed"
    [ -f "$f" ] && tracks+=(--track "${tag}=${f}:${COLOURS[$i]}:${MARKERS[$i]}")
    i=$(( (i + 1) % ${#COLOURS[@]} ))
  done
  [ ${#tracks[@]} -gt 0 ] || continue
  python "$BIN/plot_manhattan_windows.py" "${tracks[@]}" \
    --chrom-regex "$CHROM_REGEX" --title "Genome-wide sex-specific k-mer density (${mode})" \
    --ylabel "${mode} per window" --outprefix "$OUT/kmer_${mode}_manhattan"
done

# Dotplots
DOT="${OUTDIR}/10_dotplots"
if [ -d "$DOT" ] && command -v Rscript > /dev/null; then
  for coords in "$DOT"/*/*.filtered*.coords; do
    [ -f "$coords" ] || continue
    name=$(basename "$(dirname "$coords")")
    Rscript "$PIPELINE_DIR/R/dotplot_genome.R" "$coords" "${name%%_vs_*}" "${name##*_vs_}" \
      "$OUT/${name}_genome" "$CHROM_REGEX"
  done
  for entry in "${DOTPLOT_CHROMS[@]}"; do
    IFS='|' read -r pair rchr qchr label <<< "$entry"
    coords=$(ls "$DOT/$pair"/*.filtered*.coords 2>/dev/null | head -n 1)
    [ -n "$coords" ] || { echo "No coords for $pair" >&2; continue; }
    Rscript "$PIPELINE_DIR/R/dotplot_chromosome.R" "$coords" "$rchr" "$qchr" \
      "${pair%%_vs_*}" "${pair##*_vs_}" "$OUT/${pair}_${label}" 10000 5
  done
fi

log "Done. Figures in $OUT"
