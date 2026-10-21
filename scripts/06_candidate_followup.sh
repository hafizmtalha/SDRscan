#!/usr/bin/env bash
# Stage 6: follow up candidate regions picked from the stage-05 output.
#
# For each region in CANDIDATES (CHROM:START-END), using the relaxed and
# strict SNP runs:
#   a. strict-site overlap with the relaxed candidate intervals (genome-wide)
#   b. zoomed density plot of the region's chromosome, region shaded
#   c. genotype inspection: many repeated sex-biased SNPs, or a few noisy ones?
#   d. PCA of the region: do samples split by sex? any sample misplaced?
#   e. sex-label permutation test on the region
#   f. chromosome-level enrichment permutation (region's chromosome vs genome)
#
# Usage: scripts/06_candidate_followup.sh [CHROM:START-END ...]   (default: CANDIDATES)
source "$(dirname "$0")/common.sh"

REGIONS=("$@")
[ ${#REGIONS[@]} -gt 0 ] || REGIONS=("${CANDIDATES[@]}")
[ ${#REGIONS[@]} -gt 0 ] || { echo "No candidate regions: set CANDIDATES in the config or pass CHROM:START-END" >&2; exit 1; }

IN="${OUTDIR}/03_merged"
RELAXED="${OUTDIR}/05_snp_windows/${RELAXED_RUN}"
STRICT="${OUTDIR}/05_snp_windows/${STRICT_RUN}"
OUT="${OUTDIR}/06_candidates"
mkdir -p "$OUT"
TAG="${WINDOW}_${STEP}"
[ "$FOCAL_SEX" = "M" ] && FOCAL=male || FOCAL=female
require_ploidy_table
require_both_sexes
GT_ARGS=(--vcf-samples "$IN/vcf_samples.txt" --males "$IN/males.txt" --females "$IN/females.txt"
         --sample-ploidy "$PLOIDY_TABLE")
STRONG_ALT="${STRONG_ALT:-0.75}"                  # "strong" site: ALT in >= 75% of the focal sex
STRONG_MAX_OPPOSITE="${STRONG_MAX_OPPOSITE:-1}"   # ... and in <= 1 of the other sex

# a. Do strict sites fall inside the relaxed candidate intervals?
if [ -d "$STRICT" ]; then
  python "$BIN/strict_overlap_with_candidates.py" \
    --intervals "$RELAXED/merged_candidate_intervals.tsv" \
    --strict-male "$STRICT/male_specific_sites.bed" --strict-female "$STRICT/female_specific_sites.bed" \
    --out "$OUT/strict_overlap_counts.tsv"
fi

for region in "${REGIONS[@]}"; do
  CHROM="${region%%:*}"; RANGE="${region#*:}"; START="${RANGE%-*}"; END="${RANGE#*-}"
  name="${CHROM}_${START}_${END}"
  d="$OUT/$name"
  mkdir -p "$d"
  log "Candidate $region"

  # b. zoomed density plot
  python "$BIN/plot_candidate_chromosome.py" --chrom "$CHROM" --highlight "$START" "$END" \
    --male "$RELAXED/male_specific_windows_${TAG}.bed" --female "$RELAXED/female_specific_windows_${TAG}.bed" \
    --outprefix "$d/${CHROM}_density"

  # genotypes of the region
  bcftools view -r "$region" "$IN/merged.filtered.snps.vcf.gz" -Oz -o "$d/region.vcf.gz"
  tabix -f -p vcf "$d/region.vcf.gz"
  bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$d/region.vcf.gz" > "$d/region.genotypes.tsv"

  # c. genotype inspection
  python "$BIN/inspect_candidate_region.py" --genotypes "$d/region.genotypes.tsv" "${GT_ARGS[@]}" \
    --strong-alt "$STRONG_ALT" --strong-max-opposite "$STRONG_MAX_OPPOSITE" --prefix "$d/inspect"

  # d. PCA, presence and (for polyploids) dosage encoding; skipped if too few sites
  for enc in presence $([ "$(max_ploidy)" -gt 2 ] && echo dosage); do
    python "$BIN/pca_from_gt_presence.py" --genotypes "$d/region.genotypes.tsv" "${GT_ARGS[@]}" \
      --encoding "$enc" --outprefix "$d/pca_${enc}" --min-sites 100 --max-sites 20000 \
      || log "PCA ($enc) skipped for $region"
  done

  # e. sex-label permutation (strong focal-sex-biased sites in the region)
  python "$BIN/permute_sex_labels.py" --genotypes "$d/region.genotypes.tsv" "${GT_ARGS[@]}" \
    --focal-sex "$FOCAL_SEX" --chrom "$CHROM" --start "$START" --end "$END" --nperm "$PERMUTATIONS" \
    --min-called "$STRONG_ALT" --min-alt "$STRONG_ALT" --max-opposite "$STRONG_MAX_OPPOSITE" \
    --outprefix "$d/permutation"

  # f. chromosome enrichment of focal-sex-specific windows
  python "$BIN/window_enrichment_permutation.py" \
    --windows "$RELAXED/${FOCAL}_specific_windows_${TAG}.bed" --chrom "$CHROM" \
    --nperm "$PERMUTATIONS" --outprefix "$d/chrom_enrichment"
done

log "Done. Outputs in $OUT"
