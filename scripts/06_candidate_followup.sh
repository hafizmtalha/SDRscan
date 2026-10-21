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
  python "$BIN/inspect_candidate_region.py" --genotypes "$d/region.genotypes.tsv" \
    --vcf-samples "$IN/vcf_samples.txt" --males "$IN/males.txt" --females "$IN/females.txt" \
    --prefix "$d/inspect"

  # d. PCA (needs enough sites; skipped with a message otherwise)
  python "$BIN/pca_from_gt_presence.py" --genotypes "$d/region.genotypes.tsv" \
    --vcf-samples "$IN/vcf_samples.txt" --males "$IN/males.txt" --females "$IN/females.txt" \
    --outprefix "$d/pca" --min-sites 100 --max-sites 20000 || log "PCA skipped for $region"

  # e. sex-label permutation (strong focal-sex-biased sites in the region)
  if [ "$FOCAL_SEX" = "M" ]; then
    PM="$IN/males.txt"; PF="$IN/females.txt"
  else
    PM="$IN/females.txt"; PF="$IN/males.txt"     # swap roles to test female bias
  fi
  n_focal=$(wc -l < "$PM"); n_other=$(wc -l < "$PF")
  python "$BIN/permute_sex_labels.py" --genotypes "$d/region.genotypes.tsv" \
    --vcf-samples "$IN/vcf_samples.txt" --males "$PM" --females "$PF" \
    --chrom "$CHROM" --start "$START" --end "$END" --nperm "$PERMUTATIONS" \
    --min-called-male $((n_focal * 3 / 4)) --min-called-female $((n_other * 3 / 4)) \
    --min-male-alt $((n_focal * 3 / 4)) --max-female-alt 1 \
    --outprefix "$d/permutation"

  # f. chromosome enrichment of focal-sex-specific windows
  python "$BIN/window_enrichment_permutation.py" \
    --windows "$RELAXED/${FOCAL}_specific_windows_${TAG}.bed" --chrom "$CHROM" \
    --nperm "$PERMUTATIONS" --outprefix "$d/chrom_enrichment"
done

log "Done. Outputs in $OUT"
