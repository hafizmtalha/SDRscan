#!/usr/bin/env bash
# Run the pipeline stages in order. Every stage skips work whose output
# already exists, so it is safe to re-run after fixing a failure.
#
# Usage:
#   ./run_pipeline.sh                    # stages 1-5 and 12 (everything that needs no decisions)
#   ./run_pipeline.sh all                # 1-5, 8-12; 6-7 if CANDIDATES is set; 13 if any sex is "U"
#   ./run_pipeline.sh 13                 # assign sex to unknown samples / check recorded sexes
#   ./run_pipeline.sh 3 5 6              # only these stages
#   CONFIG=config/male_ref.env ./run_pipeline.sh 1 2 3 5
#
# Stages
#   1  QC, trimming, alignment            scripts/01_qc_align.sh
#   2  per-sample variant calling          scripts/02_call_variants.sh
#   3  merge + filter VCFs                 scripts/03_merge_vcfs.sh
#   4  findZX runs (XY / ZW / shuffled)    scripts/04_findzx.sh
#   5  sex-specific SNP windows (SNP_RUNS) scripts/05_sex_specific_snp_windows.sh
#   6  candidate follow-up (CANDIDATES)    scripts/06_candidate_followup.sh
#   7  leave-one-out (CANDIDATES)          scripts/07_leave_one_out.sh
#   8  sex-specific k-mers (KMER_RUNS)     scripts/08_kmer_sex_specific.sh
#   9  k-mer run overlaps                  scripts/09_compare_kmer_runs.sh
#   10 assembly dotplots (DOTPLOT_PAIRS)   scripts/10_nucmer_dotplots.sh
#   11 summary figures                     scripts/11_figures.sh
#   12 depth & heterozygosity scan         scripts/12_depth_het_scan.sh   (also the 1 + 1 pilot analysis)
#   13 sex assignment / sex check          scripts/13_assign_sex.sh
#
# With fewer than 2 samples of a sex (pilot, ALLOW_SINGLE_PER_SEX=1) only
# stages 1-4, 10 and 12 run; the others are skipped with an explanation.
set -euo pipefail
cd "$(dirname "$0")"
S=scripts
source "$S/common.sh"

sample_size_report
memory_report
STAGES=("$@")
if [ ${#STAGES[@]} -eq 0 ]; then
  STAGES=(1 2 3 4 5 12)
elif [ "${STAGES[0]}" = "all" ]; then
  STAGES=(1 2 3 4 5)
  [ ${#CANDIDATES[@]} -gt 0 ] && STAGES+=(6 7)
  STAGES+=(8 9 10 11 12)
  [ "$(n_of U)" -gt 0 ] && STAGES+=(13)
fi

# Which stages the design supports
#   >= 2 per sex : all
#   pilot (1 of a sex) : no stages 5-9, 11 (4 and 12 need ALLOW_SINGLE_PER_SEX=1)
#   a sex missing : only stages that do not compare sexes (1-3, 10, 13)
SKIP=""
if [ "$(n_of M)" -lt 1 ] || [ "$(n_of F)" -lt 1 ]; then
  SKIP=" 4 5 6 7 8 9 11 12 "; WHY="needs samples of both sexes"
elif is_pilot; then
  SKIP=" 5 6 7 8 9 11 "; WHY="needs >= $MIN_PER_SEX samples per sex (pilot design)"
fi

for stage in "${STAGES[@]}"; do
  if [[ "$SKIP" == *" $stage "* ]]; then
    log "===== Stage $stage skipped: $WHY ====="
    continue
  fi
  log "===== Stage $stage ====="
  case "$stage" in
    1) "$S/01_qc_align.sh" ;;
    2) "$S/02_call_variants.sh" ;;
    3) "$S/03_merge_vcfs.sh" ;;
    4) "$S/04_findzx.sh" ;;
    5)
      for entry in "${SNP_RUNS[@]}"; do
        # shellcheck disable=SC2086
        "$S/05_sex_specific_snp_windows.sh" "${entry%%|*}" ${entry#*|}
      done
      # Dosage-class runs, derived from the ploidy the genotypes were called with
      classes="${SNP_DOSAGE_CLASSES:-}"
      P=$(max_ploidy)
      [ "$classes" = "auto" ] && classes=$(seq -s ' ' 1 $((P - 1)))
      if [ -n "$classes" ]; then
        base_opts=""
        for entry in "${SNP_RUNS[@]}"; do
          [ "${entry%%|*}" = "$SNP_DOSAGE_BASE" ] && base_opts="${entry#*|}"
        done
        [ -n "$base_opts" ] || die "SNP_DOSAGE_BASE '$SNP_DOSAGE_BASE' is not a SNP_RUNS name"
        for c in $classes; do
          [ "$c" -ge 1 ] && [ "$c" -le "$P" ] || die "dosage class $c outside 1..$P (ploidy)"
          # shellcheck disable=SC2086
          "$S/05_sex_specific_snp_windows.sh" "${SNP_DOSAGE_BASE}_dose${c}" $base_opts --alt-dosage "$c"
        done
      fi
      log "Now inspect ${OUTDIR}/05_snp_windows/*/sex_specific_snp_windows.png and"
      log "merged_candidate_intervals.tsv, then set CANDIDATES in the config and run stages 6 7."
      ;;
    6) "$S/06_candidate_followup.sh" ;;
    7)
      if [ "$(n_of "$FOCAL_SEX")" -lt 3 ]; then
        log "Skipping leave-one-out: needs >= 3 samples of the focal sex ($FOCAL_SEX)"
        continue
      fi
      relaxed_opts=""
      for entry in "${SNP_RUNS[@]}"; do
        [ "${entry%%|*}" = "$RELAXED_RUN" ] && relaxed_opts="${entry#*|}"
      done
      for region in "${CANDIDATES[@]}"; do
        c="${region%%:*}"; r="${region#*:}"
        # shellcheck disable=SC2086
        "$S/07_leave_one_out.sh" "$FOCAL_SEX" "$c" "${r%-*}" "${r#*-}" $relaxed_opts
        dest="${OUTDIR}/07_leave_one_out_${c}_${r%-*}_${r#*-}"
        rm -rf "$dest"
        mv "${OUTDIR}/07_leave_one_out" "$dest"
      done
      ;;
    8)
      for entry in "${KMER_RUNS[@]}"; do
        "$S/08_kmer_sex_specific.sh" "${entry%%|*}" "${entry#*|}" "${KMER_REFS[@]}"
      done
      ;;
    9)
      beds=()
      for entry in "${KMER_RUNS[@]}"; do
        mode="${entry%%|*}"
        for spec in "${KMER_REFS[@]}"; do
          tag="${spec%%=*}"
          f="${OUTDIR}/08_kmer/${mode}/windows/kmer_${mode}_candidates.${tag}.bed"
          [ -f "$f" ] && beds+=("${mode}_${tag}=${f}")
        done
      done
      if [ -n "${KMER_NAME_NORMALISE:-}" ]; then
        "$S/09_compare_kmer_runs.sh" "${OUTDIR}/09_kmer_overlap" --rename "$KMER_NAME_NORMALISE" "${beds[@]}"
      else
        "$S/09_compare_kmer_runs.sh" "${OUTDIR}/09_kmer_overlap" "${beds[@]}"
      fi
      ;;
    10)
      THREADS="$THREADS" MIN_IDENTITY="$DOTPLOT_MIN_IDENTITY" MIN_LENGTH="$DOTPLOT_MIN_LENGTH" \
        CHROM_REGEX="$CHROM_REGEX" "$S/10_nucmer_dotplots.sh" "$DOTPLOT_PAIRS" "${OUTDIR}/10_dotplots"
      ;;
    11) "$S/11_figures.sh" ;;
    12) "$S/12_depth_het_scan.sh" ;;
    13) "$S/13_assign_sex.sh" ;;
    *) echo "Unknown stage: $stage" >&2; exit 1 ;;
  esac
done
log "Finished stages: ${STAGES[*]}"
