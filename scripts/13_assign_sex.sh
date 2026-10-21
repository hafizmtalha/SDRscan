#!/usr/bin/env bash
# Stage 13: assign sex to samples (sex "U" in the sample sheet) and check
# recorded sexes, using known sex-linked loci.
#
# Evidence used when available:
#   depth  depth over SEXING_REGIONS (default: CANDIDATES), from stage 12
#   snp    genotypes at SEXING_MARKERS (default: the strict stage-5 run's sites)
#   kmer   share of the stage-8 MSK / FSK k-mer sets found in each sample
#
# The loci can come from this project or a previous one: point SEXING_REGIONS,
# SEXING_MARKERS, SEXING_MSK_DB and SEXING_FSK_DB at earlier results to sex new
# samples of the same species without rediscovering the region.
source "$(dirname "$0")/common.sh"

OUT="${OUTDIR}/13_sexing"
mkdir -p "$OUT"
ARGS=(--sample-sheet "$SAMPLES" --outprefix "$OUT/sex_assignment")
[ -n "${SEX_SYSTEM:-}" ] && ARGS+=(--system "$SEX_SYSTEM")

# depth
REGIONS="${SEXING_REGIONS:-}"
if [ -z "$REGIONS" ] && [ ${#CANDIDATES[@]} -gt 0 ]; then
  REGIONS="$OUT/regions_from_candidates.bed"
  for r in "${CANDIDATES[@]}"; do
    c="${r%%:*}"; x="${r#*:}"; printf '%s\t%s\t%s\n' "$c" "${x%-*}" "${x#*-}"
  done > "$REGIONS"
fi
DEPTH="${OUTDIR}/12_depth_het/depth_het.per_sample_window_depth.tsv"
if [ -n "$REGIONS" ] && [ -s "$DEPTH" ]; then
  ARGS+=(--regions "$REGIONS" --depth "$DEPTH")
else
  log "depth evidence skipped (needs SEXING_REGIONS or CANDIDATES, and stage 12)"
fi

# snp markers
MARKERS="${SEXING_MARKERS:-${OUTDIR}/05_snp_windows/${STRICT_RUN}/sex_specific_sites.tsv}"
if [ -s "$MARKERS" ] && [ -s "${OUTDIR}/03_merged/merged.genotypes.tsv" ]; then
  ARGS+=(--markers "$MARKERS" --genotypes "${OUTDIR}/03_merged/merged.genotypes.tsv"
         --vcf-samples "${OUTDIR}/03_merged/vcf_samples.txt" --sample-ploidy "$PLOIDY_TABLE")
else
  log "SNP marker evidence skipped (needs SEXING_MARKERS or stage 5, and stage 3)"
fi

# k-mer sets
SFX="k${K}.ci${KMER_MIN_COUNT}"
KDIR="${OUTDIR}/08_kmer/kmc"
MSK="${SEXING_MSK_DB:-${OUTDIR}/08_kmer/sets/MSK/MSK}"
FSK="${SEXING_FSK_DB:-${OUTDIR}/08_kmer/sets/FSK/FSK}"
kmer_total() { kmc_tools info "$1" 2>/dev/null | awk -F: '/total k-mers/{gsub(/ /,"",$2); print $2}'; }
if command -v kmc_tools > /dev/null && { [ -f "${MSK}.kmc_pre" ] || [ -f "${FSK}.kmc_pre" ]; }; then
  echo -e "sample\tmsk_fraction\tfsk_fraction" > "$OUT/kmer_fractions.tsv"
  nm=$( [ -f "${MSK}.kmc_pre" ] && kmer_total "$MSK" || echo 0 )
  nf=$( [ -f "${FSK}.kmc_pre" ] && kmer_total "$FSK" || echo 0 )
  for S in $(samples_of); do
    [ -f "$KDIR/${S}.${SFX}.kmc_pre" ] || { log "no k-mer database for $S (run stage 8)"; continue; }
    fm=NA; ff=NA
    if [ "${nm:-0}" -gt 0 ]; then
      kmc_tools -t"$THREADS" simple "$MSK" "$KDIR/${S}.${SFX}" intersect "$OUT/_tmp_m" > /dev/null
      fm=$(awk -v a="$(kmer_total "$OUT/_tmp_m")" -v b="$nm" 'BEGIN{printf "%.4f", a/b}')
    fi
    if [ "${nf:-0}" -gt 0 ]; then
      kmc_tools -t"$THREADS" simple "$FSK" "$KDIR/${S}.${SFX}" intersect "$OUT/_tmp_f" > /dev/null
      ff=$(awk -v a="$(kmer_total "$OUT/_tmp_f")" -v b="$nf" 'BEGIN{printf "%.4f", a/b}')
    fi
    echo -e "${S}\t${fm}\t${ff}" >> "$OUT/kmer_fractions.tsv"
  done
  rm -f "$OUT"/_tmp_*
  ARGS+=(--kmer-fractions "$OUT/kmer_fractions.tsv")
else
  log "k-mer evidence skipped (needs stage 8 or SEXING_MSK_DB/SEXING_FSK_DB, and kmc_tools)"
fi

python "$BIN/assign_sex.py" "${ARGS[@]}" | tee "$OUT/sex_assignment.log"
log "Done. Outputs in $OUT"
