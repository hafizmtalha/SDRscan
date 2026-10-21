#!/usr/bin/env bash
# Stage 4: findZX (read depth + heterozygosity between sexes), using the
# findZX submodule in external/findZX.
#
# For every FINDZX_RUNS entry (SYSTEM|REFERENCE) this writes a findZX units
# file and config under $OUTDIR/04_findzx/<run>/ and runs the findZX Snakefile.
#   XY        males heterogametic       (pair with a female reference)
#   ZW        females heterogametic     (pair with a male reference)
#   SHUFFLED  XY with sex labels randomly reassigned - negative control
#
# findZX needs its own environment (Snakemake 6.x); set FINDZX_CONDA_ENV to
# its name, or activate it before running this script.
#
# Usage: scripts/04_findzx.sh [SYSTEM|REFERENCE ...]   (default: FINDZX_RUNS)
source "$(dirname "$0")/common.sh"

FINDZX_DIR="${PIPELINE_DIR}/external/findZX"
if [ ! -f "$FINDZX_DIR/workflow/findZX" ]; then
  echo "findZX submodule missing. Run: git submodule update --init" >&2
  exit 1
fi

RUNS=("$@")
[ ${#RUNS[@]} -gt 0 ] || RUNS=("${FINDZX_RUNS[@]}")
OUT="$(mkdir -p "${OUTDIR}/04_findzx" && cd "${OUTDIR}/04_findzx" && pwd)"
SAMPLES_ABS="$(cd "$(dirname "$SAMPLES")" && pwd)/$(basename "$SAMPLES")"

if [ -n "${FINDZX_CONDA_ENV:-}" ]; then
  SNAKE=(conda run --no-capture-output -n "$FINDZX_CONDA_ENV" snakemake)
else
  SNAKE=(snakemake)
fi

for entry in "${RUNS[@]}"; do
  SYSTEM="${entry%%|*}"
  FREF="${entry#*|}"
  RUN_NAME="${SYSTEM}_$(basename "${FREF%.*}")"
  d="$OUT/$RUN_NAME"
  mkdir -p "$d"
  log "findZX run $RUN_NAME"

  case "$SYSTEM" in
    XY|ZW)    python "$BIN/../findzx/make_units.py" --samples "$SAMPLES_ABS" --system "$SYSTEM" --out "$d/units.tsv" ;;
    SHUFFLED) python "$BIN/../findzx/make_units.py" --samples "$SAMPLES_ABS" --system XY --shuffle --seed 7 --out "$d/units.tsv" ;;
    *) echo "Unknown findZX system: $SYSTEM" >&2; exit 1 ;;
  esac

  [ -f "${FREF}.fai" ] || samtools faidx "$FREF"
  cut -f1 "${FREF}.fai" | grep -E "$CHROM_REGEX" > "$d/main_chromosomes.txt" || true
  CHR_FILE="$d/main_chromosomes.txt"
  [ -s "$CHR_FILE" ] || CHR_FILE=None

  sed -e "s#^run_name:.*#run_name: ${RUN_NAME}#" \
      -e "s#^units:.*#units: ${d}/units.tsv#" \
      -e "s#^ref_genome:.*#ref_genome: ${FREF}#" \
      -e "s#^chr_file:.*#chr_file: ${CHR_FILE}#" \
      -e "s#^threads_max:.*#threads_max: ${FINDZX_CORES}#" \
      -e "s#^window_sizes:.*#window_sizes: ${FINDZX_WINDOWS}#" \
      "$PIPELINE_DIR/findzx/config.example.yml" > "$d/config.yml"

  n_units=$(($(wc -l < "$d/units.tsv") - 1))
  log "$RUN_NAME: $n_units samples in units.tsv"

  (
    cd "$FINDZX_DIR"
    "${SNAKE[@]}" -s workflow/findZX --configfile "$d/config.yml" --unlock > /dev/null 2>&1 || true
    CMD=("${SNAKE[@]}" -s workflow/findZX --configfile "$d/config.yml" --cores "$FINDZX_CORES"
         --use-conda --conda-prefix "$FINDZX_CONDA_PREFIX" --conda-frontend mamba
         --rerun-incomplete -k -p)
    if [ "${FINDZX_BACKGROUND:-0}" = "1" ]; then
      nohup "${CMD[@]}" > "$d/findzx.log" 2>&1 &
      log "$RUN_NAME running in background (PID $!), log: $d/findzx.log"
    else
      "${CMD[@]}" > "$d/findzx.log" 2>&1
    fi
  )
  ln -sfn "$FINDZX_DIR/results/$RUN_NAME" "$d/results"
done

log "findZX outputs: $OUT/<run>/results/output/no_synteny/{plots,tables}"
log "Check that the plots show all samples before using them."
