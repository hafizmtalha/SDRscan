#!/usr/bin/env bash
# Run findZX (https://github.com/hsigeman/findZX) in the background with the
# settings that worked for us. Run from inside a findZX checkout, with the
# findZX conda environment active.
#
# Usage: findzx/run_findzx.sh config/my_config.yml [cores] [log]
#
# Notes from our runs:
#  - Pin the conda prefix: re-creating the per-rule environments with newer
#    package versions broke several rules.
#  - mamba as the conda frontend was much faster and more reliable.
#  - Snakemake only re-runs missing outputs, so after fixing an error just run
#    this again (with --unlock first if a run was killed).
#  - Before trusting the plots, check the sample count in the outputs matches
#    units.tsv; one of our runs silently used an older units file.
set -euo pipefail

CFG="${1:?usage: $0 config.yml [cores] [log]}"
CORES="${2:-16}"
LOG="${3:-findzx_$(basename "${CFG%.*}").log}"
CONDA_PREFIX_DIR="${FINDZX_CONDA_PREFIX:-$HOME/snakemake_conda_store}"

echo "Units file used: $(grep -E '^units:' "$CFG")"
snakemake -s workflow/findZX --configfile "$CFG" --unlock > /dev/null 2>&1 || true

nohup snakemake -s workflow/findZX \
  --configfile "$CFG" \
  --cores "$CORES" \
  --use-conda --conda-prefix "$CONDA_PREFIX_DIR" --conda-frontend mamba \
  --rerun-incomplete -k -p \
  > "$LOG" 2>&1 &

echo "findZX started (PID $!). Follow progress with: tail -f $LOG"
