# Sourced by every stage script.
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${PIPELINE_DIR}/bin"
CONFIG="${CONFIG:-${PIPELINE_DIR}/config/pipeline.env}"

if [ ! -f "$CONFIG" ]; then
  echo "Config not found: $CONFIG (copy config/pipeline.env.example to config/pipeline.env)" >&2
  exit 1
fi
# shellcheck disable=SC1090
source "$CONFIG"

log() { echo "[$(date '+%F %T')] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

# ---- Sample sheet: sample, sex, fq1, fq2 [, ploidy] (columns found by header name) ----
_col() { awk -F'\t' -v name="$1" 'NR==1{for(i=1;i<=NF;i++) if($i==name){print i; exit}}' "$SAMPLES"; }

# Print sample names, optionally for one sex (M or F).
samples_of() {
  local sex="${1:-}" c
  c=$(_col sex)
  awk -F'\t' -v s="$sex" -v c="$c" 'NR>1 && $1!="" && (s=="" || $c==s) {print $1}' "$SAMPLES"
}

# Field of one sample by column name: sample_field <sample> <column name>
sample_field() {
  local c
  c=$(_col "$2")
  [ -n "$c" ] || return 0
  awk -F'\t' -v id="$1" -v c="$c" 'NR>1 && $1==id {print $c}' "$SAMPLES"
}

# ---- Ploidy ----
# PLOIDY must be set in the config; a "ploidy" column in the sample sheet
# overrides it per sample (mixed-ploidy populations). The value used for
# variant calling is recorded and every later stage reads that record, so the
# whole analysis uses the ploidy the genotypes were called with.
[[ "${PLOIDY:-}" =~ ^[1-9][0-9]*$ ]] || die "PLOIDY must be set to a positive integer in $CONFIG (2 diploid, 4 tetraploid, 6 hexaploid, ...)"

sample_ploidy() {
  local p
  p=$(sample_field "$1" ploidy)
  if [ -z "$p" ] || [ "$p" = "NA" ]; then p="$PLOIDY"; fi
  [[ "$p" =~ ^[1-9][0-9]*$ ]] || die "invalid ploidy '$p' for sample $1"
  echo "$p"
}

# Write sample<TAB>ploidy for every sample to $1
write_sample_ploidy() {
  local S
  : > "$1"
  for S in $(samples_of); do printf '%s\t%s\n' "$S" "$(sample_ploidy "$S")" >> "$1"; done
}

# Highest ploidy among the samples
max_ploidy() {
  local S m=0 p
  for S in $(samples_of); do p=$(sample_ploidy "$S"); [ "$p" -gt "$m" ] && m=$p; done
  echo "$m"
}

# Path of the ploidy table written by stage 3; later stages pass it to bin/ tools
PLOIDY_TABLE="${OUTDIR}/03_merged/sample_ploidy.tsv"
require_ploidy_table() {
  [ -s "$PLOIDY_TABLE" ] || die "$PLOIDY_TABLE missing; run stage 3 first"
}
