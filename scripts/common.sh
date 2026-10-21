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

# Print sample names from the sample sheet, optionally for one sex (M or F).
samples_of() {
  local sex="${1:-}"
  awk -v s="$sex" 'NR>1 && $1!="" && (s=="" || $2==s) {print $1}' "$SAMPLES"
}

# Column lookup from the sample sheet: sample_field <sample> <column number>
sample_field() {
  awk -v id="$1" -v c="$2" 'NR>1 && $1==id {print $c}' "$SAMPLES"
}
