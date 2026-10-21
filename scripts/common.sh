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

# ---- Sample sizes ----
# The pipeline works with any number of samples, but comparing sexes needs at
# least MIN_PER_SEX (2) of each, and some analyses need more. See README,
# "How many samples do I need?".
MIN_PER_SEX=2

n_of() { samples_of "$1" | wc -l | tr -d ' '; }

# Number of distinct ways to label n_focal of n_total samples as the focal sex
n_labellings() {
  awk -v n="$1" -v k="$2" 'BEGIN{c=1; for(i=1;i<=k;i++) c=c*(n-k+i)/i; printf "%.0f", c}'
}

# Sex values allowed in the sample sheet: M, F, or U (unknown; such samples
# are left out of every male-vs-female comparison and can be sexed in stage 13)
bad_sex=$(awk -F'\t' -v c="$(_col sex)" 'NR>1 && $1!="" && $c!="M" && $c!="F" && $c!="U" {print $1"="$c}' "$SAMPLES")
[ -z "$bad_sex" ] || die "sex must be M, F or U in $SAMPLES: $bad_sex"

# A pilot design has at least one but fewer than MIN_PER_SEX of a sex
is_pilot() {
  local m f
  m=$(n_of M); f=$(n_of F)
  [ "$m" -ge 1 ] && [ "$f" -ge 1 ] && { [ "$m" -lt "$MIN_PER_SEX" ] || [ "$f" -lt "$MIN_PER_SEX" ]; }
}

# Stop unless both sexes have at least MIN_PER_SEX samples. Stages that are
# meaningful for a 1 + 1 pilot call it with "pilot_ok" and then only need
# one of each, provided ALLOW_SINGLE_PER_SEX=1.
require_both_sexes() {
  local m f
  m=$(n_of M); f=$(n_of F)
  [ "$m" -ge "$MIN_PER_SEX" ] && [ "$f" -ge "$MIN_PER_SEX" ] && return 0
  if [ "${1:-}" = "pilot_ok" ] && [ "${ALLOW_SINGLE_PER_SEX:-0}" = "1" ] && is_pilot; then
    log "PILOT MODE: $m male / $f female. Results show differences between individuals, not replicated sex differences."
    return 0
  fi
  if is_pilot && [ "${1:-}" = "pilot_ok" ]; then
    die "need at least $MIN_PER_SEX males and $MIN_PER_SEX females (found $m M, $f F). Set ALLOW_SINGLE_PER_SEX=1 for a pilot run."
  fi
  die "need at least $MIN_PER_SEX males and $MIN_PER_SEX females in $SAMPLES (found $m M, $f F); this stage has no pilot mode"
}

# Print which analyses the current design supports
sample_size_report() {
  local m f n u lab minp
  m=$(n_of M); f=$(n_of F); u=$(n_of U); n=$((m + f))
  log "Samples: $m male, $f female, $u unknown sex (ploidy $(samples_of | while read -r s; do sample_ploidy "$s"; done | sort -un | tr '\n' ' '))"
  [ "$u" -gt 0 ] && log "  unknown-sex samples are left out of sex comparisons; assign them with stage 13"
  if [ "$m" -lt 1 ] || [ "$f" -lt 1 ]; then
    log "  ! no samples of one sex: only stages 1-3 and 10, and stage 13 with markers from an earlier project"
    return 0
  fi
  if is_pilot; then
    log "  ! PILOT design (< $MIN_PER_SEX of a sex): with ALLOW_SINGLE_PER_SEX=1 stages 1-4, 10 and 12 run;"
    log "    SNP density, k-mers, permutation tests, PCA and leave-one-out need >= $MIN_PER_SEX per sex"
    return 0
  fi
  lab=$(n_labellings "$n" "$m")
  minp=$(awk -v l="$lab" 'BEGIN{printf "%.3g", 1/l}')
  log "  sex-label permutation: $lab distinct labellings, smallest possible p = $minp"
  awk -v p="$minp" 'BEGIN{exit !(p >= 0.05)}' && \
    log "  ! permutation tests cannot reach p < 0.05 with this design (need more samples)"
  for sex in M F; do
    [ "$(n_of "$sex")" -lt 3 ] && log "  ! leave-one-out needs >= 3 of the focal sex (have $(n_of "$sex") $sex)"
  done
  [ "$n" -lt 6 ] && log "  ! PCA with < 6 samples is descriptive only"
  return 0
}

# Resolve a "how many samples" setting: a count (7), a fraction (0.8),
# "all" or "all-N", against n samples. Result is clamped to 1..n.
resolve_count() {
  local spec="$1" n="$2" v
  case "$spec" in
    all)     v=$n ;;
    all-*)   v=$((n - ${spec#all-})) ;;
    *.*)     v=$(awk -v f="$spec" -v n="$n" 'BEGIN{v=f*n; c=int(v); if (c<v) c++; print c}') ;;
    *)       v=$spec ;;
  esac
  [[ "$v" =~ ^-?[0-9]+$ ]] || die "invalid sample count '$spec'"
  [ "$v" -lt 1 ] && v=1
  [ "$v" -gt "$n" ] && v=$n
  echo "$v"
}

# Path of the ploidy table written by stage 3; later stages pass it to bin/ tools
PLOIDY_TABLE="${OUTDIR}/03_merged/sample_ploidy.tsv"
require_ploidy_table() {
  [ -s "$PLOIDY_TABLE" ] || die "$PLOIDY_TABLE missing; run stage 3 first"
}
