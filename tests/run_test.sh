#!/usr/bin/env bash
# End-to-end test of SDRscan on simulated data.
#
# Simulates an XY dataset with a known sex-determining region (SDR), runs the
# pipeline and checks that each analysis finds what was planted.
#
# Usage (from the repository root, with the sdrscan environment active):
#   tests/run_test.sh NAME PLOIDY MALES FEMALES [UNKNOWN] [COVERAGE]
# e.g.
#   tests/run_test.sh diploid_2v2    2 2 2
#   tests/run_test.sh tetraploid_6v6 4 6 6 2 20
#   tests/run_test.sh pilot_1v1      2 1 1
#
# Environment: THREADS (default 8), JOBS (parallel freebayes chunks, default 4),
# KMER_MEM (GB for KMC, default 3), RUN_FINDZX=1 to include stage 4
# (needs the findZX conda env; slow), TEST_ROOT (default tests/output).
set -euo pipefail

NAME="${1:?usage: $0 NAME PLOIDY MALES FEMALES [UNKNOWN] [COVERAGE]}"
PLOIDY="$2"; MALES="$3"; FEMALES="$4"; UNKNOWN="${5:-0}"; COV="${6:-15}"
REPO="$(cd "$(dirname "$0")/.." && pwd)"
ROOT="${TEST_ROOT:-$REPO/tests/output}"
D="$ROOT/$NAME"
THREADS="${THREADS:-8}"
mkdir -p "$D"
cd "$REPO"

echo "### [$NAME] simulating: ploidy $PLOIDY, $MALES M, $FEMALES F, $UNKNOWN U, ${COV}x"
[ -s "$D/data/samples.tsv" ] || python tests/simulate_dataset.py --outdir "$D/data" --ploidy "$PLOIDY" \
  --males "$MALES" --females "$FEMALES" --unknown "$UNKNOWN" --coverage "$COV" --seed 7

PILOT=0
{ [ "$MALES" -lt 2 ] || [ "$FEMALES" -lt 2 ]; } && PILOT=1

cat > "$D/pipeline.env" <<EOF
SAMPLES=$D/data/samples.tsv
REF=$D/data/refs/female.fa
OUTDIR=$D/results
THREADS=$THREADS
SORT_MEM=1G
MAX_PARALLEL_JOBS=${JOBS:-4}
CALL_CHUNK_BP=250000
MARK_DUPLICATES=0
PLOIDY=$PLOIDY
MIN_MAPQ=20
MIN_BASEQ=20
MIN_ALT_COUNT=2
MIN_QUAL=20
MIN_DP=3
ALLOW_SINGLE_PER_SEX=$PILOT
CHROM_REGEX='chr_[0-9]+\$'
FINDZX_RUNS=("XY|$D/data/refs/female.fa")
FINDZX_CORES=$THREADS
FINDZX_CONDA_ENV=findZX
FINDZX_BACKGROUND=0
FINDZX_WINDOWS="[20000, 50000]"
FINDZX_CONDA_PREFIX=$ROOT/findzx_conda
WINDOW=20000
STEP=10000
SNP_RUNS=(
  "relaxed|--min-called-males 0.75 --min-called-females 0.75 --min-alt-males 0.75 --min-alt-females 0.75 --max-alt-opposite 0.2"
  "strict|--min-called-males 0.75 --min-called-females 0.75 --min-alt-males 1.0 --min-alt-females 1.0 --max-alt-opposite 0"
)
RELAXED_RUN=relaxed
STRICT_RUN=strict
SNP_DOSAGE_CLASSES=auto
SNP_DOSAGE_BASE=relaxed
CANDIDATES=("chr_03:300000-700000")
FOCAL_SEX=M
PERMUTATIONS=2000
STRONG_ALT=0.75
STRONG_MAX_OPPOSITE=0.2
K=31
KMER_MIN_COUNT=2
KMER_MEM_GB=${KMER_MEM:-3}
KMER_SUBSAMPLE=200000
KMER_WINDOW=10000
KMER_STEP=5000
KMER_RUNS=("MSK|auto" "FSK|auto")
KMER_REFS=("F=$D/data/refs/female.fa" "M=$D/data/refs/male.fa")
KMER_NAME_NORMALISE=''
DOTPLOT_PAIRS=$D/data/pairs.tsv
DOTPLOT_MIN_IDENTITY=95
DOTPLOT_MIN_LENGTH=1000
DOTPLOT_CHROMS=("female_vs_male|chr_03|chr_03|chr3")
DEPTH_WINDOW=20000
DEPTH_LOG2_THRESHOLD=auto
HET_DIFF_THRESHOLD=0.1
SEX_SYSTEM=XY
SEXING_REGIONS=
SEXING_MARKERS=
SEXING_MSK_DB=
SEXING_FSK_DB=
EOF

STAGES="1 2 3 5 6 7 8 9 10 11 12 13"
[ "${RUN_FINDZX:-0}" = "1" ] && STAGES="1 2 3 4 5 6 7 8 9 10 11 12 13"
export CONFIG="$D/pipeline.env"
echo "### [$NAME] running stages: $STAGES"
start=$(date +%s)
# shellcheck disable=SC2086
./run_pipeline.sh $STAGES 2>&1 | tee "$D/pipeline.log"
echo "### [$NAME] pipeline finished in $(( $(date +%s) - start )) s"

python tests/check_results.py --results "$D/results" --data "$D/data" --pilot "$PILOT" \
  --males "$MALES" --females "$FEMALES" | tee "$D/checks.txt"
