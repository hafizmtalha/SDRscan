#!/usr/bin/env bash
# Pairwise overlap of k-mer candidate intervals between runs (MSK/FSK x
# references). Only runs on the same coordinate system can be compared, so
# give the runs in groups that share chromosome names, or normalise names
# first with --rename 'regex=replacement' (applied with sed -E to column 1).
#
# Usage:
#   scripts/07_compare_kmer_runs.sh OUTDIR [--rename 'hap[12]_=hap_'] NAME=file.bed NAME=file.bed ...
# e.g. (two haplotypes of the same assembly, names differ only by haplotype tag)
#   scripts/07_compare_kmer_runs.sh results/07_kmer_overlap --rename '_H[12]_chr_=_chr_' \
#     MSK_A1=.../kmer_MSK_candidates.A1.bed MSK_A2=.../kmer_MSK_candidates.A2.bed \
#     FSK_A1=.../kmer_FSK_candidates.A1.bed FSK_A2=.../kmer_FSK_candidates.A2.bed
set -euo pipefail

OUT="${1:?usage: $0 OUTDIR [--rename 'regex=repl'] NAME=file.bed ...}"
shift
RENAME=""
if [ "${1:-}" = "--rename" ]; then RENAME="$2"; shift 2; fi
mkdir -p "$OUT"

NAMES=()
for spec in "$@"; do
  name="${spec%%=*}"; file="${spec#*=}"
  if [ -n "$RENAME" ]; then
    sed -E "s/${RENAME%%=*}/${RENAME#*=}/" "$file" | sort -k1,1 -k2,2n > "$OUT/${name}.bed"
  else
    sort -k1,1 -k2,2n "$file" > "$OUT/${name}.bed"
  fi
  NAMES+=("$name")
done

: > "$OUT/overlap_summary.tsv"
for ((i = 0; i < ${#NAMES[@]}; i++)); do
  for ((j = i + 1; j < ${#NAMES[@]}; j++)); do
    a="${NAMES[$i]}"; b="${NAMES[$j]}"
    bedtools intersect -a "$OUT/$a.bed" -b "$OUT/$b.bed" -wa -wb > "$OUT/overlap_${a}__${b}.tsv"
    echo -e "$(wc -l < "$OUT/overlap_${a}__${b}.tsv")\t${a}\t${b}" >> "$OUT/overlap_summary.tsv"
  done
done
sort -k1,1nr "$OUT/overlap_summary.tsv" -o "$OUT/overlap_summary.tsv"
echo -e "overlaps\trun_a\trun_b"
cat "$OUT/overlap_summary.tsv"
