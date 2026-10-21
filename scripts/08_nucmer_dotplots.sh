#!/usr/bin/env bash
# Whole-genome alignments between reference assemblies (e.g. the two
# haplotypes of a female and a male assembly) with MUMmer4, for dotplots and
# chromosome correspondence.
#
#   (optional) keep only main chromosomes -> nucmer --mum
#   -> delta-filter (identity, length, 1-to-1 query) -> show-coords -> mummerplot
#
# The .coords files feed R/dotplot_genome.R and R/dotplot_chromosome.R, which
# make the publication versions of the plots.
#
# Usage:
#   scripts/08_nucmer_dotplots.sh pairs.tsv [OUTDIR]
# pairs.tsv: one line per comparison: <ref.fa> <query.fa> <name>
#
# Environment overrides: THREADS MIN_IDENTITY MIN_LENGTH MIN_CLUSTER MIN_MATCH CHROM_REGEX
set -euo pipefail

PAIRS="${1:?usage: $0 pairs.tsv [OUTDIR]}"
OUT="${2:-results/08_dotplots}"
THREADS="${THREADS:-16}"
MIN_IDENTITY="${MIN_IDENTITY:-95}"
MIN_LENGTH="${MIN_LENGTH:-1000}"
MIN_CLUSTER="${MIN_CLUSTER:-100}"
MIN_MATCH="${MIN_MATCH:-20}"
CHROM_REGEX="${CHROM_REGEX:-}"      # e.g. 'chr_?[0-9]+$' to drop unplaced scaffolds
mkdir -p "$OUT/fasta"

main_chroms() {
  local fa="$1" out
  out="$OUT/fasta/$(basename "${fa%.*}").main.fa"
  if [ -z "$CHROM_REGEX" ]; then echo "$fa"; return; fi
  if [ ! -s "$out" ]; then
    [ -f "${fa}.fai" ] || samtools faidx "$fa"
    cut -f1 "${fa}.fai" | grep -E "$CHROM_REGEX" | xargs samtools faidx "$fa" > "$out"
  fi
  echo "$out"
}

while read -r REF QRY NAME; do
  [ -z "${REF:-}" ] || [[ "$REF" == \#* ]] && continue
  d="$OUT/$NAME"
  mkdir -p "$d"
  R=$(main_chroms "$REF"); Q=$(main_chroms "$QRY")
  echo "== $NAME: $R vs $Q"

  nucmer --mum --mincluster="$MIN_CLUSTER" --minmatch="$MIN_MATCH" --threads="$THREADS" \
    --prefix="$d/${NAME}" "$R" "$Q"
  delta-filter -i "$MIN_IDENTITY" -l "$MIN_LENGTH" -q "$d/${NAME}.delta" \
    > "$d/${NAME}.filtered${MIN_IDENTITY}.delta"
  show-coords -T -r -c -l "$d/${NAME}.filtered${MIN_IDENTITY}.delta" \
    > "$d/${NAME}.filtered${MIN_IDENTITY}.coords"
  mummerplot --png --large --filter --fat \
    -R "$R" -Q "$Q" --title "$NAME" --prefix "$d/${NAME}" "$d/${NAME}.filtered${MIN_IDENTITY}.delta" \
    > /dev/null 2>&1 || echo "mummerplot failed for $NAME (gnuplot missing?) - use the R scripts instead"

  cat > "$d/summary.txt" <<EOF
pair        $NAME
reference   $REF
query       $QRY
identity    >= ${MIN_IDENTITY}%
min_length  ${MIN_LENGTH} bp
alignments  $(grep -c '^[0-9]' "$d/${NAME}.filtered${MIN_IDENTITY}.coords" || true)
EOF
done < "$PAIRS"
echo "Done. Results in $OUT"
