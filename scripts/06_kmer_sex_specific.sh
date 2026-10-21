#!/usr/bin/env bash
# Reference-free sex-specific k-mers, mapped back to one or more references.
#
#   1. Count k-mers per sample (KMC, k=K, -ci KMER_MIN_COUNT)
#   2. Focal-sex set: k-mers present in at least MIN_PRESENT focal samples
#   3. Subtract the union of all opposite-sex k-mers
#      -> MSK (male-specific k-mers) or FSK (female-specific k-mers)
#   4. Dump, subsample KMER_SUBSAMPLE k-mers, map to each reference with bwa mem
#   5. Keep MAPQ >= 20, count k-mers per window and per chromosome
#
# Step 2 counts, for every k-mer, how many focal samples contain it (each
# sample database is set to count 1 and the databases are unioned with summed
# counters). For MIN_PRESENT = n-1 this gives the same set as unioning all
# leave-one-out intersections, which is how the original runs built it.
#
# Usage:
#   scripts/06_kmer_sex_specific.sh MSK|FSK MIN_PRESENT TAG=ref.fa [TAG=ref.fa ...]
# e.g.
#   scripts/06_kmer_sex_specific.sh MSK 7 REF1=/refs/hapA.fa REF2=/refs/hapB.fa
source "$(dirname "$0")/common.sh"

MODE="${1:?usage: $0 MSK|FSK MIN_PRESENT TAG=ref.fa ...}"
MIN_PRESENT="${2:?MIN_PRESENT required}"
shift 2
[ $# -gt 0 ] || { echo "Give at least one TAG=ref.fa" >&2; exit 1; }
case "$MODE" in
  MSK) FOCAL=M; OTHER=F ;;
  FSK) FOCAL=F; OTHER=M ;;
  *) echo "MODE must be MSK or FSK" >&2; exit 1 ;;
esac

OUT="${OUTDIR}/06_kmer"
KDIR="$OUT/kmc"; SETDIR="$OUT/sets/${MODE}"; TMP="$OUT/tmp"
mkdir -p "$KDIR" "$SETDIR" "$TMP" "$OUT/${MODE}/fa" "$OUT/${MODE}/map" "$OUT/${MODE}/windows"
SFX="k${K}.ci${KMER_MIN_COUNT}"

# 1. Per-sample k-mer databases (shared between MSK and FSK runs)
for S in $(samples_of); do
  if [ -f "$KDIR/${S}.${SFX}.kmc_pre" ]; then continue; fi
  log "KMC: $S"
  printf '%s\n%s\n' "$(sample_field "$S" 3)" "$(sample_field "$S" 4)" > "$TMP/${S}.files"
  kmc -k"$K" -ci"$KMER_MIN_COUNT" -t"$THREADS" @"$TMP/${S}.files" "$KDIR/${S}.${SFX}" "$TMP" > "$TMP/${S}.kmc.log"
done

# Union a list of databases, summing counters; result written to $1
union_all() {
  local out="$1"; shift
  local acc="$1"; shift
  local i=0
  for db in "$@"; do
    i=$((i + 1))
    kmc_tools -t"$THREADS" simple "$acc" "$db" union "$SETDIR/_u${i}" -ocsum
    acc="$SETDIR/_u${i}"
  done
  cp "${acc}.kmc_pre" "${out}.kmc_pre"; cp "${acc}.kmc_suf" "${out}.kmc_suf"
}

# 2. Focal-sex presence set
mapfile -t FOC < <(samples_of "$FOCAL")
mapfile -t OPP < <(samples_of "$OTHER")
log "${MODE}: ${#FOC[@]} focal samples, require presence in >= ${MIN_PRESENT}"
PRES=()
for S in "${FOC[@]}"; do
  kmc_tools -t"$THREADS" transform "$KDIR/${S}.${SFX}" set_counts 1 "$SETDIR/${S}.presence"
  PRES+=("$SETDIR/${S}.presence")
done
union_all "$SETDIR/focal.count" "${PRES[@]}"
kmc_tools -t"$THREADS" transform "$SETDIR/focal.count" -ci"$MIN_PRESENT" reduce "$SETDIR/focal.present_ge${MIN_PRESENT}"

# 3. Opposite-sex union and subtraction
OPP_DBS=()
for S in "${OPP[@]}"; do OPP_DBS+=("$KDIR/${S}.${SFX}"); done
union_all "$SETDIR/opposite.union" "${OPP_DBS[@]}"
kmc_tools -t"$THREADS" simple "$SETDIR/focal.present_ge${MIN_PRESENT}" "$SETDIR/opposite.union" \
  kmers_subtract "$SETDIR/${MODE}"
rm -f "$SETDIR"/_u* "$SETDIR"/*.presence.*
kmc_tools info "$SETDIR/${MODE}" | tee "$OUT/${MODE}/${MODE}.info.txt"

# 4. Dump + subsample
FA="$OUT/${MODE}/fa"
kmc_dump "$SETDIR/${MODE}" "$FA/${MODE}.txt"
awk '{print ">kmer_" NR "\n" $1}' "$FA/${MODE}.txt" > "$FA/${MODE}.fasta"
seqtk sample -s100 "$FA/${MODE}.fasta" "$KMER_SUBSAMPLE" > "$FA/${MODE}.sub.fasta"

# 5. Map to each reference and count per window
for spec in "$@"; do
  TAGNAME="${spec%%=*}"; R="${spec#*=}"
  M="$OUT/${MODE}/map"; WD="$OUT/${MODE}/windows"
  log "${MODE}: mapping to ${TAGNAME}"
  [ -f "${R}.bwt" ] || bwa index "$R"
  [ -f "${R}.fai" ] || samtools faidx "$R"
  bwa mem -t "$THREADS" "$R" "$FA/${MODE}.sub.fasta" 2> "$M/${TAGNAME}.bwa.log" \
    | samtools view -b -q 20 - | samtools sort -@ "$THREADS" -o "$M/${MODE}_${TAGNAME}.mapq20.bam"
  samtools index "$M/${MODE}_${TAGNAME}.mapq20.bam"
  bedtools bamtobed -i "$M/${MODE}_${TAGNAME}.mapq20.bam" > "$M/${MODE}_${TAGNAME}.mapq20.bed"

  bedtools makewindows -g "${R}.fai" -w "$KMER_WINDOW" -s "$KMER_STEP" \
    | awk -v re="$CHROM_REGEX" '$1 ~ re' \
    | bedtools intersect -a - -b "$M/${MODE}_${TAGNAME}.mapq20.bed" -c > "$WD/genome.${MODE}_density.${TAGNAME}.bed"
  awk '$4>0' "$WD/genome.${MODE}_density.${TAGNAME}.bed" | sort -k4,4nr | head -n 200 \
    > "$WD/genome.${MODE}_density.${TAGNAME}.top200.tsv"
  # Candidate peaks: merge the top windows into intervals
  sort -k1,1 -k2,2n "$WD/genome.${MODE}_density.${TAGNAME}.top200.tsv" \
    | bedtools merge -i - -c 4 -o max > "$WD/kmer_${MODE}_candidates.${TAGNAME}.bed"
  cut -f1 "$M/${MODE}_${TAGNAME}.mapq20.bed" | sort | uniq -c | sort -nr > "$WD/${MODE}_${TAGNAME}.chrom_counts.tsv"
  log "${MODE} ${TAGNAME}: $(wc -l < "$M/${MODE}_${TAGNAME}.mapq20.bed") k-mers mapped at MAPQ>=20"
done
log "Done. Outputs in $OUT/${MODE}"
