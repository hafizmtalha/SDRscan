#!/usr/bin/env bash
# Stage 4: sex-specific SNP density.
#
#   classify sites (ALT present in most of one sex, absent/rare in the other)
#   -> BED -> sliding windows -> counts per window -> top windows
#   -> merged candidate intervals -> per-chromosome counts -> genome-wide plot
#
# Usage:
#   scripts/05_sex_specific_snp_windows.sh RUN_NAME [classifier options]
#
# Classifier options are passed to bin/classify_sex_specific_snps.py, e.g.
#   scripts/05_sex_specific_snp_windows.sh strict \
#     --min-called-males 7 --min-called-females 5 \
#     --min-alt-males 7 --min-alt-females 5 --max-alt-opposite 1
#   scripts/05_sex_specific_snp_windows.sh dosage1 ... --alt-dosage 1   # 0/0/0/1 only
#   scripts/05_sex_specific_snp_windows.sh dosage2 ... --alt-dosage 2   # 0/0/1/1 only
source "$(dirname "$0")/common.sh"

RUN="${1:?usage: $0 RUN_NAME [classifier options]}"
shift
IN="${OUTDIR}/03_merged"
OUT="${OUTDIR}/05_snp_windows/${RUN}"
mkdir -p "$OUT"
W="$WINDOW"; S="$STEP"
TAG="${W}_${S}"

echo "$*" > "$OUT/classifier_options.txt"
python "$BIN/classify_sex_specific_snps.py" \
  --genotypes "$IN/merged.genotypes.tsv" --vcf-samples "$IN/vcf_samples.txt" \
  --males "$IN/males.txt" --females "$IN/females.txt" \
  --out "$OUT/sex_specific_sites.tsv" "$@" | tee "$OUT/classifier.log"

for sex in male female; do
  awk -v c="${sex}_specific" 'BEGIN{OFS="\t"} NR>1 && $3==c {print $1,$2-1,$2}' \
    "$OUT/sex_specific_sites.tsv" > "$OUT/${sex}_specific_sites.bed"
done

bedtools makewindows -g "${REF}.fai" -w "$W" -s "$S" \
  | awk -v re="$CHROM_REGEX" '$1 ~ re' > "$OUT/windows_${TAG}.bed"
for sex in male female; do
  bedtools intersect -a "$OUT/windows_${TAG}.bed" -b "$OUT/${sex}_specific_sites.bed" -c \
    > "$OUT/${sex}_specific_windows_${TAG}.bed"
done

# Top 50 windows per sex
{
  echo -e "SEX\tCHROM\tSTART\tEND\tCOUNT"
  for sex in male female; do
    awk -v OFS="\t" -v s="${sex}_specific" '$4>=1{print s,$1,$2,$3,$4}' "$OUT/${sex}_specific_windows_${TAG}.bed" \
      | sort -k5,5nr -k2,2 -k3,3n | head -n 50
  done
} > "$OUT/top_candidate_windows.tsv"

# Merge overlapping top windows (bedtools) ...
for sex in male female; do
  awk -v s="${sex}_specific" 'NR>1 && $1==s {print $2"\t"$3"\t"$4"\t"$5}' "$OUT/top_candidate_windows.tsv" \
    | sort -k1,1 -k2,2n > "$OUT/${sex}_top_windows.bed"
  bedtools merge -i "$OUT/${sex}_top_windows.bed" -c 4 -o max > "$OUT/${sex}_top_windows.merged.bed"
done
# ... and the ranked interval table with window counts (MERGE_GAP = one step)
python "$BIN/merge_top_candidate_windows.py" --input "$OUT/top_candidate_windows.tsv" \
  --output "$OUT/merged_candidate_intervals.tsv" --merge-gap "$S" --min-count 2

# Sex-specific sites per chromosome
for sex in male female; do
  cut -f1 "$OUT/${sex}_specific_sites.bed" | sort | uniq -c | sort -nr > "$OUT/${sex}_specific_chr_counts.tsv"
done

python "$BIN/plot_sex_specific_windows.py" \
  --male "$OUT/male_specific_windows_${TAG}.bed" --female "$OUT/female_specific_windows_${TAG}.bed" \
  --outprefix "$OUT/sex_specific_snp_windows" --top-windows "$OUT/top_windows_plot.tsv"

log "Done. Outputs in $OUT"
