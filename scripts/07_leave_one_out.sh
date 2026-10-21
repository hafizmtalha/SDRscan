#!/usr/bin/env bash
# Robustness check: re-run the sex-specific SNP window analysis dropping one
# sample of the focal sex at a time, and report whether a candidate region
# still appears among the merged candidate intervals.
#
# Usage:
#   scripts/07_leave_one_out.sh M CHROM START END [classifier options]
# e.g.
#   scripts/07_leave_one_out.sh M chr_07 8200000 9000000 \
#     --min-called-males 4 --min-called-females 4 --min-alt-males 4 --min-alt-females 4 --max-alt-opposite 1
source "$(dirname "$0")/common.sh"

SEX="${1:?usage: $0 M|F CHROM START END [classifier options]}"
CHROM="$2"; START="$3"; END="$4"
shift 4
IN="${OUTDIR}/03_merged"
OUT="${OUTDIR}/07_leave_one_out"
mkdir -p "$OUT"
[ "$SEX" = "M" ] && FOCAL="male" || FOCAL="female"
[ "$SEX" = "M" ] && OTHER="female" || OTHER="male"

mapfile -t FOCAL_SAMPLES < "$IN/${FOCAL}s.txt"
echo -e "dropped\t${FOCAL}_sites\tregion_found\tmatching_intervals" > "$OUT/summary.tsv"

for drop in "${FOCAL_SAMPLES[@]}"; do
  d="$OUT/drop_${drop}"
  mkdir -p "$d"
  grep -vx "$drop" "$IN/${FOCAL}s.txt" > "$d/${FOCAL}s.txt"
  cp "$IN/${OTHER}s.txt" "$d/${OTHER}s.txt"

  python "$BIN/classify_sex_specific_snps.py" \
    --genotypes "$IN/merged.genotypes.tsv" --vcf-samples "$IN/vcf_samples.txt" \
    --males "$d/males.txt" --females "$d/females.txt" \
    --out "$d/sex_specific_sites.tsv" "$@" > "$d/classifier.log"

  awk -v c="${FOCAL}_specific" 'BEGIN{OFS="\t"} NR>1 && $3==c {print $1,$2-1,$2}' \
    "$d/sex_specific_sites.tsv" > "$d/sites.bed"
  bedtools makewindows -g "${REF}.fai" -w "$WINDOW" -s "$STEP" \
    | bedtools intersect -a - -b "$d/sites.bed" -c > "$d/windows.bed"
  {
    echo -e "SEX\tCHROM\tSTART\tEND\tCOUNT"
    awk -v OFS="\t" -v s="${FOCAL}_specific" '$4>=1{print s,$1,$2,$3,$4}' "$d/windows.bed" \
      | sort -k5,5nr | head -n 50
  } > "$d/top_candidate_windows.tsv"
  python "$BIN/merge_top_candidate_windows.py" --input "$d/top_candidate_windows.tsv" \
    --output "$d/merged_candidate_intervals.tsv" --merge-gap "$STEP" --min-count 2 > /dev/null

  hits=$(awk -v c="$CHROM" -v s="$START" -v e="$END" 'NR>1 && $2==c && $3<e && $4>s' \
           "$d/merged_candidate_intervals.tsv")
  found=no; [ -n "$hits" ] && found=yes
  echo -e "${drop}\t$(wc -l < "$d/sites.bed")\t${found}\t$(echo "$hits" | awk '{printf "%s:%s-%s(max=%s);", $2,$3,$4,$6}')" \
    >> "$OUT/summary.tsv"
done

column -t "$OUT/summary.tsv"
