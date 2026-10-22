#!/usr/bin/env bash
# Stage 1: trim, QC, align and add read groups for every sample in the sheet.
#   fastp -> FastQC (raw + clean) -> bwa mem -> samtools sort/index
#   -> GATK AddOrReplaceReadGroups [-> MarkDuplicatesSpark] -> MultiQC
#
# Usage: scripts/01_qc_align.sh [sample ...]     (default: all samples)
source "$(dirname "$0")/common.sh"

OUT="${OUTDIR}/01_align"
mkdir -p "$OUT"/{fastp,fastqc_raw,fastqc_clean,bam,logs}

[ -f "${REF}.fai" ] || samtools faidx "$REF"
[ -f "${REF}.bwt" ] || bwa index "$REF"
DICT="${REF%.*}.dict"
[ -f "$DICT" ] || gatk CreateSequenceDictionary -R "$REF" -O "$DICT"

SAMPLE_LIST=("$@")
[ ${#SAMPLE_LIST[@]} -gt 0 ] || mapfile -t SAMPLE_LIST < <(samples_of)

for S in "${SAMPLE_LIST[@]}"; do
  R1=$(sample_field "$S" fq1)
  R2=$(sample_field "$S" fq2)
  FINAL="$OUT/bam/${S}.RG.bam"
  [ "$MARK_DUPLICATES" = "1" ] && FINAL="$OUT/bam/${S}.RG.dedup.bam"
  if [ -s "$FINAL" ]; then
    log "$S: $FINAL exists, skipping"
    continue
  fi
  log "$S: fastp"
  fastp -i "$R1" -I "$R2" \
    -o "$OUT/fastp/${S}_clean_R1.fastq.gz" -O "$OUT/fastp/${S}_clean_R2.fastq.gz" \
    -h "$OUT/fastp/${S}.fastp.html" -j "$OUT/fastp/${S}.fastp.json" \
    --thread "$THREADS" 2> "$OUT/logs/${S}.fastp.log"

  log "$S: FastQC"
  fastqc --quiet --threads 2 --outdir "$OUT/fastqc_raw" "$R1" "$R2"
  fastqc --quiet --threads 2 --outdir "$OUT/fastqc_clean" "$OUT/fastp/${S}_clean_R1.fastq.gz" "$OUT/fastp/${S}_clean_R2.fastq.gz"

  log "$S: bwa mem + sort"
  bwa mem -t "$THREADS" "$REF" "$OUT/fastp/${S}_clean_R1.fastq.gz" "$OUT/fastp/${S}_clean_R2.fastq.gz" \
      2> "$OUT/logs/${S}.bwa.log" \
    | samtools sort -@ "$THREADS" -m "$SORT_MEM" -o "$OUT/bam/${S}.sorted.bam" -
  samtools index "$OUT/bam/${S}.sorted.bam"

  log "$S: add read groups"
  gatk AddOrReplaceReadGroups \
    -I "$OUT/bam/${S}.sorted.bam" -O "$OUT/bam/${S}.RG.bam" \
    -RGID "$S" -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM "$S" \
    2> "$OUT/logs/${S}.readgroups.log"
  samtools index "$OUT/bam/${S}.RG.bam"

  if [ "$MARK_DUPLICATES" = "1" ]; then
    log "$S: mark duplicates"
    mkdir -p "$OUT/tmp_${S}"
    gatk MarkDuplicatesSpark -I "$OUT/bam/${S}.RG.bam" -O "$OUT/bam/${S}.RG.dedup.bam" \
      --metrics-file "$OUT/bam/${S}.dedup.metrics" --tmp-dir "$OUT/tmp_${S}" \
      2> "$OUT/logs/${S}.markdup.log"
    samtools index "$OUT/bam/${S}.RG.dedup.bam"
    rm -rf "$OUT/tmp_${S}"
  fi
  rm -f "$OUT/bam/${S}.sorted.bam" "$OUT/bam/${S}.sorted.bam.bai"
done

multiqc -q -f "$OUT" -o "$OUT/multiqc"
log "Done. BAMs in $OUT/bam"
