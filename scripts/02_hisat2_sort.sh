#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/f/RNAseq_SRP043694"
TRIMMED="${PROJECT}/data/trimmed_fastq"
ALIGN="${PROJECT}/data/aligned_bam"
LOGS="${PROJECT}/logs/hisat2"

INDEX="/mnt/e/bioinformatics_work/digital_home/rna_project/reference/Homo_sapiens_index"
THREADS=4

mkdir -p "$ALIGN" "$LOGS"

shopt -s nullglob

for r1 in "$TRIMMED"/*_1.trimmed.fastq.gz; do
    sample="$(basename "$r1" _1.trimmed.fastq.gz)"
    r2="${TRIMMED}/${sample}_2.trimmed.fastq.gz"

    sam="${ALIGN}/${sample}.sam"
    bam="${ALIGN}/${sample}.sorted.bam"
    log="${LOGS}/${sample}.log"

    if [[ ! -s "$r2" ]]; then
        echo "ERROR: missing mate for $sample" | tee -a "$log"
        exit 1
    fi

    if [[ -s "$bam" ]]; then
        echo "SKIP: $sample" | tee -a "$log"
        continue
    fi

    echo "RUN: $sample" | tee "$log"

    hisat2 \
      -x "$INDEX" \
      -1 "$r1" \
      -2 "$r2" \
      -p "$THREADS" \
      --summary-file "${ALIGN}/${sample}.hisat2.summary.txt" \
      2>> "$log" \
    | samtools sort -@ "$THREADS" -o "$bam" -

    samtools index "$bam"

    echo "DONE: $sample" | tee -a "$log"
done
