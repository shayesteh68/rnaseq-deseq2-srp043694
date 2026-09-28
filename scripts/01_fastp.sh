#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/f/RNAseq_SRP043694"
RAW="${PROJECT}/raw_fastq"
TRIMMED="${PROJECT}/data/trimmed_fastq"
FASTP_QC="${PROJECT}/qc/fastp"
LOGS="${PROJECT}/logs"

mkdir -p "$TRIMMED" "$FASTP_QC" "$LOGS"

shopt -s nullglob

for r1 in "$RAW"/*_1.fastq.gz; do
    sample="$(basename "$r1" _1.fastq.gz)"
    r2="${RAW}/${sample}_2.fastq.gz"

    out1="${TRIMMED}/${sample}_1.trimmed.fastq.gz"
    out2="${TRIMMED}/${sample}_2.trimmed.fastq.gz"
    html="${FASTP_QC}/${sample}.html"
    json="${FASTP_QC}/${sample}.json"
    log="${LOGS}/${sample}.fastp.log"

    if [[ ! -s "$r2" ]]; then
        echo "ERROR: missing mate for $sample" | tee -a "$log"
        exit 1
    fi

    if [[ -s "$out1" && -s "$out2" && -s "$json" ]]; then
        echo "SKIP: $sample"
        continue
    fi

    echo "RUN: $sample" | tee "$log"

    fastp \
      --in1 "$r1" \
      --in2 "$r2" \
      --out1 "$out1" \
      --out2 "$out2" \
      --html "$html" \
      --json "$json" \
      --detect_adapter_for_pe \
      --qualified_quality_phred 20 \
      --length_required 35 \
      --thread 4 \
      >> "$log" 2>&1

    echo "DONE: $sample" | tee -a "$log"
done
