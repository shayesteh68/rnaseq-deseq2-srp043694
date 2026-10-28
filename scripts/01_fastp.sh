#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
RAW="${RAW_FASTQ_DIR:-$PROJECT_ROOT/raw_fastq}"
TRIMMED="${TRIMMED_FASTQ_DIR:-$PROJECT_ROOT/data/trimmed_fastq}"
FASTP_QC="${FASTP_QC_DIR:-$PROJECT_ROOT/qc/fastp}"
LOGS="${LOG_DIR:-$PROJECT_ROOT/logs}"
THREADS="${THREADS:-2}"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: THREADS must be a positive integer." >&2; exit 2; }
for tool in fastp gzip; do command -v "$tool" >/dev/null || { echo "ERROR: $tool is required." >&2; exit 2; }; done
mkdir -p "$TRIMMED" "$FASTP_QC" "$LOGS"
working=""
cleanup() {
  if [[ -n "$working" && -d "$working" ]]; then
    rm -f -- "$working/r1.fastq.gz" "$working/r2.fastq.gz" "$working/fastp.html" "$working/fastp.json"
    rmdir -- "$working" 2>/dev/null || true
  fi
}
trap cleanup EXIT
shopt -s nullglob
r1_files=("$RAW"/*_1.fastq.gz)
(( ${#r1_files[@]} > 0 )) || { echo "ERROR: no paired-end R1 files in $RAW" >&2; exit 1; }
for r1 in "${r1_files[@]}"; do
  sample="$(basename "$r1" _1.fastq.gz)"
  r2="$RAW/${sample}_2.fastq.gz"
  out1="$TRIMMED/${sample}_1.trimmed.fastq.gz"
  out2="$TRIMMED/${sample}_2.trimmed.fastq.gz"
  html="$FASTP_QC/$sample.html"
  json="$FASTP_QC/$sample.json"
  marker="$FASTP_QC/$sample.complete"
  log="$LOGS/$sample.fastp.log"
  [[ -s "$r1" && -s "$r2" ]] || { echo "ERROR: missing or empty mate for $sample" >&2; exit 1; }
  if [[ -e "$out1" || -e "$out2" || -e "$json" || -e "$html" || -e "$marker" ]]; then
    if [[ -f "$marker" && -s "$out1" && -s "$out2" && -s "$json" && -s "$html" ]] && gzip -t "$out1" "$out2"; then
      echo "SKIP: completed trimming for $sample"
      continue
    fi
    echo "ERROR: incomplete or unverified outputs for $sample; inspect them or use fresh output directories." >&2
    exit 1
  fi
  working="$(mktemp -d "$TRIMMED/.$sample.fastp.XXXXXX")"
  fastp --in1 "$r1" --in2 "$r2" --out1 "$working/r1.fastq.gz" --out2 "$working/r2.fastq.gz" --html "$working/fastp.html" --json "$working/fastp.json" --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 35 --thread "$THREADS" > "$log" 2>&1
  gzip -t "$working/r1.fastq.gz" "$working/r2.fastq.gz"
  [[ -s "$working/fastp.json" && -s "$working/fastp.html" ]] || { echo "ERROR: missing fastp reports for $sample" >&2; exit 1; }
  mv -- "$working/r1.fastq.gz" "$out1"
  mv -- "$working/r2.fastq.gz" "$out2"
  mv -- "$working/fastp.html" "$html"
  mv -- "$working/fastp.json" "$json"
  printf 'fastp finished and both gzip streams validated\n' > "$marker"
  cleanup
  working=""
  echo "DONE: $sample"
done
