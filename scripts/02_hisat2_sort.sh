#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
TRIMMED="${TRIMMED_FASTQ_DIR:-$PROJECT_ROOT/data/trimmed_fastq}"
ALIGN="${ALIGN_BAM_DIR:-$PROJECT_ROOT/data/aligned_bam}"
LOGS="${HISAT2_LOG_DIR:-$PROJECT_ROOT/logs/hisat2}"
HISAT2_INDEX="${HISAT2_INDEX:-}"
THREADS="${THREADS:-2}"
SORT_THREADS="${SORT_THREADS:-1}"
SORT_MEMORY="${SORT_MEMORY:-512M}"
for value in "$THREADS" "$SORT_THREADS"; do
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: thread counts must be positive integers." >&2; exit 2; }
done
for tool in hisat2 samtools; do command -v "$tool" >/dev/null || { echo "ERROR: $tool is required." >&2; exit 2; }; done
[[ -n "$HISAT2_INDEX" ]] || { echo "ERROR: set HISAT2_INDEX to the GRCh38 index prefix." >&2; exit 2; }
index_ready=false
for suffix in ht2 ht2l; do
  complete=true
  for part in {1..8}; do [[ -s "$HISAT2_INDEX.$part.$suffix" ]] || complete=false; done
  if "$complete"; then index_ready=true; break; fi
done
"$index_ready" || { echo "ERROR: expected all eight HISAT2 index files for $HISAT2_INDEX" >&2; exit 2; }
mkdir -p "$ALIGN" "$LOGS"
temporary=""
cleanup() { if [[ -n "$temporary" ]]; then rm -f -- "$temporary" "$temporary.bai"; fi; }
trap cleanup EXIT
shopt -s nullglob
r1_files=("$TRIMMED"/*_1.trimmed.fastq.gz)
(( ${#r1_files[@]} > 0 )) || { echo "ERROR: no trimmed R1 files in $TRIMMED" >&2; exit 1; }
for r1 in "${r1_files[@]}"; do
  sample="$(basename "$r1" _1.trimmed.fastq.gz)"
  r2="$TRIMMED/${sample}_2.trimmed.fastq.gz"
  bam="$ALIGN/$sample.sorted.bam"
  log="$LOGS/$sample.log"
  [[ -s "$r1" && -s "$r2" ]] || { echo "ERROR: missing trimmed mate for $sample" >&2; exit 1; }
  if [[ -e "$bam" ]]; then
    samtools quickcheck -v "$bam" || { echo "ERROR: invalid existing BAM: $bam" >&2; exit 1; }
    [[ -s "$bam.bai" ]] || samtools index "$bam"
    echo "SKIP: checked existing BAM $sample"
    continue
  fi
  temporary="$(mktemp "$ALIGN/.$sample.XXXXXX.bam")"
  hisat2 -x "$HISAT2_INDEX" -1 "$r1" -2 "$r2" -p "$THREADS" --summary-file "$ALIGN/$sample.hisat2.summary.txt" 2> "$log" | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEMORY" -o "$temporary" -
  samtools quickcheck -v "$temporary"
  samtools index "$temporary"
  mv -- "$temporary" "$bam"
  mv -- "$temporary.bai" "$bam.bai"
  temporary=""
  echo "DONE: $sample"
done
