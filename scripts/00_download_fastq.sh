#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
META="${SAMPLE_SHEET:-$PROJECT_ROOT/metadata/sample_sheet.tsv}"
RAW="${RAW_FASTQ_DIR:-$PROJECT_ROOT/raw_fastq}"

[[ -s "$META" ]] || { echo "ERROR: sample sheet not found: $META" >&2; exit 2; }
command -v curl >/dev/null || { echo "ERROR: curl is required." >&2; exit 2; }
command -v gzip >/dev/null || { echo "ERROR: gzip is required." >&2; exit 2; }
mkdir -p "$RAW"

download_one() {
  local url="$1" destination="$2" temporary
  temporary="${destination}.part"
  if [[ -s "$destination" ]]; then
    if gzip -t "$destination"; then
      echo "SKIP: valid file exists: $destination"
      return 0
    fi
    echo "ERROR: existing file is not a valid gzip; preserve and inspect it: $destination" >&2
    return 1
  fi
  [[ "$url" == http://* || "$url" == https://* ]] || url="https://$url"
  echo "GET: $url"
  if ! curl --fail --location --retry 3 --retry-delay 2 --connect-timeout 30 --output "$temporary" "$url"; then
    rm -f -- "$temporary"
    return 1
  fi
  if ! gzip -t "$temporary"; then
    rm -f -- "$temporary"
    echo "ERROR: downloaded file failed gzip integrity check: $url" >&2
    return 1
  fi
  mv -- "$temporary" "$destination"
  echo "SAVED: $destination"
}

header="$(head -n 1 "$META")"
header="${header%$'\r'}"
[[ "$header" == $'sample_id\tcondition\trun_accession\tfastq_1\tfastq_2' ]] || {
  echo "ERROR: unexpected sample-sheet header in $META" >&2
  exit 2
}

sample_count=0
seen_sample_ids=$'\n'
while IFS=$'\t' read -r sample_id condition run_accession fastq_1 fastq_2; do
  sample_id="${sample_id%$'\r'}"
  condition="${condition%$'\r'}"
  run_accession="${run_accession%$'\r'}"
  fastq_1="${fastq_1%$'\r'}"
  fastq_2="${fastq_2%$'\r'}"
  [[ -n "${sample_id:-}" ]] || continue
  [[ "$sample_id" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "ERROR: unsafe sample_id: $sample_id" >&2; exit 2; }
  [[ "$seen_sample_ids" != *$'\n'"$sample_id"$'\n'* ]] || { echo "ERROR: duplicate sample_id: $sample_id" >&2; exit 2; }
  seen_sample_ids+="$sample_id"$'\n'
  [[ -n "${condition:-}" && -n "${run_accession:-}" ]] || { echo "ERROR: condition or run accession missing for $sample_id" >&2; exit 2; }
  [[ -n "${fastq_1:-}" && -n "${fastq_2:-}" ]] || { echo "ERROR: FASTQ paths are missing for $sample_id" >&2; exit 2; }
  download_one "$fastq_1" "$RAW/${sample_id}_1.fastq.gz"
  download_one "$fastq_2" "$RAW/${sample_id}_2.fastq.gz"
  sample_count=$((sample_count + 1))
done < <(tail -n +2 "$META")

(( sample_count > 0 )) || { echo "ERROR: sample sheet contains no samples." >&2; exit 2; }
echo "FASTQ inputs ready for $sample_count samples."
