#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
META="${SAMPLE_SHEET:-$PROJECT_ROOT/metadata/sample_sheet.tsv}"
MANIFEST="${FASTQ_MANIFEST:-$PROJECT_ROOT/metadata/fastq_manifest.tsv}"
RAW="${RAW_FASTQ_DIR:-$PROJECT_ROOT/raw_fastq}"

[[ -s "$META" ]] || { echo "ERROR: sample sheet not found: $META" >&2; exit 2; }
[[ -s "$MANIFEST" ]] || { echo "ERROR: FASTQ checksum manifest not found: $MANIFEST" >&2; exit 2; }
for tool in curl gzip md5sum stat cut; do
  command -v "$tool" >/dev/null || { echo "ERROR: $tool is required." >&2; exit 2; }
done
mkdir -p "$RAW"

manifest_header="$(head -n 1 "$MANIFEST")"
manifest_header="${manifest_header#$'\xef\xbb\xbf'}"
manifest_header="${manifest_header%$'\r'}"
[[ "$manifest_header" == $'sample_id\tmate\texpected_bytes\tmd5' ]] || {
  echo "ERROR: unexpected FASTQ manifest header in $MANIFEST" >&2
  exit 2
}
declare -A expected_md5 expected_bytes
while IFS=$'\t' read -r sample_id mate expected_size expected_digest; do
  sample_id="${sample_id%$'\r'}"
  mate="${mate%$'\r'}"
  expected_size="${expected_size%$'\r'}"
  expected_digest="${expected_digest%$'\r'}"
  [[ -n "${sample_id:-}" ]] || continue
  [[ "$sample_id" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "ERROR: unsafe sample_id in manifest: $sample_id" >&2; exit 2; }
  [[ "$mate" == 1 || "$mate" == 2 ]] || { echo "ERROR: invalid mate in manifest for $sample_id: $mate" >&2; exit 2; }
  [[ "$expected_size" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid expected byte count for $sample_id mate $mate" >&2; exit 2; }
  [[ "$expected_digest" =~ ^[[:xdigit:]]{32}$ ]] || { echo "ERROR: invalid MD5 for $sample_id mate $mate" >&2; exit 2; }
  key="$sample_id:$mate"
  [[ -z "${expected_md5[$key]+x}" ]] || { echo "ERROR: duplicate checksum entry for $sample_id mate $mate" >&2; exit 2; }
  expected_md5[$key]="${expected_digest,,}"
  expected_bytes[$key]="$expected_size"
done < <(tail -n +2 "$MANIFEST")

download_one() {
  local sample_id="$1" mate="$2" url="$3" destination="$4" temporary key expected_digest expected_size actual_size actual_digest
  key="$sample_id:$mate"
  expected_digest="${expected_md5[$key]:-}"
  expected_size="${expected_bytes[$key]:-}"
  [[ -n "$expected_digest" && -n "$expected_size" ]] || {
    echo "ERROR: checksum manifest has no entry for $sample_id mate $mate" >&2
    return 1
  }
  temporary="${destination}.part"
  if [[ -s "$destination" ]]; then
    actual_size="$(stat -L -c '%s' -- "$destination")"
    actual_digest="$(md5sum -- "$destination" | cut -d ' ' -f 1)"
    if [[ "$actual_size" == "$expected_size" && "${actual_digest,,}" == "$expected_digest" ]] && gzip -t "$destination"; then
      echo "SKIP: ENA size, MD5, and gzip checks passed: $destination"
      return 0
    fi
    echo "ERROR: existing file failed ENA size, MD5, or gzip verification; preserve it and use a clean RAW_FASTQ_DIR: $destination" >&2
    return 1
  fi
  [[ "$url" == http://* || "$url" == https://* ]] || url="https://$url"
  echo "GET: $url"
  if ! curl --fail --location --retry 3 --retry-delay 2 --connect-timeout 30 --output "$temporary" "$url"; then
    rm -f -- "$temporary"
    return 1
  fi
  actual_size="$(stat -c '%s' -- "$temporary")"
  actual_digest="$(md5sum -- "$temporary" | cut -d ' ' -f 1)"
  if [[ "$actual_size" != "$expected_size" || "${actual_digest,,}" != "$expected_digest" ]]; then
    rm -f -- "$temporary"
    echo "ERROR: downloaded FASTQ failed ENA byte-count or MD5 verification: $url" >&2
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
header="${header#$'\xef\xbb\xbf'}"
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
  download_one "$sample_id" 1 "$fastq_1" "$RAW/${sample_id}_1.fastq.gz"
  download_one "$sample_id" 2 "$fastq_2" "$RAW/${sample_id}_2.fastq.gz"
  sample_count=$((sample_count + 1))
done < <(tail -n +2 "$META")

(( sample_count > 0 )) || { echo "ERROR: sample sheet contains no samples." >&2; exit 2; }
echo "FASTQ inputs ready for $sample_count samples."
