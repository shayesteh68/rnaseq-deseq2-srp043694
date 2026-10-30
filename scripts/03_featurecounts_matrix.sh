#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
ALIGN="${ALIGN_BAM_DIR:-$PROJECT_ROOT/data/aligned_bam}"
COUNTS_DIR="${COUNTS_DIR:-$PROJECT_ROOT/results/counts}"
REF_GTF="${REF_GTF:-}"
THREADS="${THREADS:-2}"
STRANDNESS="${STRANDNESS:-0}"
[[ "$STRANDNESS" =~ ^[012]$ ]] || { echo "ERROR: STRANDNESS must be 0, 1, or 2." >&2; exit 2; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: THREADS must be a positive integer." >&2; exit 2; }
[[ -n "$REF_GTF" && -s "$REF_GTF" ]] || { echo "ERROR: set REF_GTF to the matching annotation GTF." >&2; exit 2; }
command -v featureCounts >/dev/null || { echo "ERROR: featureCounts is required." >&2; exit 2; }
mkdir -p "$COUNTS_DIR"
shopt -s nullglob
bams=("$ALIGN"/*.sorted.bam)
(( ${#bams[@]} > 0 )) || { echo "ERROR: no sorted BAMs found in $ALIGN" >&2; exit 1; }
sample_csv=""
for bam in "${bams[@]}"; do
  sample="$(basename "$bam" .sorted.bam)"
  [[ "$sample" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "ERROR: unsupported sample name: $sample" >&2; exit 2; }
  sample_csv="${sample_csv:+$sample_csv,}$sample"
done
work="$(mktemp -d "$COUNTS_DIR/.featurecounts.XXXXXX")"
cleanup() { rm -f -- "$work/counts.txt" "$work/counts.txt.summary" "$work/matrix.tsv"; rmdir -- "$work" 2>/dev/null || true; }
trap cleanup EXIT
featureCounts -T "$THREADS" -a "$REF_GTF" -o "$work/counts.txt" -t exon -g gene_id -s "$STRANDNESS" -p --countReadPairs -B -C "${bams[@]}"
fc_version="$(sed -n '1s/^# Program:featureCounts v\([^;]*\);.*/\1/p' "$work/counts.txt")"
[[ -n "$fc_version" ]] || { echo "ERROR: could not identify featureCounts version in its output." >&2; exit 1; }
# Keep useful provenance while removing machine-specific reference, output, and BAM paths.
{
  printf '# Program:featureCounts v%s; Command: -T %s -t exon -g gene_id -s %s -p --countReadPairs -B -C; Annotation: %s; Samples: %s\n' \
    "$fc_version" "$THREADS" "$STRANDNESS" "$(basename "$REF_GTF")" "$sample_csv"
  tail -n +2 "$work/counts.txt"
} > "$work/counts.portable.txt"
mv -- "$work/counts.portable.txt" "$work/counts.txt"
awk -F '\t' -v OFS='\t' -v sample_csv="$sample_csv" '
BEGIN {
  n=split(sample_csv,samples,","); printf "gene_id"
  for(i=1;i<=n;i++) printf OFS "%s",samples[i]
  printf "\n"
}
/^#/ || $1=="Geneid" || $1 ~ /^__/ {next}
{
  if(NF!=n+6) {bad=1; exit 2}
  for(i=7;i<=NF;i++) if($i !~ /^[0-9]+$/) {bad=1; exit 2}
  printf "%s",$1
  for(i=7;i<=NF;i++) printf OFS "%s",$i
  printf "\n"; rows++
}
END {if(bad || rows==0) exit 2}
' "$work/counts.txt" > "$work/matrix.tsv"
raw_counts="$COUNTS_DIR/SRP043694_featureCounts_s$STRANDNESS.txt"
matrix="$COUNTS_DIR/SRP043694_count_matrix_s$STRANDNESS.tsv"
mv -- "$work/counts.txt" "$raw_counts"
mv -- "$work/counts.txt.summary" "$raw_counts.summary"
mv -- "$work/matrix.tsv" "$matrix"
echo "Saved matrix: $matrix"
echo "Review sample IDs and biological replication against metadata/sample_sheet.tsv."
