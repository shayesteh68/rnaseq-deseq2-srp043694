# Copy to scripts/config.local.sh, edit paths, then source scripts/config.local.sh.
export PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export REF_DIR="/path/to/reference"
export REF_GTF="$REF_DIR/Homo_sapiens.GRCh38.110.gtf"
export HISAT2_INDEX="$REF_DIR/Homo_sapiens_index"
export THREADS=2
export SORT_THREADS=1
export SORT_MEMORY=512M
export STRANDNESS=0
