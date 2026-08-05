#!/bin/bash
#
# NUMTs_detection_fixed.sh
#
# Detects candidate NUMTs (nuclear-embedded mitochondrial DNA segments) in a
# single whole-genome sequencing BAM by extracting MT-nuclear discordant and
# split-read pairs and clustering them into candidate insertion sites.
#
# Pipeline:
#   1. Index the BAM if no .bai is present.
#   2. Extract read pairs where either mate maps to the MT contig
#      (chrM/MT/M), then split them into discordant vs. split-read SAM
#      files with samblaster.
#   3. Cluster the nuclear-genome anchors of those reads and score each
#      cluster for read support via numtAnchorCluster.py.
#
# Requirements: samtools, samblaster, python3 (pandas) on PATH.
#
# Usage:
#   ./NUMTs_detection_fixed.sh [options] <input.bam> <output_dir>
#
# Options (all thresholds forwarded to numtAnchorCluster.py):
#   --min-mapq N          Minimum mapping quality for anchor/split reads (default: 20)
#   --min-disc-reads N    Minimum discordant read pairs to call a candidate (default: 2)
#   --min-split-reads N   Minimum split reads to mark PASS_DISC_AND_SPLIT (default: 1)
#   --max-cluster-gap N   Max bp gap between anchors in the same cluster (default: 500)
#   -h, --help            Show this help message
#
set -euo pipefail

MIN_MAPQ=20
MIN_DISC_READS=2
MIN_SPLIT_READS=1
MAX_CLUSTER_GAP=500

usage() {
  echo "Usage: $0 [options] <input.bam> <output_dir>" >&2
  echo "  --min-mapq N          (default: ${MIN_MAPQ})" >&2
  echo "  --min-disc-reads N    (default: ${MIN_DISC_READS})" >&2
  echo "  --min-split-reads N   (default: ${MIN_SPLIT_READS})" >&2
  echo "  --max-cluster-gap N   (default: ${MAX_CLUSTER_GAP})" >&2
  exit 1
}

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
    --min-disc-reads) MIN_DISC_READS="$2"; shift 2 ;;
    --min-split-reads) MIN_SPLIT_READS="$2"; shift 2 ;;
    --max-cluster-gap) MAX_CLUSTER_GAP="$2"; shift 2 ;;
    -h|--help) usage ;;
    --) shift; POSITIONAL+=("$@"); break ;;
    -*) echo "Unknown option: $1" >&2; usage ;;
    *) POSITIONAL+=("$1"); shift ;;
  esac
done

[[ ${#POSITIONAL[@]} -eq 2 ]] || usage

INPUTBAM="${POSITIONAL[0]}"
OUTPUTDIR="${POSITIONAL[1]}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLUSTERSCRIPT="${SCRIPT_DIR}/numtAnchorCluster.py"

SAMPLEID=$(basename "${INPUTBAM}" .bam)

mkdir -p "${OUTPUTDIR}"

INPUTDISC="${OUTPUTDIR}/${SAMPLEID}.mt.disc.sam"
INPUTSPLIT="${OUTPUTDIR}/${SAMPLEID}.mt.split.sam"
OUTTSV="${OUTPUTDIR}/${SAMPLEID}.NUMTs_candidates.tsv"

[[ -f "${INPUTBAM}" ]] || { echo "ERROR: BAM not found: ${INPUTBAM}" >&2; exit 1; }
[[ -f "${CLUSTERSCRIPT}" ]] || { echo "ERROR: CLUSTERSCRIPT not found: ${CLUSTERSCRIPT}" >&2; exit 1; }

echo "=== NUMTs detection: ${SAMPLEID} ==="

echo "[1] Checking BAM index"
if [[ ! -f "${INPUTBAM}.bai" ]]; then
  echo "BAM index not found. Creating index..."
  samtools index "${INPUTBAM}"
fi

echo "[2] Extracting MT-related discordant and split reads"

# -F 2 drops reads flagged as a "properly paired" alignment, since those are
# not candidates for a nuclear-MT discordant/split pair.
# The awk filter keeps only the header and read pairs where either the read
# itself (col 3, RNAME) or its mate (col 7, RNEXT) maps to the MT contig,
# under any of its common names (chrM/MT/M) -- this must stay in sync with
# the MT_NAMES set in numtAnchorCluster.py.
# samblaster then re-derives, from those pairs, which are discordant
# (mates on different chromosomes) vs. split (soft-clipped/supplementary
# alignments) — writing each category to its own SAM file. -o /dev/null
# discards the (unneeded) copy of all reads that samblaster normally passes
# through.
samtools view -h -F 2 "${INPUTBAM}" \
  | awk 'BEGIN{OFS="\t"} /^@/ {print; next} ($3=="chrM" || $3=="MT" || $3=="M" || $7=="chrM" || $7=="MT" || $7=="M") {print}' \
  | samtools sort -n - \
  | samtools view -h - \
  | samblaster \
      --ignoreUnmated \
      -e \
      -d "${INPUTDISC}" \
      -s "${INPUTSPLIT}" \
      -o /dev/null

DISC_LINES=$(grep -vc '^@' "${INPUTDISC}" || true)
SPLIT_LINES=$(grep -vc '^@' "${INPUTSPLIT}" || true)

echo "Discordant SAM records: ${DISC_LINES}"
echo "Split SAM records: ${SPLIT_LINES}"

if [[ "${DISC_LINES}" -eq 0 ]]; then
  echo "No discordant MT-related reads. Creating empty result."
  echo -e "SampleID\tchr\tstart\tend\tMT_start\tMT_end\tNUMT_discordant_reads\tNUMT_split_reads\tNUMT_total_reads\tmean_MAPQ\tfilter_status\tdiscFile\tsplitFile\twgsBAM" > "${OUTTSV}"
  exit 0
fi

echo "[3] Clustering and filtering NUMTs candidates"

python3 "${CLUSTERSCRIPT}" \
  --sample "${SAMPLEID}" \
  --bam "${INPUTBAM}" \
  --disc "${INPUTDISC}" \
  --split "${INPUTSPLIT}" \
  --out "${OUTTSV}" \
  --min-mapq "${MIN_MAPQ}" \
  --min-disc-reads "${MIN_DISC_READS}" \
  --min-split-reads "${MIN_SPLIT_READS}" \
  --max-cluster-gap "${MAX_CLUSTER_GAP}"

echo "[4] Done"
echo "Output: ${OUTTSV}"
