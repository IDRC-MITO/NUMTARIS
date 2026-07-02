#!/bin/bash
#
# run_NUMTs_All_fixed.sh
#
# Batch driver for NUMTs_detection_fixed.sh: runs NUMT candidate detection
# on every *.bam file found directly under <bam_dir>, writing each sample's
# results to <output_root>/<sample>/.
#
# Usage:
#   ./run_NUMTs_All_fixed.sh <bam_dir> <output_root> [-- extra options]
#
# Any arguments after the two required positional ones (or after a literal
# --) are forwarded verbatim to NUMTs_detection_fixed.sh for every sample,
# e.g. to override its detection thresholds:
#   ./run_NUMTs_All_fixed.sh bams/ out/ --min-disc-reads 3
#
set -euo pipefail
shopt -s nullglob

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <bam_dir> <output_root> [-- extra options for NUMTs_detection_fixed.sh]" >&2
  exit 1
fi

BAM_BASE="$1"
OUT_ROOT="$2"
shift 2
[[ "${1:-}" == "--" ]] && shift
EXTRA_OPTS=("$@")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NUMT_SCRIPT="${SCRIPT_DIR}/NUMTs_detection_fixed.sh"

mkdir -p "${OUT_ROOT}"

[[ -d "${BAM_BASE}" ]] || { echo "ERROR: BAM_BASE not found: ${BAM_BASE}" >&2; exit 1; }
[[ -f "${NUMT_SCRIPT}" ]] || { echo "ERROR: NUMT_SCRIPT not found: ${NUMT_SCRIPT}" >&2; exit 1; }

BAMS=( "${BAM_BASE}"/*.bam )

if (( ${#BAMS[@]} == 0 )); then
  echo "No BAM files found in ${BAM_BASE}"
  exit 0
fi

for BAM in "${BAMS[@]}"; do
  SAMPLE=$(basename "${BAM}" .bam)
  OUTDIR="${OUT_ROOT}/${SAMPLE}"

  echo "======================================"
  echo "Processing sample: ${SAMPLE}"
  echo "BAM: ${BAM}"
  echo "OUTDIR: ${OUTDIR}"
  echo "======================================"

  mkdir -p "${OUTDIR}"

  bash "${NUMT_SCRIPT}" \
    "${EXTRA_OPTS[@]}" \
    "${BAM}" \
    "${OUTDIR}"

  echo "Finished: ${SAMPLE}"
done

echo "All samples completed."
