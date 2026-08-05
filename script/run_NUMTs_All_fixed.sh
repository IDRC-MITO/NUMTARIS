#!/bin/bash
#
# run_NUMTs_All_fixed.sh
#
# Batch driver for NUMTs_detection_fixed.sh: runs NUMT candidate detection
# on every *.bam file found directly under <bam_dir>, writing each sample's
# results to <output_root>/<sample>/. Once every sample has been processed,
# it calls run_NUMTs_Fig_html.py to aggregate all samples'
# *.NUMTs_candidates.tsv outputs into circos-style figures (PNG/SVG/PDF)
# and an interactive HTML viewer under <output_root>/output_<label>/.
#
# Usage:
#   ./run_NUMTs_All_fixed.sh <bam_dir> <output_root> [options] [-- extra options]
#
# Options:
#   --skip-figures         Skip the final run_NUMTs_Fig_html.py aggregation step.
#
# Any other arguments (or anything after a literal --) are forwarded
# verbatim to NUMTs_detection_fixed.sh for every sample, e.g. to override
# its detection thresholds:
#   ./run_NUMTs_All_fixed.sh bams/ out/ --min-disc-reads 3
#
set -euo pipefail
shopt -s nullglob

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <bam_dir> <output_root> [--skip-figures] [-- extra options for NUMTs_detection_fixed.sh]" >&2
  exit 1
fi

BAM_BASE="$1"
OUT_ROOT="$2"
shift 2

# Pull --skip-figures out of the remaining args (it can appear anywhere
# before or after a literal --); everything else is forwarded to
# NUMTs_detection_fixed.sh unchanged.
SKIP_FIGURES=0
ARGS=()
for a in "$@"; do
  if [[ "$a" == "--skip-figures" ]]; then
    SKIP_FIGURES=1
  else
    ARGS+=("$a")
  fi
done
set -- "${ARGS[@]+"${ARGS[@]}"}"

[[ "${1:-}" == "--" ]] && shift
EXTRA_OPTS=("$@")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NUMT_SCRIPT="${SCRIPT_DIR}/NUMTs_detection_fixed.sh"
FIG_SCRIPT="${SCRIPT_DIR}/run_NUMTs_Fig_html.py"

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
    "${EXTRA_OPTS[@]+"${EXTRA_OPTS[@]}"}" \
    "${BAM}" \
    "${OUTDIR}"

  echo "Finished: ${SAMPLE}"
done

echo "All samples completed."

if [[ "${SKIP_FIGURES}" -eq 1 ]]; then
  echo "Skipping figure generation (--skip-figures)."
  exit 0
fi

echo "======================================"
echo "Generating aggregated NUMTs figures"
echo "======================================"

if [[ ! -f "${FIG_SCRIPT}" ]]; then
  echo "WARNING: ${FIG_SCRIPT} not found; skipping figure generation." >&2
  exit 0
fi

if ! command -v python3 >/dev/null 2>&1; then
  echo "WARNING: python3 not found on PATH; skipping figure generation." >&2
  exit 0
fi

# Figure generation aggregates across all samples' *.NUMTs_candidates.tsv
# files under OUT_ROOT. Its failure (e.g. missing pandas/matplotlib) should
# not be treated as a failure of the detection pipeline itself, since every
# sample's TSV output is already safely on disk at this point.
if ! python3 "${FIG_SCRIPT}" --input-dir "${OUT_ROOT}" --output-dir "${OUT_ROOT}"; then
  echo "WARNING: figure generation failed. Per-sample NUMTs_candidates.tsv outputs are still available under ${OUT_ROOT}." >&2
  exit 0
fi

echo "All done. Figures written under ${OUT_ROOT}/output_<label>/."
