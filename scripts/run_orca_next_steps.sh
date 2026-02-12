#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

ORCA_BIN="${ORCA_BIN:-orca}"
PYTHON_BIN="${PYTHON_BIN:-python}"
WORKERS="${WORKERS:-2}"
NPROCS_PER_JOB="${NPROCS_PER_JOB:-12}"
MAXCORE="${MAXCORE:-3000}"
STAGES="${STAGES:-vertical,s1opt}"
DRY_RUN="${DRY_RUN:-0}"
MERGED_CSV="${MERGED_CSV:-runs/moe_mace_reopt_all_conformers.csv}"
SDF_PATH="${SDF_PATH:-moe_conformers/cnph_th_cf3_3d_conformers.sdf}"

if [[ ! -f "${MERGED_CSV}" ]]; then
  echo "Missing ${MERGED_CSV}" >&2
  echo "Finish and merge MACE run before starting ORCA phase." >&2
  exit 1
fi

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
  if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
  else
    echo "Neither ${PYTHON_BIN} nor python3 is available in PATH." >&2
    exit 1
  fi
fi

echo "== ORCA Phase: candidate selection =="
"${PYTHON_BIN}" scripts/orca_select_candidates.py \
  --input-csv "${MERGED_CSV}" \
  --output-csv runs/orca_candidates.csv \
  --top-per-molecule 2 \
  --max-total 96 \
  --max-theta 50 \
  --max-distance 4.6 \
  --min-overlap 30

echo
echo "== ORCA Phase: input preparation =="
"${PYTHON_BIN}" scripts/orca_prepare_inputs.py \
  --candidates-csv runs/orca_candidates.csv \
  --sdf "${SDF_PATH}" \
  --output-dir runs/orca_jobs \
  --stages "${STAGES}" \
  --nprocs "${NPROCS_PER_JOB}" \
  --maxcore "${MAXCORE}" \
  --functional CAM-B3LYP \
  --basis def2-SVP \
  --dispersion D3BJ \
  --solvent toluene \
  --nroots 10 \
  --iroot 1 \
  --tda

echo
echo "== ORCA Phase: queue execution =="
QUEUE_ARGS=(
  --manifest runs/orca_jobs/orca_job_manifest.csv
  --orca-bin "${ORCA_BIN}"
  --max-workers "${WORKERS}"
  --resume
  --summary-csv runs/orca_jobs/orca_run_summary.csv
)
if [[ "${DRY_RUN}" == "1" ]]; then
  QUEUE_ARGS+=(--dry-run)
fi
"${PYTHON_BIN}" scripts/orca_run_queue.py "${QUEUE_ARGS[@]}"

if [[ "${DRY_RUN}" == "1" ]]; then
  echo
  echo "Dry-run mode: skipping ORCA output collection step."
  echo "Review:"
  echo "  runs/orca_jobs/orca_job_manifest.csv"
  echo "  runs/orca_jobs/orca_run_summary.csv"
  exit 0
fi

echo
echo "== ORCA Phase: collect results =="
"${PYTHON_BIN}" scripts/orca_collect_results.py \
  --jobs-root runs/orca_jobs \
  --glob "**/*.out" \
  --output-csv runs/orca_jobs/orca_results_summary.csv

echo
echo "Done. Review:"
echo "  runs/orca_jobs/orca_run_summary.csv"
echo "  runs/orca_jobs/orca_results_summary.csv"
