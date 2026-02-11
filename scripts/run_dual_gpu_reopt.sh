#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

SDF="moe_conformers/cnph_th_cf3_3d_conformers.sdf"
SMILES_CSV="binaph_dimer_smiles.csv"
OUTPUT_DIR="runs"
OUTPUT_PREFIX="moe_mace_reopt"
MODEL="small"
MACE_DTYPE="float64"
MAX_STEPS=500
FMAX=0.005
SAVE_EVERY=25
GPU_IDS="0,1"
PYTHON_BIN="python"
DRY_RUN=0
EXTRA_ARGS=()

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_dual_gpu_reopt.sh [options] [-- <extra reopt args>]

Options:
  --sdf PATH              Input SDF (default: moe_conformers/cnph_th_cf3_3d_conformers.sdf)
  --smiles-csv PATH       Molecule map CSV (default: binaph_dimer_smiles.csv)
  --output-dir DIR        Output directory (default: runs)
  --output-prefix NAME    Output prefix (default: moe_mace_reopt)
  --model NAME            small|medium|large (default: small)
  --mace-dtype TYPE       float64|float32 (default: float64)
  --max-steps N           Max optimizer steps (default: 500)
  --fmax VALUE            Force threshold (default: 0.005)
  --save-every N          Checkpoint interval (default: 25)
  --gpu-ids A,B           Two physical GPU ids (default: 0,1)
  --python BIN            Python executable (default: python)
  --dry-run               Print commands without launching
  -h, --help              Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sdf) SDF="$2"; shift 2 ;;
    --smiles-csv) SMILES_CSV="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --output-prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
    --model) MODEL="$2"; shift 2 ;;
    --mace-dtype) MACE_DTYPE="$2"; shift 2 ;;
    --max-steps) MAX_STEPS="$2"; shift 2 ;;
    --fmax) FMAX="$2"; shift 2 ;;
    --save-every) SAVE_EVERY="$2"; shift 2 ;;
    --gpu-ids) GPU_IDS="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    --)
      shift
      EXTRA_ARGS=("$@")
      break
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ ! -f "${SDF}" ]]; then
  echo "Missing SDF: ${SDF}" >&2
  exit 1
fi

if [[ ! -f "${SMILES_CSV}" ]]; then
  echo "Missing CSV: ${SMILES_CSV}" >&2
  exit 1
fi

IFS=',' read -r GPU0 GPU1 <<< "${GPU_IDS}"
if [[ -z "${GPU0:-}" || -z "${GPU1:-}" ]]; then
  echo "Expected exactly two GPU IDs, got: ${GPU_IDS}" >&2
  exit 1
fi

GPU_COUNT="$("${PYTHON_BIN}" - <<'PY'
import torch
print(torch.cuda.device_count() if torch.cuda.is_available() else 0)
PY
)"
if [[ "${GPU_COUNT}" -lt 2 ]]; then
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "WARNING: dry-run with <2 CUDA devices (found ${GPU_COUNT})"
  else
    echo "Need >=2 visible CUDA devices, found ${GPU_COUNT}" >&2
    exit 1
  fi
fi

TOTAL_CONFS="$("${PYTHON_BIN}" - "${SDF}" <<'PY'
import sys
from rdkit import Chem

sdf = sys.argv[1]
sup = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)
if sup is None:
    raise SystemExit(f"Unable to read SDF: {sdf}")
count = sum(1 for mol in sup if mol is not None)
print(count)
PY
)"

if [[ "${TOTAL_CONFS}" -lt 2 ]]; then
  echo "Not enough conformers in ${SDF}: ${TOTAL_CONFS}" >&2
  exit 1
fi

HALF=$(( (TOTAL_CONFS + 1) / 2 ))
START0=0
END0="${HALF}"
START1="${HALF}"
END1="${TOTAL_CONFS}"

mkdir -p "${OUTPUT_DIR}"

COMMON_ARGS=(
  reoptimize_moe_conformers.py
  --sdf "${SDF}"
  --smiles-csv "${SMILES_CSV}"
  --model "${MODEL}"
  --device cuda
  --mace-dtype "${MACE_DTYPE}"
  --max-steps "${MAX_STEPS}"
  --fmax "${FMAX}"
  --save-every "${SAVE_EVERY}"
  --resume
)

if [[ "${#EXTRA_ARGS[@]}" -gt 0 ]]; then
  COMMON_ARGS+=("${EXTRA_ARGS[@]}")
fi

CMD0=(
  "${PYTHON_BIN}" "${COMMON_ARGS[@]}"
  --start-index "${START0}" --end-index "${END0}"
  --output-prefix "${OUTPUT_DIR}/${OUTPUT_PREFIX}_gpu0"
)
CMD1=(
  "${PYTHON_BIN}" "${COMMON_ARGS[@]}"
  --start-index "${START1}" --end-index "${END1}"
  --output-prefix "${OUTPUT_DIR}/${OUTPUT_PREFIX}_gpu1"
)

LOG0="${OUTPUT_DIR}/${OUTPUT_PREFIX}_gpu0.log"
LOG1="${OUTPUT_DIR}/${OUTPUT_PREFIX}_gpu1.log"

echo "Total conformers: ${TOTAL_CONFS}"
echo "GPU0=${GPU0}: [${START0}, ${END0}) -> ${OUTPUT_PREFIX}_gpu0"
echo "GPU1=${GPU1}: [${START1}, ${END1}) -> ${OUTPUT_PREFIX}_gpu1"
echo "Output dir: ${OUTPUT_DIR}"
echo

echo "Command GPU0:"
printf 'CUDA_VISIBLE_DEVICES=%s ' "${GPU0}"
printf '%q ' "${CMD0[@]}"
echo
echo "Command GPU1:"
printf 'CUDA_VISIBLE_DEVICES=%s ' "${GPU1}"
printf '%q ' "${CMD1[@]}"
echo
echo

if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo "Dry-run mode: not launching jobs."
  exit 0
fi

CUDA_VISIBLE_DEVICES="${GPU0}" "${CMD0[@]}" >"${LOG0}" 2>&1 &
PID0=$!
CUDA_VISIBLE_DEVICES="${GPU1}" "${CMD1[@]}" >"${LOG1}" 2>&1 &
PID1=$!

echo "Launched PID0=${PID0}, PID1=${PID1}"
echo "Logs:"
echo "  ${LOG0}"
echo "  ${LOG1}"

wait "${PID0}"
STATUS0=$?
wait "${PID1}"
STATUS1=$?

if [[ "${STATUS0}" -ne 0 || "${STATUS1}" -ne 0 ]]; then
  echo "At least one GPU job failed (status0=${STATUS0}, status1=${STATUS1})." >&2
  exit 1
fi

echo
echo "Both GPU chunks finished. Merging outputs..."
"${PYTHON_BIN}" scripts/merge_reopt_chunks.py \
  --input-glob "${OUTPUT_DIR}/${OUTPUT_PREFIX}_gpu*_all_conformers.csv" \
  --output-dir "${OUTPUT_DIR}" \
  --output-prefix "${OUTPUT_PREFIX}"

echo "Dual-GPU run complete."
