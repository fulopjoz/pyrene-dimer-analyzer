#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

echo "== Cluster Preflight =="
echo "repo: ${ROOT_DIR}"
echo

echo "-- Python and package check --"
python - <<'PY'
import importlib
import platform
import sys

print("python_version:", sys.version.replace("\n", " "))
print("platform:", platform.platform())

mods = ["torch", "rdkit", "ase", "mace", "xtb", "pandas", "numpy"]
for name in mods:
    try:
        mod = importlib.import_module(name)
        ver = getattr(mod, "__version__", "unknown")
        print(f"{name}: OK ({ver})")
    except Exception as exc:
        print(f"{name}: MISSING ({exc})")

try:
    import torch
    print("torch.cuda.is_available:", torch.cuda.is_available())
    print("torch.cuda.device_count:", torch.cuda.device_count())
    if torch.cuda.is_available():
        for idx in range(torch.cuda.device_count()):
            props = torch.cuda.get_device_properties(idx)
            vram_gb = props.total_memory / (1024**3)
            print(f"gpu[{idx}]: {props.name}, vram_gb={vram_gb:.2f}")
except Exception as exc:
    print("torch cuda probe failed:", exc)
PY
echo

echo "-- nvidia-smi --"
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi
else
  echo "nvidia-smi not found"
fi
echo

echo "-- Required input files --"
for p in \
  "binaph_dimer_smiles.csv" \
  "moe_conformers/cnph_th_cf3_3d_conformers.sdf" \
  "validate_mace.py" \
  "reoptimize_moe_conformers.py"
do
  if [[ -f "${p}" ]]; then
    echo "OK: ${p}"
  else
    echo "MISSING: ${p}"
  fi
done
echo

echo "-- Quick sanity commands --"
echo "python -m pytest tests/test_validate_mace.py -q"
echo "python validate_mace.py --all-methods --output validation_results.csv"
echo
echo "Preflight complete."
