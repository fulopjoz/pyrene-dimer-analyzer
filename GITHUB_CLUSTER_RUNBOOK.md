# GitHub + Cluster Runbook

## Scope
This runbook describes how to:
1. Stop local long-running jobs.
2. Prepare a clean, portable GitHub repository.
3. Clone and run on a GPU cluster with faster throughput.
4. Run large re-optimization workloads safely (resume + chunking).

The repo was updated for performance and portability:
- `reoptimize_moe_conformers.py` now supports:
  - `--device auto|cpu|cuda` (default `auto`)
  - default model `small`
  - `--fmax`
  - `--mace-dtype float64|float32`
  - `--save-every`
  - `--start-index` / `--end-index` (chunked jobs)
- `pyrene_analyzer/mace_optimizer.py` now supports `default_dtype`.
- `.gitignore` updated to avoid committing generated outputs/local binaries.

## 1. Stop Local Jobs (if running)
From Windows PowerShell:

```powershell
wsl -d Ubuntu-24.04 -- bash -lc "ps -ef | rg -n 'reoptimize_moe_conformers.py|validate_mace.py|run_screening_xtb.py' || true"
```

If jobs are found and should be terminated:

```powershell
wsl -d Ubuntu-24.04 -- bash -lc "pkill -f 'reoptimize_moe_conformers.py|validate_mace.py|run_screening_xtb.py' || true"
```

Verify no residual jobs:

```powershell
wsl -d Ubuntu-24.04 -- bash -lc "ps -ef | rg -n 'reoptimize_moe_conformers.py|validate_mace.py|run_screening_xtb.py' || true"
```

## 2. Validate Local Repo Before Push

### 2.1 Confirm no file exceeds 49 MB

```powershell
@'
import os
from pathlib import Path
limit = 49 * 1024 * 1024
big = []
for dp, _, fns in os.walk('.'):
    if '.git' in dp.split(os.sep):
        continue
    for fn in fns:
        p = Path(dp) / fn
        try:
            sz = p.stat().st_size
        except OSError:
            continue
        if sz > limit:
            big.append((sz, str(p)))
print('count>', len(big))
for sz, p in sorted(big, reverse=True):
    print(round(sz/1024/1024,2), 'MB', p)
'@ | python -
```

### 2.2 Run quick correctness checks

```powershell
python -m pytest tests/test_validate_mace.py
wsl -d Ubuntu-24.04 -- bash -lc "cd /mnt/c/Users/dodo/Documents/projects/pyrene-dimer-analyzer && ./bin/micromamba run -n pyrene-xtb python validate_mace.py --all-methods --output validation_results.csv"
```

## 3. Prepare Clean Commit

### 3.1 Review what will be committed

```powershell
git status
git diff -- . ':!*.csv' ':!*.png'
```

### 3.2 Stage recommended files for compute portability
This stages code, configs, docs, and required input data only.

```powershell
git add `
  .gitignore `
  pyproject.toml `
  requirements.txt `
  requirements-dev.txt `
  pyrene_analyzer `
  validate_mace.py `
  reoptimize_moe_conformers.py `
  run_screening_xtb.py `
  compare_methods.py `
  tests/test_validate_mace.py `
  binaph_dimer_smiles.csv `
  moe_conformers/cnph_th_cf3_3d_conformers.sdf `
  GITHUB_CLUSTER_RUNBOOK.md `
  README.md
```

Then inspect staged files:

```powershell
git diff --cached --name-status
```

Commit:

```powershell
git commit -m "Optimize MACE runtime workflow, add chunked GPU execution controls, and cluster runbook"
```

## 4. Create and Push GitHub Repository

If repo is not initialized remotely yet:

```powershell
git remote add origin https://github.com/<your-org-or-user>/<repo-name>.git
git branch -M main
git push -u origin main
```

If remote already exists:

```powershell
git push
```

## 5. Clone and Set Up on Cluster

## 5.1 Clone

```bash
git clone https://github.com/<your-org-or-user>/<repo-name>.git
cd <repo-name>
```

## 5.2 Create environment
Using micromamba (recommended):

```bash
micromamba create -n pyrene-mace python=3.11 -y
micromamba activate pyrene-mace
python -m pip install --upgrade pip
python -m pip install -e ".[mace]"
```

If your cluster requires explicit CUDA torch wheel, install torch first per cluster policy, then:

```bash
python -m pip install -e .
python -m pip install mace-torch ase
```

## 5.3 Verify GPU is visible

```bash
nvidia-smi
python - <<'PY'
import torch
print("torch", torch.__version__)
print("cuda_available", torch.cuda.is_available())
if torch.cuda.is_available():
    print("device", torch.cuda.get_device_name(0))
PY
```

## 6. Smoke Tests on Cluster

```bash
python -m pytest tests/test_validate_mace.py
python validate_mace.py --all-methods --output validation_results.csv
python reoptimize_moe_conformers.py --test --model small --device auto --max-steps 500 --fmax 0.005 --output-prefix smoke
```

## 7. Production Runs

### 7.1 Single-job (simple)

```bash
python reoptimize_moe_conformers.py \
  --model small \
  --device auto \
  --mace-dtype float64 \
  --max-steps 500 \
  --fmax 0.005 \
  --save-every 25 \
  --resume \
  --output-prefix moe_mace_reopt
```

### 7.2 Faster (if acceptable): float32 on GPU

```bash
python reoptimize_moe_conformers.py \
  --model small \
  --device cuda \
  --mace-dtype float32 \
  --max-steps 500 \
  --fmax 0.005 \
  --save-every 25 \
  --resume \
  --output-prefix moe_mace_reopt_f32
```

### 7.3 Chunked array jobs (recommended for clusters)
Use index windows to split the 3347 conformers.

Example shell logic per chunk:

```bash
CHUNK=200
START=$((TASK_ID * CHUNK))
END=$((START + CHUNK))

python reoptimize_moe_conformers.py \
  --model small \
  --device cuda \
  --mace-dtype float64 \
  --max-steps 500 \
  --fmax 0.005 \
  --save-every 25 \
  --start-index ${START} \
  --end-index ${END} \
  --output-prefix runs/moe_mace_chunk_${TASK_ID}
```

## 8. Merge Chunk Outputs

```bash
python - <<'PY'
from pathlib import Path
import pandas as pd

parts = sorted(Path("runs").glob("moe_mace_chunk_*_all_conformers.csv"))
dfs = [pd.read_csv(p) for p in parts if p.exists()]
if not dfs:
    raise SystemExit("No chunk files found")
df = pd.concat(dfs, ignore_index=True).drop_duplicates(subset=["name", "conformer_id"])
df.to_csv("moe_mace_reopt_all_conformers.csv", index=False)
print("merged rows:", len(df))
PY
```

Then regenerate summaries if needed by rerunning aggregation/analysis scripts.

## 9. Recommended Cluster Practices
1. Use local scratch for intermediate outputs, then copy final CSVs to project storage.
2. Keep `--save-every` modest (e.g., 25-100) to reduce I/O overhead.
3. Prefer `small + cuda` for bulk, then re-run top candidates with `medium`.
4. Keep `--resume` enabled for long jobs.

## 10. Notes
1. `cuequivariance` acceleration is optional; absence does not block runs.
2. `float64` is scientifically safer; `float32` is faster.
3. Generated output files are intentionally ignored by `.gitignore` to keep GitHub clean.
