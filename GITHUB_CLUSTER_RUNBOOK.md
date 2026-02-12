# GitHub + Cluster Runbook

## Scope
This runbook is the operational checklist to:
1. Finalize and push this repository to GitHub.
2. Run reproducible validation locally.
3. Move production compute to a stronger Linux server (24 CPU + 2 GPUs, no scheduler yet).
4. Execute full MACE re-optimization efficiently with dual-GPU chunking.

## Current Technical Status (Verified)
1. `validate_mace.py` is wired correctly:
   - xTB path calls `optimize_with_gfn2xtb` (no stale `optimize_conformer` import).
   - Benchmark starts at `initial_distance=4.0 A`.
   - Tight convergence is used (`fmax=0.005`, higher step limits).
2. Validation benchmark output is saved in `validation_results.csv`.
3. Large-file check passed: no files >49 MB in the repo.
4. Local WSL environment detects CUDA (`torch.cuda.is_available() == True`).

## 1. Stop Local Jobs Before Packaging
From PowerShell:

```powershell
wsl -d Ubuntu-24.04 -- bash -lc "ps -ef | rg -n 'reoptimize_moe_conformers.py|validate_mace.py|run_screening_xtb.py' || true"
```

If jobs should be terminated:

```powershell
wsl -d Ubuntu-24.04 -- bash -lc "pkill -f 'reoptimize_moe_conformers.py|validate_mace.py|run_screening_xtb.py' || true"
```

Recheck:

```powershell
wsl -d Ubuntu-24.04 -- bash -lc "ps -ef | rg -n 'reoptimize_moe_conformers.py|validate_mace.py|run_screening_xtb.py' || true"
```

## 2. Pre-Push Validation

### 2.1 Verify File Size Policy (<49 MB)

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

Expected output for a clean repo is exactly `count> 0` (meaning zero files exceeded 49 MB).

### 2.2 Run Regression + Benchmark Checks

```powershell
python -m pytest tests/test_validate_mace.py -q
wsl -d Ubuntu-24.04 -- bash -lc "cd /mnt/c/Users/dodo/Documents/projects/pyrene-dimer-analyzer && ./bin/micromamba run -n pyrene-xtb python validate_mace.py --all-methods --output validation_results.csv"
```

Expected:
1. No xTB import failure.
2. Numeric GFN2-xTB rows in `validation_results.csv`.
3. MACE small pyrene error within pass threshold (`|error| < 0.2 A`).

## 3. Commit + Push to GitHub

### 3.1 Review Changes

```powershell
git status -sb
git diff -- . ':!*.csv' ':!*.png'
```

### 3.2 Stage Recommended Files

```powershell
git add `
  validate_mace.py `
  tests/test_validate_mace.py `
  reoptimize_moe_conformers.py `
  pyrene_analyzer/mace_optimizer.py `
  GITHUB_CLUSTER_RUNBOOK.md `
  docs/scientific-analysis/mace_validation.md `
  docs/COMPUTE_EXECUTION_PLAN.md `
  scripts/cluster_preflight.sh `
  scripts/run_dual_gpu_reopt.sh `
  scripts/merge_reopt_chunks.py `
  NEXT_CHAT_HANDOFF_PROMPT.md
```

Check staged set:

```powershell
git diff --cached --name-status
```

Commit:

```powershell
git commit -m "Finalize MACE validation wiring, add cluster run automation, and document compute migration"
```

Push:

```powershell
git push
```

## 4. Clone on the New Linux Server

```bash
git clone https://github.com/fulopjoz/pyrene-dimer-analyzer.git
cd pyrene-dimer-analyzer
```

## 5. Environment Setup on Server

```bash
micromamba create -n pyrene-xtb python=3.11 -y
micromamba activate pyrene-xtb
python -m pip install --upgrade pip
python -m pip install -e ".[mace]"
python -m pip install xtb-python ase
```

If your server policy requires manual torch install, install torch first, then:

```bash
python -m pip install -e .
python -m pip install mace-torch ase xtb-python
```

## 6. Server Preflight (No Scheduler)
Run the shipped preflight script:

```bash
bash scripts/cluster_preflight.sh
```

This checks:
1. Python + package imports (`torch`, `rdkit`, `ase`, `mace`, `xtb`).
2. GPU visibility and memory.
3. Required input files.

## 7. Smoke Tests on Server

```bash
python -m pytest tests/test_validate_mace.py -q
python validate_mace.py --all-methods --output validation_results.csv
python reoptimize_moe_conformers.py --test --model small --device cuda --mace-dtype float64 --max-steps 500 --fmax 0.005 --save-every 1 --output-prefix smoke
```

## 8. Production Runs (No Batch System Yet)

### 8.1 Dual-GPU Recommended Path

```bash
bash scripts/run_dual_gpu_reopt.sh \
  --model small \
  --mace-dtype float64 \
  --max-steps 500 \
  --fmax 0.005 \
  --save-every 25 \
  --output-dir runs \
  --output-prefix moe_mace_reopt
```

This automatically:
1. Counts conformers from the SDF.
2. Splits into two index ranges.
3. Launches one process per GPU (`CUDA_VISIBLE_DEVICES=<gpu>`).
4. Writes logs + per-GPU CSVs.
5. Merges chunk outputs into final `_all_conformers.csv` and `_summary.csv`.

### 8.2 Faster Throughput Mode (if acceptable)

```bash
bash scripts/run_dual_gpu_reopt.sh \
  --model small \
  --mace-dtype float32 \
  --max-steps 500 \
  --fmax 0.005 \
  --save-every 25 \
  --output-dir runs_f32 \
  --output-prefix moe_mace_reopt_f32
```

## 9. Manual Merge (if needed)
If you rerun chunks manually, merge them with:

```bash
python scripts/merge_reopt_chunks.py \
  --input-glob "runs/moe_mace_reopt_gpu*_all_conformers.csv" \
  --output-dir runs \
  --output-prefix moe_mace_reopt
```

## 10. Runtime Guidance
Based on local measurements in this repo:
1. CPU-only MACE (float64) is too slow for production.
2. GPU float64 is scientifically safest for geometry optimization.
3. GPU float32 is much faster and suitable for pre-screening.
4. Best practice:
   - Bulk screening: `small + float32 + dual GPU`
   - Final shortlist: rerun with `small/medium + float64`

## 11. Troubleshooting
1. `cuequivariance ... not available` is non-fatal; runs continue.
2. `TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD` warnings are non-fatal.
3. If one GPU process fails, inspect:
   - `runs/<prefix>_gpu0.log`
   - `runs/<prefix>_gpu1.log`
4. Resume failed chunk manually:

```bash
CUDA_VISIBLE_DEVICES=0 python reoptimize_moe_conformers.py \
  --resume \
  --start-index <start> \
  --end-index <end> \
  --output-prefix runs/<prefix>_gpu0 \
  --model small \
  --device cuda \
  --mace-dtype float64 \
  --max-steps 500 \
  --fmax 0.005 \
  --save-every 25
```

## 12. Next Phase: ORCA Excited-State Jobs
After MACE re-optimization is complete and merged:
1. Build ORCA shortlist and inputs.
2. Run TD-DFT queue locally (no scheduler).
3. Collect parsed excited-state summary.

Quick command:

```bash
bash scripts/run_orca_next_steps.sh
```

Manual steps are documented in:
- `docs/ORCA_EXCITED_STATE_WORKFLOW.md`

If `rg` is not installed on server, replace process checks with `grep`:

```bash
ps -ef | grep -E "run_dual_gpu_reopt|reoptimize_moe_conformers.py|orca_run_queue.py" | grep -v grep
```
