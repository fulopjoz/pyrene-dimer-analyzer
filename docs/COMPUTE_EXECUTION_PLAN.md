# Compute Execution Plan (Local vs Server)

## Goal
Run the aromatic dimer pipeline with correct physics and practical runtime by splitting work between:
1. Local workstation (development, validation, smoke checks).
2. Stronger Linux server (full production optimization on 24 CPU + 2 GPUs).

## Current Baseline (Measured in this repo)

### Validation benchmark
Command:

```bash
python validate_mace.py --all-methods --output validation_results.csv
```

Key outcome:
1. `MACE-OFF23` small on pyrene dimer: `error = -0.0460 A` (passes `|error| < 0.2 A`).
2. xTB import-path failure is fixed (numeric GFN2-xTB rows produced).

### Local re-optimization smoke (3 conformers, same SDF pipeline)

1. `cuda + float64`:
   - Runtime reported by script: ~205 s for 3 conformers.
   - Wall time: ~230 s.
2. `cuda + float32`:
   - Runtime reported by script: ~63 s for 3 conformers.
   - Wall time: ~79 s.
3. `cpu + float64`:
   - Exceeded 10-minute timeout for the same 3-conformer smoke run.

Interpretation:
1. CPU-only path is not suitable for full 3,347-conformer production.
2. GPU path is required.
3. Float32 gives strong throughput gains and is suitable for pre-screening.
4. Float64 remains preferred for final geometry accuracy.

## Hardware-Aware Strategy

## What to run locally
1. Unit tests and wiring checks:
   - `python -m pytest tests/test_validate_mace.py -q`
2. Benchmark validation:
   - `python validate_mace.py --all-methods --output validation_results.csv`
3. Small smoke jobs:
   - `python reoptimize_moe_conformers.py --test ...`

## What to run on the stronger server
1. Full 3,347-conformer MACE re-optimization.
2. Medium-model reruns on shortlisted candidates.
3. Optional xTB verification on top-ranked molecules.

Reason:
1. Full production jobs are throughput-bound and benefit directly from dual-GPU concurrency.
2. Server has enough CPU to support I/O, RDKit preprocessing, and two GPU worker processes.

## Recommended Production Modes

## Mode A: Accuracy-first
Use when final geometry fidelity matters most.

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

## Mode B: Throughput-first
Use for rapid ranking/pre-screening before final verification.

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

## Expected Runtime Scaling (order-of-magnitude)
Using local measured per-conformer timings as baseline:
1. Single GPU float64: roughly multi-day scale for full run.
2. Single GPU float32: substantially faster (roughly 2x-3x speedup observed locally).
3. Dual GPU: near 2x further speedup if both GPUs are saturated.

Exact server timing depends on GPU model, driver, and thermal limits. Run a `--test` smoke job on server first and extrapolate from actual wall time.

## Verification and Acceptance Criteria
1. `validate_mace.py` completes with no xTB import errors.
2. `validation_results.csv` has numeric rows for all selected methods.
3. Full run creates:
   - `<prefix>_all_conformers.csv`
   - `<prefix>_summary.csv`
   - optional Boltzmann CSVs
4. No unrecoverable chunk failures in GPU logs.

## Research Basis / Method Rationale
1. MACE-OFF23: trained on DFT-level data suitable for noncovalent interactions (Batatia et al., 2024).
2. xTB (GFN2-xTB): practical semi-empirical verifier for pi-stacking-sensitive geometries (Bannwarth et al., 2019).
3. Pyrene dimer reference target: DLPNO-CCSD(T)/CBS benchmark from PMC11476719 (2024).
4. Validation benchmark references are documented in:
   - `docs/scientific-analysis/mace_validation.md`

## Practical Next Step on Server
1. Clone repo.
2. Set up environment.
3. Run preflight:
   - `bash scripts/cluster_preflight.sh`
4. Run dual-GPU production:
   - `bash scripts/run_dual_gpu_reopt.sh ...`
5. Review logs and merged outputs in `runs/`.
