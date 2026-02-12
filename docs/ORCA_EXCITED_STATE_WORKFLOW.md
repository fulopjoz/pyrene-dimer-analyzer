# ORCA Excited-State Workflow (Post-MACE)

## Goal
Add excited-state photophysics on top of completed MACE geometry screening.

This phase estimates:
1. Vertical excitation energies and oscillator strengths.
2. S1-relaxed geometries for top candidates.
3. Optional vibronic absorption/fluorescence lineshapes with `ORCA_ESD`.

## Scope and Compute Strategy
1. Do **not** run ORCA on all 3347 conformers.
2. Use staged filtering:
   - MACE geometry run (broad screening).
   - ORCA TD-DFT on shortlist.
   - ORCA_ESD only on final few.

Recommended on 24 CPU server:
1. `nprocs=12` per ORCA job.
2. `max_workers=2` queue parallelism.

## Prerequisites
1. Complete and merge MACE outputs:
   - `runs/moe_mace_reopt_all_conformers.csv`
2. ORCA installed and visible in PATH:
   - `orca` command works
3. Python env contains RDKit + pandas.

## Files Added for This Workflow
1. Candidate selection:
   - `scripts/orca_select_candidates.py`
2. ORCA input preparation:
   - `scripts/orca_prepare_inputs.py`
   - templates in `orca_templates/`
3. Local queue execution (no scheduler):
   - `scripts/orca_run_queue.py`
4. Result collection:
   - `scripts/orca_collect_results.py`
5. End-to-end helper:
   - `scripts/run_orca_next_steps.sh`

## Scientific Defaults
1. Vertical and S1 optimization defaults:
   - `CAM-B3LYP/def2-SVP`
   - `D3BJ`
   - `CPCM(toluene)`
   - `TDA true`
2. S1 optimization uses:
   - `%tddft IRoot 1`
   - `FOLLOWIROOT TRUE` (root-tracking stability).

You can change these in CLI args or templates.

## Stage Definitions
1. `vertical`
   - S0-geometry TD-DFT vertical excitations.
   - Output focus: S1 energy (eV), oscillator strength `f`.
2. `s1opt`
   - Excited-state geometry optimization on target root.
3. `s0freq`, `s1freq` (optional)
   - Frequency jobs for vibronic analysis.
4. `esd_abs`, `esd_fluor` (optional)
   - ORCA_ESD spectra/rates using Hessians from freq jobs.

## Runbook Commands
From repo root:

```bash
# 1) Candidate shortlist from merged MACE output
python scripts/orca_select_candidates.py \
  --input-csv runs/moe_mace_reopt_all_conformers.csv \
  --output-csv runs/orca_candidates.csv \
  --top-per-molecule 2 \
  --max-total 96

# 2) Generate ORCA job folders and inputs
python scripts/orca_prepare_inputs.py \
  --candidates-csv runs/orca_candidates.csv \
  --sdf moe_conformers/cnph_th_cf3_3d_conformers.sdf \
  --output-dir runs/orca_jobs \
  --stages vertical,s1opt \
  --nprocs 12 \
  --maxcore 3000 \
  --functional CAM-B3LYP \
  --basis def2-SVP \
  --dispersion D3BJ \
  --solvent toluene \
  --nroots 10 \
  --iroot 1 \
  --tda

# 3) Run queue on local server (resume-safe)
python scripts/orca_run_queue.py \
  --manifest runs/orca_jobs/orca_job_manifest.csv \
  --orca-bin orca \
  --max-workers 2 \
  --resume \
  --summary-csv runs/orca_jobs/orca_run_summary.csv

# 4) Parse outputs into summary table
python scripts/orca_collect_results.py \
  --jobs-root runs/orca_jobs \
  --glob "**/*.out" \
  --output-csv runs/orca_jobs/orca_results_summary.csv
```

## One-Command Helper
```bash
# Uses defaults above (set ORCA_BIN/WORKERS/NPROCS_PER_JOB if needed)
bash scripts/run_orca_next_steps.sh
```

Environment overrides example:
```bash
export ORCA_BIN=/opt/orca/orca
export WORKERS=2
export NPROCS_PER_JOB=12
export DRY_RUN=1
bash scripts/run_orca_next_steps.sh
```

## Validation Checklist
1. `orca_job_manifest.csv` exists and has expected number of jobs.
2. `orca_run_summary.csv` has no `failed` rows.
3. `orca_results_summary.csv` has parsed S1 energies for vertical stage.
4. Spot-check 3 outputs for `ORCA TERMINATED NORMALLY`.

## Notes
1. The current phase is still model-driven screening; final experimental validation is separate.
2. If convergence/root issues appear, increase `nroots` and adjust `iroot`.
3. Use smaller shortlist first (for example `--max-total 24`) before scaling.

## Primary References
1. ORCA parallel execution:
   - https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html
2. ORCA TD-DFT excited states:
   - https://www.faccts.de/docs/orca/6.0/manual/contents/typical/excitedstates.html
3. ORCA TD-DFT root following (`IRoot`, `FOLLOWIROOT`):
   - https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/tddft.html
4. ORCA ESD module:
   - https://www.faccts.de/docs/orca/6.0/manual/contents/typical/esd.html
