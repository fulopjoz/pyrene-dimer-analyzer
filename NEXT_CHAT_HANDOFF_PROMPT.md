# Next Chat Handoff Prompt

Use this when starting the next chat/session on the new Linux server.

## Context Prompt (paste first)
```text
We are continuing work on pyrene-dimer-analyzer after fixing validate_mace.py wiring and preparing dual-GPU execution scripts.

Repository goals:
1) Validate MACE-OFF23 for aromatic pi-stacking geometry.
2) Run full 3347-conformer re-optimization efficiently.
3) Produce reproducible CSV outputs for downstream molecular ranking.

What is already done:
- validate_mace.py now uses optimize_with_gfn2xtb (xTB import issue fixed)
- benchmark starts at initial_distance=4.0 A
- tight optimizer settings in benchmark wrappers (fmax=0.005)
- regression tests in tests/test_validate_mace.py
- runbook: GITHUB_CLUSTER_RUNBOOK.md
- cluster scripts:
  - scripts/cluster_preflight.sh
  - scripts/run_dual_gpu_reopt.sh
  - scripts/merge_reopt_chunks.py
- validation evidence doc: docs/scientific-analysis/mace_validation.md
- compute migration doc: docs/COMPUTE_EXECUTION_PLAN.md

Please first inspect current repo state, run preflight + smoke tests, then execute the best production run plan for this machine and report bottlenecks with concrete timings.
```

## First Message to Send
```text
Start by running this checklist and summarize findings with exact commands and outputs:
1) git status -sb
2) bash scripts/cluster_preflight.sh
3) python -m pytest tests/test_validate_mace.py -q
4) python validate_mace.py --all-methods --output validation_results.csv
5) python reoptimize_moe_conformers.py --test --model small --device cuda --mace-dtype float64 --max-steps 500 --fmax 0.005 --save-every 1 --output-prefix smoke_server

Then propose the exact dual-GPU command for the full run on this server and estimate runtime from measured smoke timing.
```

## Follow-up Prompt (Post-MACE Completion)
```text
MACE chunk runs are complete and merged. Continue with ORCA excited-state phase using repo scripts:
1) python scripts/orca_select_candidates.py ...
2) python scripts/orca_prepare_inputs.py ...
3) python scripts/orca_run_queue.py --resume ...
4) python scripts/orca_collect_results.py ...

Provide a concise QC report:
- total jobs
- failed jobs
- how many vertical jobs returned parsed S1 energies and oscillator strengths
- list top candidates by lowest S1 energy and strongest fosc
```
