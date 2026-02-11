# MACE-OFF23 Validation for Pi-Stacking Distances

## Purpose
Validate that the benchmark configuration in `validate_mace.py` is physically meaningful and that MACE-OFF23 reproduces known aromatic dimer stacking distances before running full re-optimization on 3,347 conformers.

## Benchmark Systems

| System | N atoms | Reference d (A) | Reference |
| --- | ---: | ---: | --- |
| Benzene dimer | 24 | 3.40 | Jurecka et al. (2006), S22 CCSD(T)/CBS |
| Naphthalene dimer | 36 | 3.45 | Crystal stacking distance |
| Pyrene dimer | 52 | 3.43 | PMC11476719 (2024), DLPNO-CCSD(T)/CBS |

## Pass Criterion
Primary gate:

`|error| < 0.2 A` for `method == MACE-OFF23` on `benchmark == pyrene_dimer`.

Interpretation bands:
1. Pass: `< 0.2 A` -> proceed with MACE small for production runs.
2. Marginal: `0.2 - 0.3 A` -> pre-screener only; verify shortlist with xTB.
3. Fail: `> 0.3 A` -> do not use for geometry optimization.

## Run Configuration Used
Command:

```bash
python validate_mace.py --all-methods --output validation_results.csv
```

Key settings now active in code:
1. Initial dimer separation is `4.0 A` (not `5.0 A`).
2. MMFF uses `MMFF94s` with `ignoreInterfragInteractions=False`.
3. MACE and xTB wrappers use tighter convergence (`fmax=0.005`) and larger step caps.
4. xTB wrapper calls `optimize_with_gfn2xtb` directly.

## Results (from `validation_results.csv`)

| System | MMFF94s (A) | MACE small (A) | MACE medium (A) | GFN2-xTB (A) | Reference (A) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Benzene dimer | 5.1680 | 3.4811 | 3.6228 | 3.4483 | 3.40 |
| Naphthalene dimer | 4.9657 | 3.4590 | 3.6293 | 3.3918 | 3.45 |
| Pyrene dimer | 4.4977 | 3.3840 | 3.6316 | 3.3565 | 3.43 |

## Error Summary

| Method | MAE (A) | Max abs error (A) | Pyrene error (A) | Pyrene pass? |
| --- | ---: | ---: | ---: | --- |
| MMFF94s | 1.4505 | 1.7680 | +1.0677 | No |
| MACE-OFF23 (small) | 0.0454 | 0.0811 | -0.0460 | **Yes** |
| MACE-OFF23 (medium) | 0.2012 | 0.2228 | +0.2016 | Borderline/No |
| GFN2-xTB | 0.0600 | 0.0735 | -0.0735 | Yes |

## Decision
1. `MACE-OFF23 (small)` is **validated** for this benchmark and meets the defined pyrene pass criterion.
2. `MACE-OFF23-medium` is slightly above the pyrene threshold in this setup and should be treated as informational unless rerun with stricter validation protocol.
3. `GFN2-xTB` remains suitable as a high-confidence verifier for top-ranked candidates.
4. `MMFF94s` is unsuitable for pi-stacking geometry benchmarking in disconnected dimers.

## Reproducibility Checks

```bash
python -m pytest tests/test_validate_mace.py -q
python validate_mace.py --all-methods --output validation_results.csv
```

The regression tests confirm:
1. Correct xTB API wiring (`optimize_with_gfn2xtb`).
2. Correct MMFF inter-fragment settings.
3. Benchmark initialization at `initial_distance=4.0`.
4. Unchanged pyrene pass/fail threshold logic.

## References
1. Batatia et al. (2024), MACE-OFF23 foundation model, arXiv:2401.00096  
   https://arxiv.org/abs/2401.00096
2. Jurecka et al. (2006), S22 noncovalent benchmark set, PCCP 8, 1985.
3. Bannwarth et al. (2019), GFN2-xTB, JCTC 15, 1652-1671  
   https://doi.org/10.1021/acs.jctc.8b01176
4. Pyrene dimer benchmark source: PMC11476719 (2024), DLPNO-CCSD(T)/CBS  
   https://pmc.ncbi.nlm.nih.gov/articles/PMC11476719/
