# MACE-OFF23 Validation for Pi-Stacking Distances

## Purpose

Before running the full 3,347-conformer re-optimization with MACE-OFF23, we validate
that the model reproduces known pi-stacking distances for aromatic dimers of increasing
size. This is critical because MACE-OFF23 was trained on the SPICE dataset (mostly
MW < 500), and our binaphthalene dimers (MW 800-1200) may represent extrapolation.

## Benchmark Molecules

| System | N atoms | Reference d (A) | Source |
| --- | --- | --- | --- |
| Benzene dimer | 24 | 3.40 | Jurecka et al. 2006 (S22, CCSD(T)/CBS) |
| Naphthalene dimer | 36 | 3.45 | Crystal structure data |
| Pyrene dimer | 52 | 3.43 | PMC11476719 (2024), DLPNO-CCSD(T)/CBS |

## Pass Criteria

**MACE-OFF23 pyrene dimer distance within 0.2 A of 3.43 A reference.**

- Pass (|error| < 0.2 A): Proceed with full re-optimization
- Marginal (0.2-0.3 A): Use as pre-screener, verify top candidates with GFN2-xTB
- Fail (|error| > 0.3 A): Reject MACE for geometry; use only for energy ranking

## Results

*To be filled after running `python validate_mace.py --all-methods`*

| System | MMFF94s (A) | MACE-OFF23 small (A) | MACE-OFF23 medium (A) | GFN2-xTB (A) | Reference (A) |
| --- | --- | --- | --- | --- | --- |
| Benzene dimer | - | - | - | - | 3.40 |
| Naphthalene dimer | - | - | - | - | 3.45 |
| Pyrene dimer | - | - | - | - | 3.43 |

### Error Analysis

*To be computed from results.*

| Method | MAE (A) | Max error (A) | Passes? |
| --- | --- | --- | --- |
| MMFF94s | - | - | Expected: NO (lacks dispersion) |
| MACE-OFF23 (small) | - | - | - |
| MACE-OFF23 (medium) | - | - | - |
| GFN2-xTB | - | - | Expected: YES (D4 dispersion) |

## Decision

*To be filled after validation.*

- [ ] Full adoption: MACE-OFF23 validated, proceed with `reoptimize_moe_conformers.py`
- [ ] Pre-screener only: MACE for initial filter, GFN2-xTB for top candidates
- [ ] Rejected: MACE not suitable for these systems

## How to Run

```bash
# Requires MACE-OFF23 installed:
pip install mace-torch ase torch

# Run validation benchmark:
python validate_mace.py --output validation_results.csv

# Include GFN2-xTB comparison (requires WSL):
python validate_mace.py --all-methods --output validation_results.csv
```

## References

- Batatia et al. (2024) arXiv:2401.00096 (MACE-OFF23 architecture and training)
- Jurecka et al. (2006) PCCP 8, 1985 (S22 benchmark set)
- Bannwarth et al. (2019) JCTC 15, 1652 (GFN2-xTB with D4 dispersion)
- PMC11476719 (2024) DLPNO-CCSD(T)/CBS pyrene dimer equilibrium distance
- Najibi & Goerigk (2020) JCTC 16, 4479 (wB97M-D3BJ accuracy for NCI)
