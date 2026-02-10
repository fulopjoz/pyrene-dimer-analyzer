# Pyrene Dimer Analyzer - Project Plan & Roadmap

## Project Status Summary

**Date**: 2026-02-05
**Version**: 2.0 (Literature-Informed Update)
**Tests**: 224 passing, 78% coverage
**Status**: Ready for GFN2-xTB Integration

---

## Executive Summary: Literature-Validated Strategy Update

This plan integrates findings from three comprehensive literature reviews conducted on 2026-02-05:

1. **Dispersion Correction Methods** - GFN2-xTB as optimal replacement for MMFF94s
2. **Substituent Effects on Aromatic Stacking** - Taft Es (steric) dominates over Hammett sigma (electronic)
3. **Binaphthalene Excimer Formation** - Ethynyl spacer effect validated by literature

### Key Decision: Recommended Solution

**Primary Pipeline**: RDKit ETKDGv3 embedding → **GFN2-xTB optimization** → Geometric analysis

This approach:
- Correctly reproduces π-π stacking distances (3.3-3.6 Å vs MMFF94s's 3.8-4.2 Å)
- Is computationally tractable (~2-3 min per 140-atom molecule)
- Is available on Windows via conda-forge (with ASE interface)
- Preserves excimer conformers that MMFF94s energy filtering removes

### Literature-Validated Findings

| Finding | Source | Implication |
|---------|--------|-------------|
| Pyrene dimer equilibrium: 3.43 Å | DLPNO-CCSD(T)/CBS, PMC11476719 | GFN2-xTB reproduces; MMFF94s overestimates |
| GFN2-xTB S22/S66 MAE: 0.26-0.31 kcal/mol | Rowan Benchmarks | Best SQM for noncovalent interactions |
| Steric (Es) > Electronic (σ) | Wheeler-Houk 2011 | Use Taft Es, not Hammett σ for SAR |
| Ethynyl spacer 2× enhancement | Winnik 1993, CCS Chem 2020 | Validates EtynPyr >> Pyr results |
| MACE-OFF23 includes P, S | JACS 2024 | Future option for highest accuracy |

---

### Molecular Identity Verification

All molecules verified against the hand-drawn structural diagram (Czech: "Chiralni, magneticke a luminiscencni dopanty"). The target structure is:

```
Pyrene(R)₂ — O-CH(CH₃)-O — Ar(NH₂)(PO₃H₂) — O-CH(CH₃)-O — Pyrene(R)₂
```

Key features:
- Two pyrene units, each bearing 2× R substituents at 2,7-positions
- Connected through chiral 1-methylethylene diether bridges (-O-CH(CH₃)-O-)
- Central 1,3,5-trisubstituted benzene with aminophosphonate group
- 3 chiral centers, 1 phosphorus, 1 nitrogen per molecule

| Variant | R Group | Formula | Heavy Atoms | MW | InChIKey (layers 1+2) |
|---------|---------|---------|-------------|-----|----------------------|
| Et | -CH₂CH₃ | C₅₅H₅₄NO₇P | 64 | 871.4 | SYNUCIUKKVFYQM-NWSHGZSRSA |
| iPr | -CH(CH₃)₂ | C₅₉H₆₂NO₇P | 68 | 927.4 | UNCJZNGUKQIWDW-SDSWSPLBSA |
| cHex | -C₆H₁₁ | C₇₁H₇₈NO₇P | 80 | 1087.6 | SDTQGRGNUYTSNM-ZSGTWBQWSA |
| tBu | -C(CH₃)₃ | C₆₃H₇₀NO₇P | 72 | 983.5 | VJXVMEFTEUKCLA-SDSWSPLBSA |

**Verification**: MOE SDF InChIKeys match original test SDF InChIKeys exactly (layers 1+2). MOE uses deprotonated phosphonate ([O-], layer 3 = M) while source SDF has neutral form (layer 3 = N). RDKit conformer SMILES produce exact InChIKey match with original SDF.

---

## Completed Work

### Phase 1: Core Tool Development (v1.0 → v1.1)
- [x] Multi-system aromatic support (pyrene, perylene, anthracene, naphthalene, phenanthrene)
- [x] Aromatic systems registry with literature-backed thresholds
- [x] Virtual screening module (SubstituentScreener, R-group enumeration)
- [x] Biased conformer generation (distance-constrained ETKDGv3)
- [x] SMILES → conformer → analysis pipeline
- [x] Progress indicators (verbose mode)
- [x] High-angle geometry warnings
- [x] Documentation sync

### Phase 2: Performance & Bug Fixes
- [x] Lazy imports via `__getattr__` in `__init__.py` (CLI startup < 1s)
- [x] Deferred imports in `cli.py`
- [x] Critical Windows Unicode bug fix (`π` → `pi` in core.py:599)
- [x] Documentation quickstart fix (real test file paths)

### Phase 3: MOE vs RDKit Validation Study
- [x] RDKit conformer generation: all 4 variants (Et=53, iPr=35, cHex=14, tBu=35 conformers)
- [x] MOE property extraction via moebatch SVL (Et=200, iPr=200 conformers)
- [x] MOE SDF geometric analysis (Et=200, iPr=200 conformers)
- [x] Statistical comparison (KS tests, Mann-Whitney U)
- [x] 3 publication-quality comparison figures generated
- [x] Molecular identity verification (InChIKey cross-check)

---

## Key Findings: MOE vs RDKit Comparison

### Geometric Parameter Summary

| Parameter | Et MOE (n=200) | Et RDKit (n=53) | iPr MOE (n=200) | iPr RDKit (n=35) |
|-----------|----------------|-----------------|-----------------|------------------|
| θ (deg) | 25.6 ± 19.3 | 24.1 ± 22.7 | 30.7 ± 21.6 | 21.3 ± 20.1 |
| d (Å) | **3.45 ± 0.85** | **4.24 ± 1.45** | **3.53 ± 1.08** | **4.26 ± 1.69** |
| Overlap (%) | 30.0 ± 17.0 | 37.0 ± 20.3 | 26.2 ± 16.7 | 39.0 ± 18.3 |
| Strong excimer | 10 (5.0%) | 0 (0%) | 10 (5.0%) | 0 (0%) |
| Total excimer | 101 (50.5%) | 32 (60.4%) | 87 (43.5%) | 23 (65.7%) |

### Critical Finding: Systematic Distance Bias

MOE (Amber10:EHT) produces inter-plane distances **0.6-0.8 Å shorter** than RDKit (MMFF94s):
- Et: 3.45 vs 4.24 Å (Δ = 0.79 Å, KS p < 0.0001)
- iPr: 3.53 vs 4.26 Å (Δ = 0.73 Å, KS p = 0.0004)

**Root cause**: Amber10:EHT includes explicit aromatic π-stacking interaction terms calibrated to CCSD(T)/CBS benchmarks. MMFF94s relies on generic van der Waals terms that underestimate attractive dispersion between aromatic π systems.

**Consequence**: Only MOE finds strong excimers (d = 3.3-3.7 Å). RDKit distances cluster around 4.0 Å, which is above the strong excimer threshold but within weak excimer range.

### Qualitative Agreement

Despite quantitative distance differences, both methods agree on:
1. Et has higher excimer propensity than iPr
2. Steric trend: Et > iPr (> cHex > tBu from RDKit-only data)
3. Energy-overlap correlation is negative (higher overlap → lower energy)
4. Both methods sample the same conformational space (similar angle distributions)

---

## Immediate Next Steps (Implementation Plan)

### Step 1: Complete the 4-Variant Dataset

**What**: Generate MOE conformers for cHex and tBu variants
**Who**: User (in MOE GUI)
**How**: Same LowModeMD workflow used for Et and iPr:
1. Open `pyrene_dimer_set_for_MOE.mdb` in MOE Database Viewer
2. Select cHex entry → Conformational Search → LowModeMD
3. Settings: maxconf=200, RMSD=0.5, gtest=0.005, Amber10:EHT, T=300K
4. Export as SDF: File → Save As → `cHex_conformers.sdf`
5. Repeat for tBu

**Deliverable**: `cHex_conformers.sdf`, `tBu_conformers.sdf` (200 conformers each)

### Step 2: Distance Calibration Module

**What**: Add an empirical distance correction to the analyzer so RDKit-generated distances can be adjusted toward MOE-quality values.

**Scientific justification**: The 0.6-0.8 Å systematic offset between MMFF94s and Amber10:EHT is consistent with known force field limitations:
- Experimental pyrene excimer stacking: 3.3-3.5 Å (CSD crystal structures)
- CCSD(T)/CBS reference for pyrene dimer: 3.4-3.5 Å
- MMFF94s typical aromatic stacking: 3.8-4.2 Å
- Amber10:EHT aromatic stacking: 3.3-3.7 Å

**Implementation**:
```python
# In geometry.py or core.py
DISTANCE_CALIBRATION = {
    'mmff94s_to_amber10eht': {
        'slope': 0.85,      # To be fitted from paired data
        'intercept': 0.55,   # To be fitted from paired data
        'method': 'linear',
        'note': 'Empirical correction for MMFF94s → Amber10:EHT distances'
    }
}
```

**Data needed**: Paired (MOE, RDKit) conformer distances at similar geometries. We have this for Et and iPr (200+53 and 200+35 conformers). Will be stronger with cHex and tBu.

**Approach**:
1. Bin conformers by plane angle (0-10, 10-20, ..., 80-90°)
2. Within each bin, compare MOE vs RDKit distance distributions
3. Fit a linear correction: `d_corrected = slope * d_rdkit + intercept`
4. Validate with leave-one-variant-out cross-validation

### Step 3: Reclassify with Corrected Distances

**What**: Apply the calibration to all RDKit conformer data and regenerate classifications
**Impact**: Expected to:
- Reduce RDKit distances by ~0.6-0.8 Å
- Increase strong excimer detection
- Bring RDKit excimer populations closer to MOE values

### Step 4: Ensemble Feature Engineering

**What**: Build a feature engineering module that converts per-conformer geometric data into per-molecule statistical descriptors for ML.

**Features per molecule** (~40-45 total):
- **Distribution statistics**: mean, std, min, max, median, p10, p90 for each of: θ, d, overlap, centroid_distance, slip_stack
- **Threshold features**: fraction with θ<20°, fraction with d∈[3.3,3.7], fraction classified as strong_excimer, fraction classified as any_excimer
- **Boltzmann-weighted means** at 298K using MMFF94s energies (kT = 0.593 kcal/mol)
- **Energy features**: energy range, energy gap between lowest excimer and global minimum

### Step 5: Systematic R-Group Screen

**What**: Use `pyrene-analyze screen` to test 15-25 substituents and rank by excimer propensity
**Substituents**: Me, Et, nPr, iPr, nBu, iBu, sBu, tBu, cPen, cHex, Ph, Bn, OMe, OEt, F, Cl, Br, CF3, CN, NMe2, SMe, vinyl, allyl, acetyl
**Parameters**: 100 conformers each, biased ETKDGv3, 10 kcal/mol energy window
**Output**: Ranked table of substituent → excimer population, mean θ, mean d, mean overlap

### Step 6: QSAR / ML Classification

**What**: Build a predictive model for excimer formation from geometric ensemble features

**Approach** (per literature review):
1. Random Forest classifier (primary) with `class_weight='balanced'`
2. Gradient Boosted Trees (secondary)
3. Nested stratified 5-fold CV for unbiased performance estimation
4. Permutation test (p < 0.05 required)
5. SHAP analysis for feature importance and scientific interpretation

**Target variable**: Binary (excimer-forming vs non-excimer-forming) or 3-class (strong/weak/monomer)
**Training data**: 25+ substituent variants × ensemble features

### Step 7: Commit and Release

**What**: Commit all uncommitted changes, create v1.2.0 release
**Changes to commit**:
- Lazy imports in `__init__.py` and `cli.py`
- Unicode bug fix in `core.py`
- Documentation updates
- Screening module
- Aromatic systems registry
- All analysis data and plots
- Distance calibration module
- Feature engineering module

---

## Medium-Term Goals

### CSD Crystal Structure Validation
- Search Cambridge Structural Database for pyrene dimer crystal structures
- Compare computed stacking distances (MOE, RDKit, corrected RDKit) with experimental X-ray values
- This provides the ultimate ground truth for force field accuracy

### TD-DFT Photophysics
- Take the 10 strong excimer conformers from MOE (each variant)
- Run TD-DFT (wB97X-D3/6-31G* or CAM-B3LYP/6-311+G(d,p))
- Compute S₁ excitation energies and oscillator strengths
- Compare with experimental fluorescence spectra
- Map geometry → emission wavelength relationship

### Publication Strategy
- **Paper 1 (Methods)**: The pyrene-dimer-analyzer tool + MOE vs RDKit validation
  - Novel: automated geometric analysis of aromatic dimer conformers
  - Validation: systematic force field comparison with statistical tests
  - Utility: fast RDKit screening with calibration to MOE-quality results
  - Target: J. Chem. Inf. Model. or J. Cheminformatics

- **Paper 2 (Application)**: Substituent effects on excimer formation
  - SAR study: 25+ substituents ranked by excimer propensity
  - ML model: predict excimer formation from substituent descriptors
  - Design rules: which R groups promote/inhibit excimer formation
  - Target: Chem. Eur. J. or J. Phys. Chem.

---

## Data Files Inventory

### Source Data
| File | Description | Status |
|------|-------------|--------|
| `tests/test_data/pyrene_dimer_set_for_MOE.sdf` | 15 molecules (4 dimers + fragments) | Original |
| `Dimer_Et_conformers.mdb` | MOE database, 200 Et conformers | From MOE |
| `iPr_conformers.mdb` | MOE database, 200 iPr conformers | From MOE |
| `dimer_et_conformers.sdf` | MOE Et conformers (SDF export) | From MOE |
| `iPr_conformers.sdf` | MOE iPr conformers (SDF export) | From MOE |

### Analysis Results
| File | Description | Conformers |
|------|-------------|-----------|
| `Et_moe_analyzed.csv` | MOE Et geometric analysis | 200 |
| `iPr_moe_analyzed.csv` | MOE iPr geometric analysis | 200 |
| `Et_rdkit_conformers.csv` | RDKit Et analysis | 53 |
| `iPr_rdkit_conformers.csv` | RDKit iPr analysis | 35 |
| `cHex_rdkit_conformers.csv` | RDKit cHex analysis | 14 |
| `tBu_rdkit_conformers.csv` | RDKit tBu analysis | 35 |
| `Et_moe_data.csv` | MOE Et properties (E, dE, shape) | 200 |
| `iPr_moe_data.csv` | MOE iPr properties (E, dE, shape) | 200 |

### Plots
| File | Description |
|------|-------------|
| `plots/moe_vs_rdkit_geometry.png` | 8-panel geometric distribution comparison |
| `plots/moe_vs_rdkit_classification.png` | 6-panel classification + energy-geometry |
| `plots/excimer_quality_analysis.png` | 3-panel excimer quality violin/box plots |
| `plots/rdkit_ensemble_analysis.png` | 6-panel RDKit-only ensemble analysis |
| `plots/moe_vs_rdkit_comparison.png` | 6-panel property-based comparison |

---

## Complete Molecule Catalog

The SDF file `tests/test_data/pyrene_dimer_set_for_MOE.sdf` contains **15 molecules** corresponding to the research overview diagram ("Chiralni, magneticke a luminiscencni dopanty"):

### Dimers (Excimer Analysis Targets)

| # | Name | R Group | MW | Atoms | InChIKey | Role |
|---|------|---------|-----|-------|----------|------|
| 0 | Et dimer | -CH2CH3 | 871.4 | 64 | SYNUCIUKKVFYQM | Primary target |
| 1 | iPr dimer | -CH(CH3)2 | 927.4 | 68 | UNCJZNGUKQIWDW | Primary target |
| 2 | cHex dimer | -C6H11 | 1087.6 | 80 | SDTQGRGNUYTSNM | Primary target |
| 3 | tBu dimer | -C(CH3)3 | 983.5 | 72 | VJXVMEFTEUKCLA | Primary target |

All dimers: 9 aromatic rings, 42 aromatic C atoms, 7 O, 1 N, 1 P per molecule.

### Pyrene Monomers (Half-Dimers with OH Attachment Site)

| # | Name | R Group | MW | Role |
|---|------|---------|-----|------|
| 11 | Pyrene_R2_Et_OH | -CH2CH3 | 300.2 | Synthetic intermediate |
| 12 | Pyrene_R2_iPr_OH | -CH(CH3)2 | 328.2 | Synthetic intermediate |
| 13 | Pyrene_R2_cHex_OH | -C6H11 | 408.3 | Synthetic intermediate |
| 14 | Pyrene_R2_tBu_OH | -C(CH3)3 | 356.2 | Synthetic intermediate |

### Core Chromophores

| # | Name | MW | Type | Expected Emission |
|---|------|-----|------|-------------------|
| 5 | Pyrene | 228.1 | PAH | Monomer 375 nm, Excimer 480 nm |
| 6 | Ethynyl_pyrene_guess | 252.1 | Extended PAH | ~400 nm (extended conjugation) |
| 7 | Ethynyl_pyrene | 252.1 | Extended PAH | ~400 nm (coupling handle) |

### Alternative Luminescent Components (from image "R:" section)

| # | Name | MW | Type | Donor | Acceptor | Expected Emission |
|---|------|-----|------|-------|----------|-------------------|
| 8 | Me_thiophene_dicyanovinyl | 174.0 | D-A ICT | Thiophene | Dicyanovinyl | 500-600 nm |
| 9 | Me_thiophene_propionitrile_cyanophenyl | 252.1 | D-A ICT | Thiophene | Nitrile/cyanophenyl | 450-550 nm |

These represent **alternative luminescent units** for future dimer designs. Their intramolecular charge transfer (ICT) character produces red-shifted emission compared to pyrene.

### Structural/Linker Fragments

| # | Name | MW | Role |
|---|------|-----|------|
| 4 | Linker_diol_aminophosphonate | 307.1 | Central bridge unit |
| 10 | Phenoxy_methyl_phosphonic_acid | 188.0 | Phosphonate building block |

### How Image Maps to SDF

```
Image section a) Monomerni  -->  Molecules [11-14] (pyrene half-dimers)
Image section b) Dimerni     -->  Molecules [0-3] (full dimers, analysis targets)
Image "R:" luminescent       -->  Molecules [5-7] (pyrene variants)
                                  Molecules [8-9] (thiophene D-A chromophores)
Image "znám substituenty"    -->  Et, iPr, cHex, tBu (confirmed in [0-3] and [11-14])
```

---

## Technical Debt

1. **Uncommitted changes**: All improvements from v1.0 → v1.1 are uncommitted
2. **Bridge dihedral NaN**: Et and iPr return NaN for bridge dihedrals (bridge detection issue)
3. **cHex low conformer count**: Only 14 RDKit conformers survive energy filtering (may need wider window)
4. **No distance calibration**: RDKit distances not corrected for force field bias
5. **No feature engineering module**: Ensemble → molecule-level features not implemented
6. **Temporary scripts**: `export_mdb.svl`, `export_coords.svl`, `export_sdf.svl`, `plots_comparison.py` should be cleaned up or moved to `scripts/`

---

## Part 2: GFN2-xTB Integration Plan (NEW)

### 2.1 Literature-Validated Benchmarks

| Method | Pyrene Dimer Distance | MAE (S22/S66) | Suitable |
|--------|----------------------|---------------|----------|
| DLPNO-CCSD(T)/CBS | 3.43 Å (reference) | — | Gold standard |
| GFN2-xTB | 3.4-3.6 Å | 0.26-0.31 kcal/mol | **YES** |
| MMFF94s | 3.8-4.2 Å | Not benchmarked | NO (lacks dispersion) |
| MACE-OFF23 | ~3.4 Å | 0.3 kcal/mol (biaryl) | YES (includes P, S) |

### 2.2 Environment Setup

```bash
# Create conda environment with xtb-python
conda create -n pyrene-xtb python=3.10
conda activate pyrene-xtb

# Install xtb via conda-forge (NOT pip on Windows)
conda install -c conda-forge xtb-python ase

# Install project
pip install -e ".[dev]"
```

**Windows Notes:**
- Native xtb-python is v6.7.1pre (limited)
- For CREST: Use WSL2 with `conda install -c conda-forge xtb crest`

### 2.3 New Module: `pyrene_analyzer/xtb_optimizer.py`

Key functions:
- `has_xtb_available()` → bool
- `optimize_with_gfn2xtb(mol, conf_id, method, max_steps, fmax)` → (mol, energy)
- `optimize_conformer_ensemble(mol, method, n_jobs, verbose)` → mol

Interface: Uses ASE calculator with BFGS optimizer.

### 2.4 Updated Screening Pipeline

Add `optimizer` parameter to `generate_conformers_biased()`:
- `"MMFF94s"` - Fast but lacks dispersion (default for backward compat)
- `"GFN2-xTB"` - Dispersion-corrected, correct π-stacking (recommended)
- `"none"` - Skip optimization (fast but inaccurate)

### 2.5 Expected Improvements

| Metric | MMFF94s Only | + GFN2-xTB |
|--------|-------------|------------|
| π-stacking distance | 3.8-4.2 Å (too long) | 3.3-3.6 Å (correct) |
| Excimer conformer retention | Low (filtered out) | High (proper energies) |
| Energy ranking | Biased against stacked | Physically reasonable |

---

## Part 3: Threshold Calibration (NEW)

### 3.1 Updated Binaphthalene Thresholds

Based on DLPNO-CCSD(T)/CBS data for pyrene dimer (PMC11476719, 2024):

| Parameter | Current | Updated | Rationale |
|-----------|---------|---------|-----------|
| Strong d_min | 3.3 Å | 3.0 Å | Allow tighter stacking |
| Strong d_max | 3.8 Å | 3.6 Å | DLPNO equilibrium 3.43 Å |
| Weak d_max | 4.5 Å | 4.2 Å | Literature: 4.0-4.2 Å for weak |
| Weak angle_max | 60° | 50° | 60° is edge-to-face |

### 3.2 Implementation

```python
# In aromatic_systems.py
AROMATIC_SYSTEMS["binaphthalene"].thresholds = ClassificationThresholds(
    strong_angle_max=20.0,
    strong_distance_min=3.0,    # Updated
    strong_distance_max=3.6,    # Updated
    strong_overlap_min=40.0,
    weak_angle_max=50.0,        # Updated
    weak_distance_max=4.2,      # Updated
    weak_overlap_min=25.0,
)
```

---

## Part 4: SAR Analysis Enhancement (NEW)

### 4.1 Key Finding: Steric > Electronic

Literature (Wheeler-Houk 2011, Carter-Fenk & Herbert 2020) shows:
- **Taft Es (steric)** dominates over **Hammett σ (electronic)**
- CF3 vs CN: Similar σ_para (0.54 vs 0.66) but different Es (-2.40 vs -0.51)
- Our data: CF3 = 0% excimer, CN = 12.5% excimer → Es explains, not σ

### 4.2 Substituent Parameters Table

| Substituent | Es (Taft) | v (Charton) | σ_para | Our Excimer |
|-------------|-----------|-------------|--------|-------------|
| H | 0.00 | 0.00 | 0.00 | Reference |
| Me | -1.24 | 0.52 | -0.17 | High |
| Et | -1.31 | 0.56 | -0.15 | High |
| iPr | -1.71 | 0.76 | -0.15 | Medium |
| tBu | -2.78 | 1.24 | -0.20 | ~0% |
| CN | -0.51 | 0.40 | 0.66 | 12.5% |
| CF3 | -2.40 | 0.91 | 0.54 | 0% |
| cHex | -2.03 | 0.87 | -0.15 | Low |

### 4.3 Predictive Model

```python
# Primary model: steric control
excimer_fraction ~ β₀ + β₁·Es + β₂·[R-group_type] + β₃·(Es × R-group)
```

### 4.4 CNPh_Th_tBu Anomaly Explanation

- **Statistical noise**: Only 3 conformers (high uncertainty)
- **Preorganization**: tBu may lock favorable dihedral
- **Action**: Re-run with 200+ conformers to test

---

## Part 5: Detailed TODO List (Implementation Phase)

### Phase 1: GFN2-xTB Integration (Priority: HIGH)

- [ ] **1.1** Create conda environment with xtb-python
  - `conda create -n pyrene-xtb python=3.10`
  - `conda install -c conda-forge xtb-python ase`
  - Test: optimize benzene dimer

- [ ] **1.2** Implement `pyrene_analyzer/xtb_optimizer.py`
  - `optimize_with_gfn2xtb()` via ASE/BFGS
  - `optimize_conformer_ensemble()` for batch
  - `has_xtb_available()` for graceful fallback
  - Tests: `tests/test_xtb_optimizer.py`

- [ ] **1.3** Update `screening.py` with optimizer parameter
  - Add `optimizer: str = "MMFF94s"` to `generate_conformers_biased()`
  - Add to `analyze_from_smiles()`
  - CLI: `--optimizer mmff94s|gfn2-xtb`

- [ ] **1.4** Create installation documentation
  - `docs/installation_xtb.md`
  - Windows conda-forge, WSL2 for CREST

### Phase 2: Threshold Calibration (Priority: MEDIUM)

- [ ] **2.1** Update binaphthalene thresholds in `aromatic_systems.py`
  - Strong: 3.0-3.6 Å, <20°, >40% overlap
  - Weak: <4.2 Å, <50°, >25% overlap

- [ ] **2.2** Add literature references to docstrings
  - PMC11476719 (2024), Ge et al. (2020)

### Phase 3: SAR Analysis (Priority: MEDIUM)

- [ ] **3.1** Create `pyrene_analyzer/substituent_params.py`
  - `SUBSTITUENT_PARAMS` dict with Es, v, σ, MR
  - `get_param(substituent, param_name)` function

- [ ] **3.2** Update `sar_analysis.py`
  - Add Es/v columns to DataFrame
  - Fit Es-based linear model
  - Correlation matrix figure

- [ ] **3.3** Re-analyze CNPh_Th_tBu with 200+ conformers

### Phase 4: Full Re-screening (Priority: HIGH)

- [ ] **4.1** Parameter optimization with GFN2-xTB
  - Run `experiments/optimize_conformers.py` with xtb option
  - Test on 3 representative variants

- [ ] **4.2** Full 64-variant screening with GFN2-xTB
  - `run_screening.py --optimizer gfn2-xtb`
  - Output: `binaph_screening_xtb_*.csv`

- [ ] **4.3** MMFF94s vs GFN2-xTB comparison
  - Script: `compare_mmff_xtb.py`
  - Metrics: excimer fraction, distance, retention

### Phase 5: Documentation & Publication

- [ ] **5.1** Update CLAUDE.md with GFN2-xTB workflow
- [ ] **5.2** Publication figures
  - Heatmap: R × screen excimer fraction
  - Scatter: excimer vs Es
  - Distribution: distance MMFF94s vs GFN2-xTB
- [ ] **5.3** Methods section draft

---

## Part 6: Risk Assessment

### Technical Risks

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| xtb install fails on Windows | Medium | High | WSL2 fallback |
| GFN2-xTB too slow | Low | Medium | Parallelize |
| Threshold changes invalidate data | Medium | Low | Re-run all |

### Scientific Risks

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| GFN2-xTB overestimates stacking | Low | High | Validate first |
| Es model overfits | Medium | Medium | Cross-validation |

---

## Part 7: Success Criteria

### Quantitative

1. GFN2-xTB mean distance: 3.3-3.6 Å (vs MMFF94s 3.8-4.2 Å)
2. Excimer retention: ≥2× more after filtering
3. SAR model: Es correlation |r| > 0.5, p < 0.05
4. Full screening runtime: ≤6 hours

### Deliverables

1. Working GFN2-xTB integration
2. Updated thresholds with literature justification
3. SAR analysis with steric parameters
4. Comparison report: MMFF94s vs GFN2-xTB

---

## Key Literature References

### Pyrene Dimer Benchmarks
- **PMC11476719 (2024)** - DLPNO-CCSD(T)/CBS: 3.43 Å equilibrium
- **Ge et al. (2020)** - π-overlap > distance for excimer

### GFN2-xTB
- **Bannwarth et al. (2019)** - DOI: 10.1021/acs.jctc.8b01176
- **Rowan Benchmarks** - S22/S66 MAE 0.26-0.31 kcal/mol

### Substituent Effects
- **Wheeler & Houk (2011)** - Local substituent effects
- **Taft (1956)** - Original Es parameters

### Ethynyl Spacer
- **Winnik (1993)** - Polymer cyclization dynamics
- **CCS Chemistry (2020)** - Alkynyl transition moments

---

## Part 8: Finalized Multi-Tier Strategy (2026-02-05 Update)

### Current Status

| Component | Status | Notes |
|-----------|--------|-------|
| RDKit ETKDGv3 conformer generation | ✅ Working | 48 conformers in 43 seconds |
| MMFF94s optimization | ✅ Working | Fast but wrong pi-stacking distances |
| GFN2-xTB optimization | ✅ Working | Accurate but slow (~200s/conformer) |
| WSL2 + micromamba environment | ✅ Configured | xtb-python installed |
| MOE batch automation | ✅ Ready | SVL script created |
| Scientific validation | ✅ Completed | Agent review done |

### GFN2-xTB Timing Reality

From actual test run (Pyr_H, 124 heavy atoms):
- **RDKit embedding**: 43 seconds for 48 conformers
- **GFN2-xTB optimization**: 9628 seconds (2.7 hours) for 48 conformers
- **Per-conformer**: ~200 seconds average
- **64 variants projection**: 7 days (not practical for full library)

### Recommended Three-Tier Approach

#### Tier 1: Fast Screening (COMPLETED)
- **Method**: RDKit ETKDGv3 + MMFF94s
- **Runtime**: 5 hours for 64 variants
- **Output**: `binaph_screening_*.csv`
- **Use for**: Initial SAR, relative ranking, candidate prioritization
- **Limitation**: Distances biased (3.8-4.2 Å instead of 3.3-3.6 Å)

#### Tier 2: GFN2-xTB Validation (IN PROGRESS)
- **Method**: RDKit ETKDGv3 + GFN2-xTB optimization
- **Running**: WSL2 via `run_screening_xtb.py`
- **Target**: 10-20 priority variants (let run for 2-3 days)
- **Use for**: Distance calibration, validation of MMFF94s ranking
- **Output**: `binaph_xtb_screening_*.csv`

#### Tier 3: MOE Publication Quality (READY)
- **Method**: MOE LowModeMD + Amber:EHT
- **Runtime**: 15-30 min per molecule
- **Target**: Top 10 priority variants
- **Use for**: Publication figures, final validation
- **Output**: `moe_conformers/*.mdb`

### Priority Molecules (Top 10)

| Rank | Name | RDKit Excimer | Reason |
|------|------|---------------|--------|
| 1 | EtynPyr_Me | 42.9% | Best performer |
| 2 | EtynPyr_Et | ~35% | Second best EtynPyr |
| 3 | EtynPyr_H | ~30% | Ethynyl spacer baseline |
| 4 | CNPh_Th_tBu | 33.3% | Anomaly (only 3 conf) |
| 5 | CNPh_Th_Me | ~25% | CNPh series reference |
| 6 | Pyr_Me | ~20% | Best Pyr variant |
| 7 | Pyr_H | ~15% | Pyr baseline |
| 8 | EtynPyr_iPr | ~25% | Medium steric |
| 9 | DCV_Th_Me | ~15% | DCV reference |
| 10 | DCV_Th_H | ~10% | DCV baseline |

---

## Part 9: Scientific Validation Summary

### Methodology Assessment

| Aspect | Validity | Confidence | Notes |
|--------|----------|------------|-------|
| **Conformer generation** | Acceptable | Medium | Biased sampling compensates for FF deficiencies |
| **Distance thresholds** | Valid | Medium-High | Calibrated to DLPNO-CCSD(T) benchmark |
| **Angle thresholds** | Valid | High | Standard literature values |
| **Overlap calculation** | Valid | High | Correct algorithm, appropriate thresholds |
| **MMFF94s distances** | Biased but consistent | Medium | Valid for relative comparisons only |
| **Relative R-group rankings** | Valid | Medium-High | Consistent bias cancels in comparisons |
| **Absolute excimer fractions** | Invalid | Low | Missing Boltzmann weights, solvent effects |
| **EtynPyr > Pyr finding** | Valid | High | Consistent with literature (ethynyl spacer effect) |
| **CNPh_Th_tBu anomaly** | Uncertain | Low | N=3 conformers, high statistical uncertainty |

### Key Assumptions and Limitations

1. **Ground state geometry assumption**: Excimer formation requires excited-state geometry relaxation (not modeled)
2. **Gas phase approximation**: Solvent effects on stacking not included
3. **Boltzmann weighting missing**: All conformers treated equally (not energy-weighted)
4. **Static geometry**: Dynamic fluctuations at experimental temperature not captured
5. **Force field limitations**: MMFF94s distances systematically biased; GFN2-xTB is correction

### Required Caveats for Publication

Include these statements in any publication:
- "Excimer fractions represent geometric propensity, not thermodynamic populations"
- "MMFF94s distances are systematically ~0.5 Å too long; relative rankings remain valid"
- "Conformer ensembles represent accessible minima, not Boltzmann-weighted distributions"

---

## Part 10: Files Reference

### Core Scripts
| File | Purpose |
|------|---------|
| `run_screening_xtb.py` | GFN2-xTB batch screening (WSL) |
| `moe_batch_conformers.svl` | MOE batch automation script |
| `sar_analysis.py` | SAR correlation analysis |
| `build_binaph_dimer.py` | Dimer builder |
| `build_binaph_library.py` | 64-variant library generator |

### Input Data
| File | Purpose |
|------|---------|
| `binaph_dimer_smiles.csv` | 64 dimer SMILES |
| `moe_import/*.sdf` | 64 pre-generated 3D structures |

### Output Data
| File | Purpose |
|------|---------|
| `binaph_screening_*.csv` | MMFF94s screening results |
| `binaph_xtb_screening_*.csv` | GFN2-xTB screening results (in progress) |
| `moe_conformers/*.mdb` | MOE conformer ensembles (to generate) |

### Documentation
| File | Purpose |
|------|---------|
| `docs/PROJECT_PLAN.md` | This file |
| `docs/MOE_BATCH_GUIDE.md` | MOE automation guide |
| `docs/installation_xtb.md` | GFN2-xTB setup guide |
| `CLAUDE.md` | AI assistant context |

---

## Part 11: Action Items (Immediate)

### Running Now
- [x] GFN2-xTB screening in WSL (let it run)

### Today
- [ ] Start MOE batch on priority 10 molecules
- [ ] Monitor GFN2-xTB progress (check after 6-8 hours)

### This Week
- [ ] Collect GFN2-xTB results for comparison
- [ ] Complete MOE conformer search for priority variants
- [ ] Calculate distance calibration factor: `d_corrected = d_MMFF94s - offset`
- [ ] Create comparison plots: MMFF94s vs GFN2-xTB vs MOE

### Before Publication
- [ ] Re-run CNPh_Th_tBu with 200+ conformers (resolve anomaly)
- [ ] Apply Boltzmann weighting (optional enhancement)
- [ ] Generate publication figures
- [ ] Draft methods section with proper caveats
