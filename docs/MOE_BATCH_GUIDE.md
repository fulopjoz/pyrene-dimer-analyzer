# MOE Batch Conformer Search Guide

## Overview

This guide describes how to run automated conformer searches for all 64 binaphthalene dimers using MOE (Molecular Operating Environment).

**Input**: 64 pre-generated 3D SDF files in `moe_import/` directory
**Output**: 64 MDB files with conformer ensembles in `moe_conformers/` directory
**Expected runtime**: 15-30 min per molecule = 16-32 hours total

---

## Quick Start

### Option A: Automated Batch (Recommended)

1. **Open MOE**
2. **Set working directory** to the project root:
   - `MOE | File | Directory...` → navigate to `pyrene-dimer-analyzer/`
3. **Load the batch script**:
   - `MOE | SVL | Commands`
   - Click `Load...` → select `moe_batch_conformers.svl`
4. **Run the batch**:
   ```svl
   batch_conformer_search []
   ```
   Or for priority molecules only:
   ```svl
   batch_priority []
   ```

### Option B: Manual Processing (One at a Time)

1. **Open MOE**
2. **Import SDF**:
   - `File | Open...` → select file from `moe_import/` (e.g., `EtynPyr_Me_3D.sdf`)
3. **Run QuickPrep**:
   - `Compute | QuickPrep` → Use defaults → OK
4. **Run Conformer Search**:
   - `Compute | Conformations | Search...`
   - Use parameters from table below
   - Set output file: `moe_conformers/{name}_conformers.mdb`
5. **Repeat** for each molecule

---

## Optimized Parameters for Large Dimers

These parameters are specifically tuned for molecules with 100-150 atoms:

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| **Method** | LowModeMD | Best for large flexible molecules |
| **Temperature** | **600 K** | Higher than default (300K) to cross barriers |
| **Iteration Limit** | 10000 | Sufficient sampling for complex potential surface |
| **Rejection Limit** | **500** | Higher than default (100) to avoid premature termination |
| **RMSD Limit** | **0.25 Å** | Tighter than default (0.5) to capture subtle variations |
| **Energy Window** | 10 kcal/mol | Standard value |
| **Conformation Limit** | 200 | Maximum conformers to keep |
| **Force Field** | Amber:EHT | Better dispersion than MMFF94 |

**Critical**: Temperature=600K and Rejection=500 are essential. Default values (300K, 100) will miss many conformers.

---

## Priority Molecules for MOE

Based on RDKit screening results, these 10 molecules should be processed first:

| # | Name | Why Priority |
|---|------|--------------|
| 1 | EtynPyr_Me | Top performer (42.9% excimer) |
| 2 | EtynPyr_Et | Second best EtynPyr |
| 3 | EtynPyr_H | Baseline EtynPyr (ethynyl effect) |
| 4 | EtynPyr_iPr | Medium steric effect |
| 5 | CNPh_Th_tBu | Anomalous (33% with bulky tBu) |
| 6 | CNPh_Th_Me | CNPh series reference |
| 7 | Pyr_Me | Best Pyr variant |
| 8 | Pyr_H | Baseline Pyr |
| 9 | DCV_Th_Me | Best DCV variant |
| 10 | DCV_Th_H | Baseline DCV |

To process only these:
```svl
batch_priority []
```

---

## Directory Structure

```
pyrene-dimer-analyzer/
├── moe_import/                    # INPUT: 64 pre-generated 3D SDFs
│   ├── Pyr_H_3D.sdf
│   ├── Pyr_Me_3D.sdf
│   ├── ...
│   └── CNPh_Th_Ph_3D.sdf
├── moe_conformers/                # OUTPUT: Conformer ensembles (created by script)
│   ├── Pyr_H_conformers.mdb
│   ├── ...
├── moe_batch_conformers.svl       # Batch automation script
└── docs/MOE_BATCH_GUIDE.md        # This file
```

---

## Workflow After MOE

After MOE conformer search completes:

### 1. Export Conformers to SDF

For each MDB file, export conformers:
- Open MDB in MOE Database Viewer
- `DBV | File | Export | SD File...`
- Save as `moe_conformers/{name}_conformers.sdf`

Or use the provided SVL export script:
```svl
run ['export_moe_conformers.svl', []]
```

### 2. Analyze with pyrene-analyzer

```bash
# Analyze MOE conformers
pyrene-analyze analyze moe_conformers/EtynPyr_Me_conformers.sdf \
    --system binaphthalene -o moe_analysis_EtynPyr_Me.csv

# Batch analysis
pyrene-analyze analyze moe_conformers/*.sdf \
    --system binaphthalene -o moe_analysis_all.csv
```

### 3. Compare with RDKit Results

```python
import pandas as pd

# Load both datasets
rdkit_df = pd.read_csv("binaph_screening_summary.csv")
moe_df = pd.read_csv("moe_analysis_summary.csv")

# Merge and compare
comparison = rdkit_df.merge(moe_df, on="name", suffixes=("_rdkit", "_moe"))
comparison["excimer_diff"] = comparison["excimer_fraction_moe"] - comparison["excimer_fraction_rdkit"]

# Correlation
print(f"Excimer fraction correlation: {comparison['excimer_fraction_rdkit'].corr(comparison['excimer_fraction_moe']):.3f}")
```

---

## Troubleshooting

### "Conformer search produces no conformers"

1. **Check QuickPrep ran**: Atoms need partial charges
2. **Increase Temperature**: Try 700K or 800K
3. **Check geometry**: Atoms may be overlapping

### "Search terminates quickly"

1. **Increase Rejection Limit**: Set to 500-1000
2. **Increase Iteration Limit**: Set to 15000-20000

### "Memory issues"

1. **Reduce Conformation Limit**: Set to 100 instead of 200
2. **Process fewer molecules at once**: Use `batch_conformer_subset`

### "Script errors"

1. **Check working directory**: Must be project root
2. **Check file paths**: Ensure `moe_import/` exists with SDF files
3. **Ensure output dir exists**: Create `moe_conformers/` manually if needed

---

## SVL Script Commands Reference

```svl
// Process all 64 molecules (full batch)
batch_conformer_search []

// Process only priority molecules (10 total)
batch_priority []

// Process specific molecules
batch_conformer_subset ['EtynPyr_Me', 'EtynPyr_Et', 'Pyr_H']

// Show help
batch_help []
```

---

## Alternative: MOE Batch Mode (Command Line)

If you have MOE batch license:

```bash
# Windows
moebatch -load moe_batch_conformers.svl -exec "batch_conformer_search []"

# Linux
$MOE/bin/moebatch -load moe_batch_conformers.svl -exec "batch_conformer_search []"
```

This allows running without GUI (e.g., overnight on server).

---

## Expected Results

For each molecule, expect:
- **50-200 conformers** (depends on flexibility)
- **Energy range**: 0-10 kcal/mol above minimum
- **Runtime**: 15-30 minutes per molecule

For the full library:
- **Total conformers**: ~3,000-10,000
- **Total runtime**: 16-32 hours
- **Disk space**: ~500 MB - 2 GB

---

## References

- MOE Conformational Search: https://www.chemcomp.com/manuals/Applications/Confsearch.htm
- LowModeMD Method: Labute, P. J. Chem. Inf. Model. 2010, 50, 792-800
- Amber:EHT Force Field: Molecular Operating Environment (MOE) Documentation
