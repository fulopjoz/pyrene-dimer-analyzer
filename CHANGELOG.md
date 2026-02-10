# Changelog

All notable changes to the `pyrene-dimer-analyzer` package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-02-01

### Added

- **Multi-system aromatic support**: The analyzer now supports 5 built-in aromatic systems
  (pyrene, perylene, anthracene, naphthalene, phenanthrene) plus custom SMARTS patterns.
- **`AromaticDimerAnalyzer` class**: New primary class name accepting `aromatic_system`,
  `custom_smarts`, and `custom_thresholds` parameters.
- **Aromatic systems registry** (`pyrene_analyzer/aromatic_systems.py`):
  `ClassificationThresholds` and `AromaticSystem` frozen dataclasses,
  `register_system()`, `get_system()`, `list_systems()` functions.
- **High-angle geometry warnings**: `calculate_interplane_distance()` now accepts an
  optional `plane_angle` parameter and emits a `UserWarning` when angle > 60 deg,
  because the inter-plane distance metric measures edge-to-face distance rather than
  pi-stacking distance at high angles.
- **`geometry_warnings` field** in conformer output dicts, collecting all geometry
  anomaly warnings per conformer.
- **CLI `--system` / `-s` option**: Select aromatic system for analysis
  (e.g., `pyrene-analyze analyze file.sdf -o out.csv --system perylene`).
- **CLI `--custom-smarts` option**: Provide an arbitrary SMARTS pattern for detection.
- **Parameterized visualization thresholds**: `plot_angle_vs_energy()`,
  `plot_distance_vs_overlap()`, `plot_variant_comparison()`, and
  `create_summary_figure()` accept an optional `thresholds` parameter.
- **`info` command** now dynamically lists all registered aromatic systems with
  their thresholds and literature references.
- **Virtual screening module** (`pyrene_analyzer/screening.py`):
  `SubstituentScreener` class for R-group enumeration, conformer generation,
  and batch geometric analysis. Includes `COMMON_SUBSTITUENTS` library (~25 entries).
- **Biased conformer generation**: `generate_conformers_biased()` uses distance
  constraints between representative aromatic atoms across 4 constraint levels
  (tight/medium/loose/unconstrained) to improve excimer geometry sampling.
- **SMILES-to-analysis pipeline**: `analyze_from_smiles()` provides a complete
  workflow from SMILES string to excimer classification without needing SDF files.
- **CLI `analyze-smiles` command**: Analyze a dimer directly from SMILES
  (e.g., `pyrene-analyze analyze-smiles "SMILES" -s naphthalene -n 100 -v`).
- **CLI `screen` command**: Virtual screening with R-group enumeration
  (e.g., `pyrene-analyze screen template.sdf "[CH2][CH3]" --subs Me Et iPr`).
- **Progress indicators**: `verbose` parameter on all screening functions with
  step-by-step timing, time estimates, and tqdm progress bars for analysis.
  CLI always shows coarse progress; `-v` enables detailed output.

### Fixed

- **Removed noisy sanitization warnings** from R-group replacement: valence
  violations from incompatible substituent/template combinations now return
  `None` silently instead of emitting `UserWarning`.
- **Removed noisy force field warnings**: Uses `MMFFHasAllMoleculeParams()`
  pre-check before MMFF94s optimization to avoid expected failures. Falls back
  to UFF silently; verbose-only messages when both fail.
- **Seaborn deprecation**: Removed explicit `orient="v"` from boxplot/stripplot
  calls in `plot_variant_comparison()` to avoid `PendingDeprecationWarning`.

### Changed

- **Strong overlap threshold corrected: 70% -> 50%**. Based on Ge et al. (2020,
  *J. Mater. Chem. C* 8, 10223: 40-80% range in crystals) and Basuroy et al.
  (2021, *J. Chem. Phys.* 155, 234304: 42% overlap confirmed as excimer geometry).
- System-specific classification thresholds loaded from registry instead of
  hardcoded constants.
- Minimum ring atom filter now uses system-specific value from registry
  (was hardcoded to >= 10).

### Backward Compatibility

- `PyreneDimerAnalyzer` is retained as an alias for `AromaticDimerAnalyzer`.
- `identify_pyrene_rings()` delegates to `identify_aromatic_rings()`.
- All existing API calls continue to work without changes.
- CLI `pyrene-analyze analyze file.sdf -o results.csv` works unchanged.
- Visualization functions work without `thresholds` parameter (defaults to pyrene).

### Scientific References

- Ge, Y. et al. (2020). *J. Mater. Chem. C*, 8, 10223-10232.
- Basuroy, K. et al. (2021). *J. Chem. Phys.*, 155, 234304.
- Mazzeo, P. et al. (2024). *ChemRxiv* (twisted excimers at 50-60 deg).
- Marazzi, M. et al. (2024). *J. Phys. Chem. Lett.* (TD-DFT d=3.24 A).

## [1.0.0] - 2026-01-15

### Added

- Initial release with pyrene-specific analysis.
- SVD-based plane fitting for accurate angle calculations.
- Shapely-based polygon intersection for pi-overlap estimation.
- Support for SDF, MOL2, and PDB file formats.
- Batch processing with parallel execution support.
- Publication-quality visualization (matplotlib/seaborn).
- Click-based CLI with `analyze`, `info`, and `preview` commands.
- Export to CSV, JSON, and Excel formats.
- 124 tests with comprehensive coverage.
