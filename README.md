# Aromatic Dimer Conformational Analysis Toolkit

[![Build Status](https://github.com/research-team/pyrene-dimer-analyzer/actions/workflows/tests.yml/badge.svg)](https://github.com/research-team/pyrene-dimer-analyzer/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/research-team/pyrene-dimer-analyzer/branch/main/graph/badge.svg)](https://codecov.io/gh/research-team/pyrene-dimer-analyzer)
[![PyPI version](https://badge.fury.io/py/pyrene-dimer-analyzer.svg)](https://badge.fury.io/py/pyrene-dimer-analyzer)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python versions](https://img.shields.io/pypi/pyversions/pyrene-dimer-analyzer.svg)](https://pypi.org/project/pyrene-dimer-analyzer)

**Automated RDKit-based analysis of aromatic dimer conformational properties**

Supports pyrene, perylene, anthracene, naphthalene, phenanthrene, and custom SMARTS patterns.

---

## Overview

`pyrene-dimer-analyzer` is a production-ready Python package for analyzing the geometric properties of aromatic dimer conformers. It is designed for computational chemistry workflows, particularly for studying how substituent groups (R-groups) affect inter-chromophore geometry and control photochemical properties like excimer vs. monomer fluorescence.

### Key Features

-   **Multi-System Support**: Built-in support for pyrene, perylene, anthracene, naphthalene, and phenanthrene with system-specific classification thresholds backed by literature data. Custom SMARTS patterns are also supported.
-   **Automatic Aromatic Detection**: Uses SMARTS patterns with overlap-aware clustering to reliably identify two distinct aromatic systems, with connectivity-based clustering as a fallback.
-   **Comprehensive Geometric Analysis**: Calculates critical parameters:
    -   **theta (Plane-plane angle)**: Angle between aromatic ring planes (0-90 deg)
    -   **d (Inter-plane distance)**: Perpendicular distance between planes (Angstrom)
    -   **pi-overlap**: Percentage of aromatic ring overlap (0-100%), calculated accurately with Shapely.
    -   **phi_L, phi_R (Bridge dihedrals)**: Conformational flexibility of linking bridges (degrees).
    -   **Energy**: Conformer relative energies (kcal/mol).
-   **High-Angle Geometry Warnings**: Flags conformers where theta > 60 deg in the `geometry_warnings` column (no console warning spam).
-   **Batch Processing**: Analyze multiple molecules and conformer files in parallel.
-   **Visualization**: Generate publication-quality plots with system-specific threshold overlays.
-   **QSAR-Ready Output**: Produces clean, structured data ready for QSAR/SAR modeling.
-   **Command-Line Interface**: A powerful CLI for easy integration into scripts and workflows.

### Supported Aromatic Systems

| System | Strong d (Angstrom) | Strong theta | Strong overlap | Key Reference |
| --- | --- | --- | --- | --- |
| Pyrene | 3.3-3.7 | < 20 deg | > 50% | Birks 1968, Ge 2020, Basuroy 2021 |
| Perylene | 3.4-3.8 | < 20 deg | > 50% | Crystal stacking data |
| Anthracene | 3.4-3.93 | < 20 deg | > 50% | High-pressure excimer literature |
| Naphthalene | 3.3-3.7 | < 20 deg | > 50% | Approximated from pyrene |
| Phenanthrene | 3.3-3.7 | < 20 deg | > 50% | Approximated from pyrene |
| Custom | User-defined | User-defined | User-defined | User-provided |

### Scientific Context

The study of aromatic dimers is crucial for understanding and designing fluorescent materials. The relative orientation of the two aromatic moieties determines whether the fluorescence emission is from the monomer (violet, ~375-400 nm) or the excimer (blue/green, ~460-500 nm). This package provides the tools to quantify the key geometric parameters that govern this behavior, based on established principles from the literature.

For pyrene (the best-characterized system), strong excimer formation requires theta < 20 deg, d = 3.3-3.7 Angstrom, and pi-overlap > 50% (Ge et al. 2020: 40-80% range; Basuroy et al. 2021: 42% overlap confirmed as excimer geometry).

---

## Installation

### Using pip

The recommended way to install `pyrene-dimer-analyzer` is from PyPI:

```bash
pip install pyrene-dimer-analyzer
```

### From Source

To install the latest development version from source:

```bash
git clone https://github.com/research-team/pyrene-dimer-analyzer
cd pyrene-dimer-analyzer
pip install -e .
```

---

## Quick Start

### As a Python Library

```python
from pyrene_analyzer import AromaticDimerAnalyzer

# Analyze pyrene dimers (default)
analyzer = AromaticDimerAnalyzer()
results_df = analyzer.analyze_file("conformers.sdf")
results_df.to_csv("analysis.csv")

# Analyze perylene dimers
analyzer = AromaticDimerAnalyzer(aromatic_system="perylene")
results_df = analyzer.analyze_file("perylene_conformers.sdf")

# Use custom SMARTS pattern
analyzer = AromaticDimerAnalyzer(custom_smarts="c1ccc2cc3ccccc3cc2c1")
results_df = analyzer.analyze_file("anthracene_dimers.sdf")

# Backward-compatible alias still works
from pyrene_analyzer import PyreneDimerAnalyzer
analyzer = PyreneDimerAnalyzer()  # defaults to pyrene
```

### From the Command Line

```bash
# Analyze pyrene dimers (default)
pyrene-analyze analyze conformers.sdf -o results.csv

# Analyze with a specific aromatic system
pyrene-analyze analyze perylene_dimers.sdf -o results.csv --system perylene

# Analyze multiple files with plots and classification
pyrene-analyze analyze Et.sdf iPr.sdf cHex.sdf -o all_results.csv --plot --classify

# Screen substituent variants on a template (R-group SMARTS must match a single
# heavy-atom attachment point; hydrogens are ignored during attachment checks)
pyrene-analyze screen dimer.sdf -r "[CH2][CH3]" -o screening_results.csv --plot

# List all supported systems and their thresholds
pyrene-analyze info

# Get help and see all options
pyrene-analyze analyze --help
```

---

## Documentation

Full documentation, including tutorials, API reference, and scientific background, is available at [pyrene-dimer-analyzer.readthedocs.io](https://pyrene-dimer-analyzer.readthedocs.io/).

---

## Citation

If you use `pyrene-dimer-analyzer` in your research, please cite it as follows:

> Research Team, UCT Prague. (2026). *pyrene-dimer-analyzer: Automated geometric analysis of aromatic dimer conformers*. GitHub. [https://github.com/research-team/pyrene-dimer-analyzer](https://github.com/research-team/pyrene-dimer-analyzer)

### Related Publications

This work is based on principles from the following key publications:

-   Birks, J.B. (1970). *Photophysics of Aromatic Molecules*.
-   Stevens, B. (1968). *Proc. Royal Soc. A*, 305, 55-70.
-   Ge, Y. et al. (2020). *J. Mater. Chem. C*, 8, 10223-10232.
-   Basuroy, K. et al. (2021). *J. Chem. Phys.*, 155, 234304.
-   Mazzeo, P. et al. (2024). *ChemRxiv* (twisted excimers).

---

## Contributing

Contributions are welcome! Please see the [Contributing Guidelines](https://github.com/research-team/pyrene-dimer-analyzer/blob/main/CONTRIBUTING.md) for more information. This project follows a [Code of Conduct](https://github.com/research-team/pyrene-dimer-analyzer/blob/main/CODE_OF_CONDUCT.md).

---

## License

This project is licensed under the MIT License. See the [LICENSE](https://github.com/research-team/pyrene-dimer-analyzer/blob/main/LICENSE) file for details.

---

## Acknowledgments

-   This project was developed by the Research Team at UCT Prague.
-   We thank the RDKit community for providing the essential cheminformatics toolkit.
