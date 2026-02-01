# Pyrene Dimer Conformational Analysis Toolkit

[![Build Status](https://github.com/research-team/pyrene-dimer-analyzer/actions/workflows/tests.yml/badge.svg)](https://github.com/research-team/pyrene-dimer-analyzer/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/research-team/pyrene-dimer-analyzer/branch/main/graph/badge.svg)](https://codecov.io/gh/research-team/pyrene-dimer-analyzer)
[![PyPI version](https://badge.fury.io/py/pyrene-dimer-analyzer.svg)](https://badge.fury.io/py/pyrene-dimer-analyzer)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python versions](https://img.shields.io/pypi/pyversions/pyrene-dimer-analyzer.svg)](https://pypi.org/project/pyrene-dimer-analyzer)

**Automated RDKit-based analysis of pyrene dimer conformational properties**

---

## 🎯 Overview

`pyrene-dimer-analyzer` is a production-ready Python package for analyzing the geometric properties of pyrene dimer conformers. It is designed for computational chemistry workflows, particularly for studying how substituent groups (R-groups) affect inter-pyrene geometry and control photochemical properties like excimer vs. monomer fluorescence.

### Key Features

-   **Automatic Pyrene Detection**: Uses SMARTS patterns to reliably identify the two pyrene ring systems.
-   **Comprehensive Geometric Analysis**: Calculates critical parameters:
    -   **θ (Plane-plane angle)**: Angle between pyrene ring planes (0-90°)
    -   **d (Inter-plane distance)**: Perpendicular distance between planes (Å)
    -   **π-overlap**: Percentage of aromatic ring overlap (0-100%), calculated accurately with Shapely.
    -   **φL, φR (Bridge dihedrals)**: Conformational flexibility of linking bridges (degrees).
    -   **Energy**: Conformer relative energies (kcal/mol).
-   **Batch Processing**: Analyze multiple molecules and conformer files in parallel.
-   **Visualization**: Generate publication-quality plots for in-depth analysis.
-   **QSAR-Ready Output**: Produces clean, structured data ready for QSAR/SAR modeling.
-   **Command-Line Interface**: A powerful CLI for easy integration into scripts and workflows.

### Scientific Context

The study of pyrene dimers is crucial for understanding and designing fluorescent materials. The relative orientation of the two pyrene moieties determines whether the fluorescence emission is from the monomer (violet, ~375-400 nm) or the excimer (blue/green, ~460-500 nm). This package provides the tools to quantify the key geometric parameters that govern this behavior, based on established principles from the literature.

---

## 📦 Installation

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

## 🚀 Quick Start

### As a Python Library

```python
from pyrene_analyzer import PyreneDimerAnalyzer

# Initialize the analyzer
analyzer = PyreneDimerAnalyzer()

# Analyze a file with conformers
results_df = analyzer.analyze_file("conformers.sdf")

# Save results to CSV
results_df.to_csv("analysis.csv")

print(results_df.head())
```

### From the Command Line

```bash
# Analyze a single SDF file and save to CSV
pyrene-analyze analyze conformers.sdf -o results.csv

# Analyze multiple files, generate plots, and run in verbose mode
pyrene-analyze analyze Et.sdf iPr.sdf cHex.sdf -o all_results.csv --plot --verbose

# Get help and see all options
pyrene-analyze analyze --help
```

---

## 📖 Documentation

Full documentation, including tutorials, API reference, and scientific background, is available at [pyrene-dimer-analyzer.readthedocs.io](https://pyrene-dimer-analyzer.readthedocs.io/).

---

## 🎓 Citation

If you use `pyrene-dimer-analyzer` in your research, please cite it as follows:

> Research Team, UCT Prague. (2026). *pyrene-dimer-analyzer: Automated geometric analysis of pyrene dimer conformers*. GitHub. [https://github.com/research-team/pyrene-dimer-analyzer](https://github.com/research-team/pyrene-dimer-analyzer)

### Related Publications

This work is based on principles from the following key publications:

-   Birks, J.B. (1970). *Photophysics of Aromatic Molecules*.
-   Stevens, B. (1968). *Proc. Royal Soc. A*, 305, 55-70.
-   Poisson, L. et al. (2017). *Phys. Chem. Chem. Phys.*, 19, 23492-23506.
-   Ge, Y. et al. (2020). *J. Mater. Chem. C*, 8, 10223-10232.

---

## 🤝 Contributing

Contributions are welcome! Please see the [Contributing Guidelines](https://github.com/research-team/pyrene-dimer-analyzer/blob/main/CONTRIBUTING.md) for more information. This project follows a [Code of Conduct](https://github.com/research-team/pyrene-dimer-analyzer/blob/main/CODE_OF_CONDUCT.md).

---

## 📝 License

This project is licensed under the MIT License. See the [LICENSE](https://github.com/research-team/pyrene-dimer-analyzer/blob/main/LICENSE) file for details.

---

## 🙏 Acknowledgments

-   This project was developed by the Research Team at UCT Prague.
-   We thank the RDKit community for providing the essential cheminformatics toolkit.
