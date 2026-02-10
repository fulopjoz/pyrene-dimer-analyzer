# Welcome to Aromatic Dimer Analyzer

**`pyrene-dimer-analyzer`** is a powerful and easy-to-use Python package for the automated geometric analysis of aromatic dimer conformers. It supports **multiple aromatic systems** including pyrene, perylene, anthracene, naphthalene, phenanthrene, and custom SMARTS-defined systems. The package is designed as a key component in computational chemistry and materials science workflows, enabling researchers to study the relationship between molecular structure and photochemical properties.

This documentation will guide you through the installation, usage, and scientific principles behind the package.

## Key Features

-   **Multi-System Support**: Analyze pyrene, perylene, anthracene, naphthalene, phenanthrene, or any aromatic system via custom SMARTS patterns.
-   **Automated Analysis**: Automatically identify aromatic moieties and calculate key geometric parameters.
-   **Comprehensive Metrics**: Compute plane-plane angle (θ), inter-plane distance (d), π-overlap, slip-stack displacement, tilt angle, and more.
-   **System-Specific Classification**: Classify conformers as strong excimer, weak excimer, or monomer using literature-backed, system-specific thresholds.
-   **Geometry Warnings**: Flag conformers where θ > 60° (inter-plane distance unreliable at high angles).
-   **High-Quality Visualization**: Generate publication-ready plots with system-specific threshold overlays.
-   **Batch Processing**: Analyze thousands of conformers from multiple files in parallel.
-   **Flexible I/O**: Supports common molecular file formats (SDF, MOL2, PDB) and exports results to CSV, JSON, and Excel.
-   **Command-Line Interface**: A full-featured CLI with `--system` and `--custom-smarts` options for easy integration into automated workflows.

## Supported Aromatic Systems

| System | Rings | Optimal d (A) | Strong Overlap |
|---|---|---|---|
| Pyrene | 4 | 3.3-3.7 | > 50% |
| Perylene | 5 | 3.4-3.8 | > 50% |
| Anthracene | 3 | 3.4-3.93 | > 50% |
| Naphthalene | 2 | 3.3-3.7 | > 50% |
| Phenanthrene | 3 | 3.3-3.7 | > 50% |
| Custom | - | User-defined | User-defined |

## Getting Started

-   **[Installation](installation.md)**: Learn how to install the package.
-   **[Quick Start](quickstart.md)**: A quick tutorial to get you started.

## User Guide

-   **[API Reference](api_reference.md)**: Detailed information about the functions and classes.
-   **[Scientific Background](scientific_background.md)**: Understand the theory behind the analysis.
-   **[Examples](examples/basic_analysis.md)**: Explore real-world usage examples.

## For Developers

-   **Contributing**: We welcome contributions! Please see our contributing guidelines on GitHub.
-   **Source Code**: The source code is available on [GitHub](https://github.com/research-team/pyrene-dimer-analyzer).
