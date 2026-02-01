"""
Pyrene Dimer Analyzer
=====================

A Python package for automated geometric analysis of pyrene dimer conformers.

This package provides tools for analyzing pyrene dimer molecular geometries,
including plane-plane angles, inter-plane distances, π-π overlap calculations,
and bridge dihedral angles. It is designed for use in computational chemistry
workflows, particularly for studying excimer formation in pyrene-based systems.

Key Features:
    - Automatic detection of pyrene ring systems using SMARTS patterns
    - SVD-based plane fitting for accurate angle calculations
    - Shapely-based polygon intersection for π-overlap estimation
    - Support for multiple file formats (SDF, MOL2, PDB)
    - Batch processing with parallel execution support
    - Publication-quality visualization
    - Command-line interface

Example:
    >>> from pyrene_analyzer import PyreneDimerAnalyzer
    >>> analyzer = PyreneDimerAnalyzer()
    >>> results = analyzer.analyze_file('conformers.sdf')
    >>> results.to_csv('analysis.csv')

Scientific Background:
    Pyrene excimer formation requires specific geometric conditions:
    - Plane-plane angle (θ) < 20° for optimal excimer formation
    - Inter-plane distance (d) of 3.3-3.7 Å for π-stacking
    - π-overlap > 50% for efficient electronic coupling

References:
    - Birks, J.B. (1970). Photophysics of Aromatic Molecules
    - Stevens, B. (1968). Proc. Royal Soc. A, 305, 55-70
    - Poisson, L. et al. (2017). Phys. Chem. Chem. Phys., 19, 23492-23506
    - Ge, Y. et al. (2020). J. Mater. Chem. C, 8, 10223-10232

Author:
    Research Team, UCT Prague

License:
    MIT License
"""

__version__ = "1.0.0"
__author__ = "Research Team, UCT Prague"
__email__ = "research@vscht.cz"
__license__ = "MIT"

from pyrene_analyzer.core import PyreneDimerAnalyzer
from pyrene_analyzer.geometry import (
    fit_plane_svd,
    calculate_plane_angle,
    calculate_interplane_distance,
    calculate_pi_overlap,
    calculate_centroid_distance,
    calculate_slip_stack_displacement,
)
from pyrene_analyzer.io import (
    load_from_sdf,
    load_from_mol2,
    load_from_pdb,
    export_to_csv,
    export_to_json,
    export_to_excel,
)
from pyrene_analyzer.visualization import (
    plot_angle_vs_energy,
    plot_distance_vs_overlap,
    plot_conformer_distribution,
    plot_variant_comparison,
)

__all__ = [
    # Version info
    "__version__",
    "__author__",
    "__email__",
    "__license__",
    # Core
    "PyreneDimerAnalyzer",
    # Geometry
    "fit_plane_svd",
    "calculate_plane_angle",
    "calculate_interplane_distance",
    "calculate_pi_overlap",
    "calculate_centroid_distance",
    "calculate_slip_stack_displacement",
    # I/O
    "load_from_sdf",
    "load_from_mol2",
    "load_from_pdb",
    "export_to_csv",
    "export_to_json",
    "export_to_excel",
    # Visualization
    "plot_angle_vs_energy",
    "plot_distance_vs_overlap",
    "plot_conformer_distribution",
    "plot_variant_comparison",
]
