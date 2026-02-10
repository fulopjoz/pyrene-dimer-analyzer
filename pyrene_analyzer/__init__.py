"""
Aromatic Dimer Analyzer
=======================

A Python package for automated geometric analysis of aromatic dimer conformers.

This package provides tools for analyzing aromatic dimer molecular geometries,
including plane-plane angles, inter-plane distances, π-π overlap calculations,
and bridge dihedral angles. It supports multiple aromatic systems (pyrene,
perylene, anthracene, naphthalene, phenanthrene) with system-specific
classification thresholds backed by literature data.

Key Features:
    - Support for multiple aromatic systems with configurable SMARTS patterns
    - Automatic detection of aromatic ring systems using SMARTS + connectivity fallback
    - SVD-based plane fitting for accurate angle calculations
    - Shapely-based polygon intersection for π-overlap estimation
    - System-specific excimer classification thresholds
    - High-angle geometry warnings for unreliable distance metrics
    - Support for multiple file formats (SDF, MOL2, PDB)
    - Batch processing with parallel execution support
    - Publication-quality visualization
    - Command-line interface

Example:
    >>> from pyrene_analyzer import AromaticDimerAnalyzer
    >>> analyzer = AromaticDimerAnalyzer(aromatic_system="pyrene")
    >>> results = analyzer.analyze_file('conformers.sdf')
    >>> results.to_csv('analysis.csv')

    >>> # Backward-compatible alias
    >>> from pyrene_analyzer import PyreneDimerAnalyzer
    >>> analyzer = PyreneDimerAnalyzer()  # defaults to pyrene

Scientific Background:
    Aromatic excimer formation requires specific geometric conditions
    that vary by chromophore. For pyrene (the best-characterized system):
    - Plane-plane angle (θ) < 20° for optimal excimer formation
    - Inter-plane distance (d) of 3.3-3.7 Å for π-stacking
    - π-overlap > 50% for efficient electronic coupling

References:
    - Birks, J.B. (1970). Photophysics of Aromatic Molecules
    - Ge, Y. et al. (2020). J. Mater. Chem. C, 8, 10223-10232
    - Basuroy et al. (2021). J. Chem. Phys., 155, 234304
    - Mazzeo et al. (2024). ChemRxiv (twisted excimers)

Author:
    Research Team, UCT Prague

License:
    MIT License
"""

import importlib

__version__ = "1.1.0"
__author__ = "Research Team, UCT Prague"
__email__ = "research@vscht.cz"
__license__ = "MIT"

__all__ = [
    # Version info
    "__version__",
    "__author__",
    "__email__",
    "__license__",
    # Core
    "AromaticDimerAnalyzer",
    "PyreneDimerAnalyzer",
    # Aromatic systems registry
    "AromaticSystem",
    "ClassificationThresholds",
    "AROMATIC_SYSTEMS",
    "get_system",
    "list_systems",
    "register_system",
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
    # Screening
    "SubstituentScreener",
    "COMMON_SUBSTITUENTS",
    "generate_conformers_biased",
    "analyze_from_smiles",
    "find_representative_atoms",
    # MACE-OFF23 optimizer
    "has_mace_available",
    "optimize_with_mace",
    "optimize_conformer_ensemble_mace",
    # Ensemble
    "compute_ensemble_features",
    "compute_distributional_features",
    "compute_threshold_features",
    "compute_boltzmann_weighted_features",
    "KT_298K",
    # Calibration
    "DistanceCalibration",
    "make_offset_calibration",
    # Visualization
    "plot_angle_vs_energy",
    "plot_distance_vs_overlap",
    "plot_conformer_distribution",
    "plot_variant_comparison",
]

# Lazy import mapping: name -> (module_path, attribute_name)
_LAZY_IMPORTS = {
    # Aromatic systems registry
    "AROMATIC_SYSTEMS": ("pyrene_analyzer.aromatic_systems", "AROMATIC_SYSTEMS"),
    "AromaticSystem": ("pyrene_analyzer.aromatic_systems", "AromaticSystem"),
    "ClassificationThresholds": (
        "pyrene_analyzer.aromatic_systems",
        "ClassificationThresholds",
    ),
    "get_system": ("pyrene_analyzer.aromatic_systems", "get_system"),
    "list_systems": ("pyrene_analyzer.aromatic_systems", "list_systems"),
    "register_system": ("pyrene_analyzer.aromatic_systems", "register_system"),
    # Core
    "AromaticDimerAnalyzer": ("pyrene_analyzer.core", "AromaticDimerAnalyzer"),
    "PyreneDimerAnalyzer": ("pyrene_analyzer.core", "PyreneDimerAnalyzer"),
    # Geometry
    "fit_plane_svd": ("pyrene_analyzer.geometry", "fit_plane_svd"),
    "calculate_plane_angle": ("pyrene_analyzer.geometry", "calculate_plane_angle"),
    "calculate_interplane_distance": (
        "pyrene_analyzer.geometry",
        "calculate_interplane_distance",
    ),
    "calculate_pi_overlap": ("pyrene_analyzer.geometry", "calculate_pi_overlap"),
    "calculate_centroid_distance": (
        "pyrene_analyzer.geometry",
        "calculate_centroid_distance",
    ),
    "calculate_slip_stack_displacement": (
        "pyrene_analyzer.geometry",
        "calculate_slip_stack_displacement",
    ),
    "DistanceCalibration": ("pyrene_analyzer.geometry", "DistanceCalibration"),
    "make_offset_calibration": (
        "pyrene_analyzer.geometry",
        "make_offset_calibration",
    ),
    # I/O
    "load_from_sdf": ("pyrene_analyzer.io", "load_from_sdf"),
    "load_from_mol2": ("pyrene_analyzer.io", "load_from_mol2"),
    "load_from_pdb": ("pyrene_analyzer.io", "load_from_pdb"),
    "export_to_csv": ("pyrene_analyzer.io", "export_to_csv"),
    "export_to_json": ("pyrene_analyzer.io", "export_to_json"),
    "export_to_excel": ("pyrene_analyzer.io", "export_to_excel"),
    # Screening
    "SubstituentScreener": ("pyrene_analyzer.screening", "SubstituentScreener"),
    "COMMON_SUBSTITUENTS": ("pyrene_analyzer.screening", "COMMON_SUBSTITUENTS"),
    "generate_conformers_biased": (
        "pyrene_analyzer.screening",
        "generate_conformers_biased",
    ),
    "analyze_from_smiles": ("pyrene_analyzer.screening", "analyze_from_smiles"),
    "find_representative_atoms": (
        "pyrene_analyzer.screening",
        "find_representative_atoms",
    ),
    # MACE-OFF23 optimizer
    "has_mace_available": ("pyrene_analyzer.mace_optimizer", "has_mace_available"),
    "optimize_with_mace": ("pyrene_analyzer.mace_optimizer", "optimize_with_mace"),
    "optimize_conformer_ensemble_mace": (
        "pyrene_analyzer.mace_optimizer",
        "optimize_conformer_ensemble_mace",
    ),
    # Ensemble
    "compute_ensemble_features": (
        "pyrene_analyzer.ensemble",
        "compute_ensemble_features",
    ),
    "compute_distributional_features": (
        "pyrene_analyzer.ensemble",
        "compute_distributional_features",
    ),
    "compute_threshold_features": (
        "pyrene_analyzer.ensemble",
        "compute_threshold_features",
    ),
    "compute_boltzmann_weighted_features": (
        "pyrene_analyzer.ensemble",
        "compute_boltzmann_weighted_features",
    ),
    "KT_298K": ("pyrene_analyzer.ensemble", "KT_298K"),
    # Visualization
    "plot_angle_vs_energy": (
        "pyrene_analyzer.visualization",
        "plot_angle_vs_energy",
    ),
    "plot_distance_vs_overlap": (
        "pyrene_analyzer.visualization",
        "plot_distance_vs_overlap",
    ),
    "plot_conformer_distribution": (
        "pyrene_analyzer.visualization",
        "plot_conformer_distribution",
    ),
    "plot_variant_comparison": (
        "pyrene_analyzer.visualization",
        "plot_variant_comparison",
    ),
}


def __getattr__(name):
    if name in _LAZY_IMPORTS:
        module_path, attr_name = _LAZY_IMPORTS[name]
        module = importlib.import_module(module_path)
        val = getattr(module, attr_name)
        globals()[name] = val  # Cache for subsequent access
        return val
    raise AttributeError(f"module 'pyrene_analyzer' has no attribute {name!r}")
