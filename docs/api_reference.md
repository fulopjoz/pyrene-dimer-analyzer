# API Reference

This section provides a detailed reference for the classes and functions in the `pyrene-dimer-analyzer` package. The documentation is organized by module.

## Aromatic Systems Module (`pyrene_analyzer.aromatic_systems`)

This module provides the registry of supported aromatic systems with their SMARTS patterns and classification thresholds.

### `ClassificationThresholds`

```python
@dataclass(frozen=True)
class ClassificationThresholds:
    strong_angle_max: float        # Max angle for strong excimer (deg)
    strong_distance_range: Tuple[float, float]  # (min, max) distance (A)
    strong_overlap_min: float      # Min overlap for strong excimer (%)
    weak_angle_max: float          # Max angle for weak excimer (deg)
    weak_distance_max: float       # Max distance for weak excimer (A)
    weak_overlap_min: float        # Min overlap for weak excimer (%)
    high_angle_warning: float = 60.0  # Angle above which distance is unreliable
```

Frozen dataclass holding all classification thresholds for an aromatic system. All values are backed by literature data.

### `AromaticSystem`

```python
@dataclass(frozen=True)
class AromaticSystem:
    name: str                      # System name (e.g., "pyrene")
    smarts: str                    # SMARTS pattern for detection
    min_ring_atoms: int            # Minimum atoms for fallback detection
    expected_ring_atoms: int       # Expected atoms in the aromatic system
    thresholds: ClassificationThresholds
    references: Tuple[str, ...]    # Literature references
    notes: str                     # Additional notes
```

Frozen dataclass representing a complete aromatic system definition.

### Functions

- `get_system(name)`: Get an `AromaticSystem` by name. Raises `ValueError` if not found.
- `list_systems()`: Return a list of all registered system names.
- `register_system(system)`: Register a new `AromaticSystem` in the global registry.

### `AROMATIC_SYSTEMS`

A `dict` mapping system names to `AromaticSystem` instances. Contains 5 built-in systems: `pyrene`, `perylene`, `anthracene`, `naphthalene`, `phenanthrene`.

## Core Module (`pyrene_analyzer.core`)

This module contains the main `AromaticDimerAnalyzer` class that orchestrates the analysis.

### `AromaticDimerAnalyzer`

```python
class AromaticDimerAnalyzer(
    verbose: bool = True,
    use_smarts: bool = True,
    use_shapely: bool = True,
    aromatic_system: str = "pyrene",
    custom_smarts: Optional[str] = None,
    custom_thresholds: Optional[ClassificationThresholds] = None
)
```

The main class for analyzing aromatic dimer geometric properties.

**Parameters:**

- `verbose` (bool): If `True`, print progress messages during analysis.
- `use_smarts` (bool): If `True`, use SMARTS patterns for aromatic detection. If `False`, use a connectivity-based clustering algorithm.
- `use_shapely` (bool): If `True`, use the Shapely library for accurate pi-overlap calculation. If `False`, a centroid-based approximation is used as a fallback.
- `aromatic_system` (str): Name of the aromatic system to use (default: `"pyrene"`). Must be a registered system name.
- `custom_smarts` (str, optional): Custom SMARTS pattern. When provided, overrides the system's SMARTS.
- `custom_thresholds` (ClassificationThresholds, optional): Custom thresholds. When provided with `custom_smarts`, overrides the default pyrene thresholds.

#### Methods

- `analyze_file(filename, show_progress=True)`: Analyze all molecules and conformers in a given file.
- `batch_analyze(filenames, n_jobs=1, show_progress=True)`: Analyze multiple files, with support for parallel processing.
- `analyze_molecule(mol, mol_name="molecule", show_progress=True)`: Analyze all conformers of a single RDKit molecule object.
- `analyze_conformer(mol, conf_id, pyrene1_atoms, pyrene2_atoms)`: Analyze a single conformer. Returns a dict including a `geometry_warnings` field.
- `identify_aromatic_rings(mol)`: Identify the two aromatic ring systems in a molecule using SMARTS or fallback.
- `identify_pyrene_rings(mol)`: Backward-compatible alias for `identify_aromatic_rings()`.
- `classify_conformer(angle, distance, overlap)`: Classify a conformer as `strong_excimer`, `weak_excimer`, or `monomer` using system-specific thresholds.
- `add_classification(df)`: Add a classification column to a results DataFrame.

### `PyreneDimerAnalyzer`

Backward-compatible alias for `AromaticDimerAnalyzer`. Defaults to pyrene system.

## Geometry Module (`pyrene_analyzer.geometry`)

This module provides functions for all geometric calculations.

### Functions

- `fit_plane_svd(coords)`: Fits a plane to a set of 3D coordinates using SVD. Returns `(centroid, normal)`.
- `calculate_plane_angle(coords1, coords2)`: Calculates the angle (0-90 deg) between two planes.
- `calculate_interplane_distance(coords1, coords2, plane_angle=None)`: Calculates the perpendicular distance (A) between two planes. When `plane_angle` is provided and > 60 deg, emits a `UserWarning` about unreliable metric.
- `calculate_pi_overlap(coords1, coords2, use_shapely=True)`: Calculates the pi-pi overlap percentage (0-100%) between two ring systems.
- `calculate_centroid_distance(coords1, coords2)`: Calculates the Euclidean distance (A) between the centroids of two coordinate sets.
- `calculate_slip_stack_displacement(coords1, coords2)`: Calculates the lateral (slip-stack) displacement (A) between two planes.
- `calculate_tilt_angle(coords1, coords2)`: Calculates the tilt angle (0-180 deg) between centroid-centroid vector and plane normal.

## I/O Module (`pyrene_analyzer.io`)

This module handles file input and output operations.

### Functions

- `load_molecules(filename, remove_hydrogens=False)`: Load molecules from a file, automatically detecting the format (SDF, MOL2, PDB).
- `export_to_csv(df, filename, float_format='%.3f')`: Export a DataFrame to a CSV file.
- `export_to_json(df, filename, orient='records', indent=2)`: Export a DataFrame to a JSON file.
- `export_to_excel(df, filename, sheet_name='Analysis')`: Export a DataFrame to an Excel file.
- `export_results(df, filename, formats=None)`: Export a DataFrame to multiple formats at once.

## Visualization Module (`pyrene_analyzer.visualization`)

This module provides functions for creating publication-quality plots. All threshold-dependent functions accept an optional `thresholds` parameter of type `ClassificationThresholds` to support different aromatic systems (defaults to pyrene).

### Functions

- `plot_angle_vs_energy(df, output=None, thresholds=None, ...)`: Creates a scatter plot of plane angle vs. relative energy. Excimer region shading uses system-specific thresholds.
- `plot_distance_vs_overlap(df, output=None, thresholds=None, ...)`: Creates a scatter plot of inter-plane distance vs. pi-overlap. Optimal region rectangle uses system-specific thresholds.
- `plot_conformer_distribution(df, output=None, ...)`: Creates histograms showing the distribution of geometric parameters.
- `plot_variant_comparison(df, output=None, thresholds=None, ...)`: Creates box plots to compare geometric parameters across different R-group variants. Reference lines use system-specific thresholds.
- `plot_energy_landscape(df, output=None, ...)`: Creates a 2D energy landscape plot.
- `plot_correlation_matrix(df, output=None, ...)`: Creates a correlation matrix heatmap of geometric parameters.
- `create_summary_figure(df, output=None, thresholds=None, ...)`: Creates a comprehensive summary figure with multiple panels.
