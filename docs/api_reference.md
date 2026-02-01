# API Reference

This section provides a detailed reference for the classes and functions in the `pyrene-dimer-analyzer` package. The documentation is organized by module.

## Core Module (`pyrene_analyzer.core`)

This module contains the main `PyreneDimerAnalyzer` class that orchestrates the analysis.

### `PyreneDimerAnalyzer`

```python
class PyreneDimerAnalyzer(
    verbose: bool = True,
    use_smarts: bool = True,
    use_shapely: bool = True
)
```

The main class for analyzing pyrene dimer geometric properties.

**Parameters:**

-   `verbose` (bool): If `True`, print progress messages during analysis.
-   `use_smarts` (bool): If `True`, use SMARTS patterns for pyrene detection. If `False`, use a connectivity-based clustering algorithm.
-   `use_shapely` (bool): If `True`, use the Shapely library for accurate π-overlap calculation. If `False`, a centroid-based approximation is used as a fallback.

#### Methods

-   `analyze_file(filename, show_progress=True)`: Analyze all molecules and conformers in a given file.
-   `batch_analyze(filenames, n_jobs=1, show_progress=True)`: Analyze multiple files, with support for parallel processing.
-   `analyze_molecule(mol, mol_name="molecule", show_progress=True)`: Analyze all conformers of a single RDKit molecule object.
-   `analyze_conformer(mol, conf_id, pyrene1_atoms, pyrene2_atoms)`: Analyze a single conformer.
-   `add_classification(df)`: Add a classification column (`strong_excimer`, `weak_excimer`, `monomer`) to a results DataFrame.

## Geometry Module (`pyrene_analyzer.geometry`)

This module provides functions for all geometric calculations.

### Functions

-   `fit_plane_svd(coords)`: Fits a plane to a set of 3D coordinates using SVD.
-   `calculate_plane_angle(coords1, coords2)`: Calculates the angle (0-90°) between two planes.
-   `calculate_interplane_distance(coords1, coords2)`: Calculates the perpendicular distance (Å) between two planes.
-   `calculate_pi_overlap(coords1, coords2, use_shapely=True)`: Calculates the π-π overlap percentage (0-100%) between two ring systems.
-   `calculate_centroid_distance(coords1, coords2)`: Calculates the Euclidean distance (Å) between the centroids of two coordinate sets.
-   `calculate_slip_stack_displacement(coords1, coords2)`: Calculates the lateral (slip-stack) displacement (Å) between two planes.

## I/O Module (`pyrene_analyzer.io`)

This module handles file input and output operations.

### Functions

-   `load_molecules(filename, remove_hydrogens=False)`: Load molecules from a file, automatically detecting the format (SDF, MOL2, PDB).
-   `export_to_csv(df, filename, float_format='%.3f')`: Export a DataFrame to a CSV file.
-   `export_to_json(df, filename, orient='records', indent=2)`: Export a DataFrame to a JSON file.
-   `export_to_excel(df, filename, sheet_name='Analysis')`: Export a DataFrame to an Excel file.
-   `export_results(df, filename, formats=None)`: Export a DataFrame to multiple formats at once.

## Visualization Module (`pyrene_analyzer.visualization`)

This module provides functions for creating publication-quality plots.

### Functions

-   `plot_angle_vs_energy(df, output=None, ...)`: Creates a scatter plot of plane angle vs. relative energy.
-   `plot_distance_vs_overlap(df, output=None, ...)`: Creates a scatter plot of inter-plane distance vs. π-overlap.
-   `plot_conformer_distribution(df, output=None, ...)`: Creates histograms showing the distribution of geometric parameters.
-   `plot_variant_comparison(df, output=None, ...)`: Creates box plots to compare geometric parameters across different R-group variants.
-   `plot_energy_landscape(df, output=None, ...)`: Creates a 2D energy landscape plot.
-   `plot_correlation_matrix(df, output=None, ...)`: Creates a correlation matrix heatmap of geometric parameters.
-   `create_summary_figure(df, output=None, ...)`: Creates a comprehensive summary figure with multiple panels.
