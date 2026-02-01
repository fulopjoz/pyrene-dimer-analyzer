# Example: Batch Processing

This example demonstrates how to analyze multiple SDF files in a batch, which is useful when comparing different R-group variants.

## Using the Python Library

This script will analyze a list of files and combine the results into a single DataFrame. It also shows how to use the visualization module to compare the variants.

```python
import pandas as pd
from pyrene_analyzer import PyreneDimerAnalyzer
from pyrene_analyzer.visualization import plot_variant_comparison

# List of input files for different R-group variants
input_files = [
    "Et_conformers.sdf",
    "iPr_conformers.sdf",
    "cHex_conformers.sdf",
    "tBu_conformers.sdf",
]

# Output file for combined results
output_csv = "all_variants_analysis.csv"
output_plot = "variant_comparison.png"

# Initialize the analyzer
analyzer = PyreneDimerAnalyzer(verbose=True)

# Use batch_analyze for parallel processing
# Set n_jobs to -1 to use all available CPU cores
all_results_df = analyzer.batch_analyze(input_files, n_jobs=-1)

if not all_results_df.empty:
    # Save the combined results
    all_results_df.to_csv(output_csv, index=False)
    print(f"\nCombined results saved to {output_csv}")

    # Create a comparison plot for the plane angle
    plot_variant_comparison(
        all_results_df,
        parameter='plane_angle_deg',
        output=output_plot
    )
    print(f"Comparison plot saved to {output_plot}")

    # Print summary statistics for each variant
    print("\nSummary by Variant:")
    summary = all_results_df.groupby('molecule')['plane_angle_deg'].describe()
    print(summary)

```

## Using the Command-Line Interface

The same batch analysis can be performed easily from the command line.

```bash
# Use a wildcard to select all SDF files
pyrene-analyze analyze *.sdf \
    -o all_variants_analysis.csv \
    --plot \
    --jobs -1
```

This command will:

-   Analyze all files ending with `.sdf` in the current directory.
-   Save the combined results to `all_variants_analysis.csv`.
-   Generate a set of plots in a `plots/` directory.
-   Use all available CPU cores for parallel processing (`--jobs -1`).
