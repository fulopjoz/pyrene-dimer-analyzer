# Quick Start

This guide provides a quick introduction to using `pyrene-dimer-analyzer` for both library and command-line usage.

## Prerequisites

Make sure you have installed the package. If not, please see the [Installation](installation.md) guide.

You will also need a molecular structure file (e.g., in SDF format) containing pyrene dimer conformers. For this guide, we will assume you have a file named `conformers.sdf`.

## Using the Python Library

Here is a simple example of how to use `pyrene-dimer-analyzer` as a Python library.

```python
# 1. Import the main analyzer class
from pyrene_analyzer import PyreneDimerAnalyzer

# 2. Initialize the analyzer
# Verbose mode will print progress messages.
analyzer = PyreneDimerAnalyzer(verbose=True)

# 3. Analyze your conformer file
print("Analyzing conformers.sdf...")
results_df = analyzer.analyze_file("conformers.sdf")

# 4. Save the results to a CSV file
output_file = "analysis_results.csv"
results_df.to_csv(output_file, index=False)
print(f"Results saved to {output_file}")

# 5. Display the first few rows of the results
print("\nAnalysis Results:")
print(results_df.head())
```

### Expected Output

The script will produce a `analysis_results.csv` file and print a summary of the results to the console, which will look something like this:

```
[INFO] Analyzing: conformers
[INFO] ============================================================
[INFO] Identifying pyrene ring systems...
[INFO] Found 2 pyrene systems via SMARTS matching
[INFO] Processing 50 conformers...

Analysis Results:
   molecule  conformer_id  plane_angle_deg  interplane_distance_A  pi_overlap_pct  ... 
0  test_mol             0         75.419640               3.835481       25.418398  ... 
1  test_mol             1         80.195799               4.178245       19.944132  ... 
2  test_mol             2         70.795811               3.999123       31.264911  ... 
3  test_mol             3         85.139938               4.315386       15.267411  ... 
4  test_mol             4         65.942821               4.012345       35.987123  ... 
```

## Using the Command-Line Interface (CLI)

The CLI provides a convenient way to perform the analysis without writing any Python code.

### Basic Analysis

To analyze a single file and save the results to a CSV file:

```bash
pyrene-analyze analyze conformers.sdf -o analysis_results.csv
```

### Batch Analysis with Plots

You can analyze multiple files at once and generate plots automatically. The `--plot` flag will create a `plots/` directory with several useful visualizations.

```bash
# Analyze multiple R-group variants and generate plots
pyrene-analyze analyze Et.sdf iPr.sdf cHex.sdf -o all_variants.csv --plot
```

This will create `all_variants.csv` and a `plots/` directory containing:
-   `angle_vs_energy.png`
-   `distance_vs_overlap.png`
-   `distributions.png`
-   `summary.png`

### Getting Help

To see all available commands and options, use the `--help` flag:

```bash
pyrene-analyze --help

# For help on a specific command
pyrene-analyze analyze --help
```

This quick start should give you a good overview of the capabilities of `pyrene-dimer-analyzer`. For more detailed information, please refer to the [API Reference](api_reference.md) and the [Examples](examples/basic_analysis.md).
