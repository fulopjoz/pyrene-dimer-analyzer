# Quick Start

This guide provides a quick introduction to using `pyrene-dimer-analyzer` for both library and command-line usage.

## Prerequisites

Make sure you have installed the package. If not, please see the [Installation](installation.md) guide.

You will also need a molecular structure file (e.g., in SDF format) containing aromatic dimer conformers. This project includes a sample file at `tests/test_data/pyrene_dimer_set_for_MOE.sdf` that you can use to follow along. Replace it with your own file for real analyses.

## Using the Python Library

Here is a simple example of how to use `pyrene-dimer-analyzer` as a Python library.

```python
# 1. Import the main analyzer class
from pyrene_analyzer import AromaticDimerAnalyzer

# 2. Initialize the analyzer (defaults to pyrene)
analyzer = AromaticDimerAnalyzer(verbose=True)

# 3. Analyze your conformer file
print("Analyzing pyrene dimer conformers...")
results_df = analyzer.analyze_file("tests/test_data/pyrene_dimer_set_for_MOE.sdf")

# 4. Save the results to a CSV file
output_file = "analysis_results.csv"
results_df.to_csv(output_file, index=False)
print(f"Results saved to {output_file}")

# 5. Display the first few rows of the results
print("\nAnalysis Results:")
print(results_df.head())
```

### Multi-System Analysis

You can analyze different aromatic systems by specifying the `aromatic_system` parameter:

```python
from pyrene_analyzer import AromaticDimerAnalyzer

# Analyze with a different aromatic system (uses the same file, different SMARTS)
analyzer = AromaticDimerAnalyzer(aromatic_system="perylene", verbose=True)
results_df = analyzer.analyze_file("tests/test_data/pyrene_dimer_set_for_MOE.sdf")

# Use a custom SMARTS pattern for any aromatic system
analyzer = AromaticDimerAnalyzer(custom_smarts="c1ccc2cc3ccccc3cc2c1")
results_df = analyzer.analyze_file("tests/test_data/pyrene_dimer_set_for_MOE.sdf")
```

Supported built-in systems: `pyrene`, `perylene`, `anthracene`, `naphthalene`, `phenanthrene`.

### Backward Compatibility

The original `PyreneDimerAnalyzer` class name still works:

```python
from pyrene_analyzer import PyreneDimerAnalyzer
analyzer = PyreneDimerAnalyzer()  # defaults to pyrene
```

### Expected Output

The script will produce a `analysis_results.csv` file and print a summary of the results to the console, which will look something like this:

```
[INFO] Analyzing: pyrene_dimer_set_for_MOE
[INFO] ============================================================
[INFO] Identifying aromatic ring systems...
[INFO] Found 2 aromatic systems via SMARTS matching
[INFO] Processing 50 conformers...

Analysis Results:
   molecule  conformer_id  plane_angle_deg  interplane_distance_A  pi_overlap_pct  geometry_warnings ...
0  test_mol             0         75.419640               3.835481       25.418398  High angle...     ...
1  test_mol             1         15.195799               3.478245       69.944132  None              ...
```

Note the `geometry_warnings` column, which flags conformers where theta > 60 deg (inter-plane distance unreliable at high angles). These warnings are recorded in the results rather than printed repeatedly during analysis.

## Using the Command-Line Interface (CLI)

The CLI provides a convenient way to perform the analysis without writing any Python code.

### Basic Analysis

To analyze a single file and save the results to a CSV file:

```bash
pyrene-analyze analyze tests/test_data/pyrene_dimer_set_for_MOE.sdf -o analysis_results.csv
```

### Specifying an Aromatic System

Use `--system` to select a different aromatic system (same file, different detection SMARTS):

```bash
pyrene-analyze analyze tests/test_data/pyrene_dimer_set_for_MOE.sdf -o results.csv --system perylene
```

### Analysis with Plots

The `--plot` flag generates a `plots/` directory with visualizations. You can also analyze multiple files at once if you have them.

```bash
# Analyze with plots
pyrene-analyze analyze tests/test_data/pyrene_dimer_set_for_MOE.sdf -o results.csv --plot

# Batch analysis with multiple files (replace with your own file paths)
# pyrene-analyze analyze Et.sdf iPr.sdf cHex.sdf -o all_variants.csv --plot
```

The `--plot` flag creates a `plots/` directory containing:

- `angle_vs_energy.png`
- `distance_vs_overlap.png`
- `distributions.png`
- `summary.png`

### Listing Supported Systems

To see all registered aromatic systems and their thresholds:

```bash
pyrene-analyze info
```

### Getting Help

To see all available commands and options, use the `--help` flag:

```bash
pyrene-analyze --help

# For help on a specific command
pyrene-analyze analyze --help

# Analyze directly from SMILES (no SDF file needed)
pyrene-analyze analyze-smiles "c1ccc2ccccc2c1CCc1ccc2ccccc2c1" -s naphthalene -n 20

# Screening templates with R-groups (replace with your own template SDF)
# pyrene-analyze screen template.sdf -r "[CH2][CH3]" -o screening_results.csv --plot
```

Note: the `screen` command requires a template SDF file containing your dimer with an R-group to replace. The R-group SMARTS must match a single heavy-atom attachment point.

This quick start should give you a good overview of the capabilities of `pyrene-dimer-analyzer`. For more detailed information, please refer to the [API Reference](api_reference.md) and the [Examples](examples/basic_analysis.md).
