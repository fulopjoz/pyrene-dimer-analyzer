# Example: Basic Analysis

This example demonstrates how to perform a basic analysis of a single SDF file containing pyrene dimer conformers.

## 1. Prepare your Python script

Create a new Python file, for example `run_analysis.py`, and add the following code:

```python
from pyrene_analyzer import PyreneDimerAnalyzer
from pyrene_analyzer.visualization import create_summary_figure

# Define the input file containing your conformers
input_file = "path/to/your/conformers.sdf"

# Define the output file for the results
output_csv = "analysis_results.csv"
output_plot = "analysis_summary.png"

# Initialize the analyzer
analyzer = PyreneDimerAnalyzer(verbose=True)

# Run the analysis
try:
    results_df = analyzer.analyze_file(input_file)

    if not results_df.empty:
        # Save the numerical results to a CSV file
        results_df.to_csv(output_csv, index=False)
        print(f"\nResults saved to {output_csv}")

        # Create and save a summary plot
        create_summary_figure(results_df, output=output_plot)
        print(f"Summary plot saved to {output_plot}")

        # Print a summary of the results
        print("\nAnalysis Summary:")
        print(results_df.describe())

except FileNotFoundError:
    print(f"Error: Input file not found at {input_file}")
except Exception as e:
    print(f"An error occurred: {e}")

```

## 2. Run the script

Execute the script from your terminal:

```bash
python run_analysis.py
```

## 3. Review the output

The script will generate two files:

-   `analysis_results.csv`: A CSV file containing the calculated geometric parameters for every conformer.
-   `analysis_summary.png`: A summary plot showing distributions and correlations of the key parameters.

You can now use this data for further analysis, for example, to build QSAR models or to understand the conformational landscape of your pyrene dimer.
