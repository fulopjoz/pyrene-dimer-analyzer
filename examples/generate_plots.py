#!/usr/bin/env python3
"""
Generate Plots
==============

This example demonstrates how to generate publication-quality plots
from analysis results.

Usage:
    python generate_plots.py <results.csv>
"""

import sys
from pathlib import Path

import pandas as pd

from pyrene_analyzer.visualization import (
    plot_angle_vs_energy,
    plot_distance_vs_overlap,
    plot_conformer_distribution,
    plot_variant_comparison,
    create_summary_figure,
)


def main():
    # Check command line arguments
    if len(sys.argv) < 2:
        print("Usage: python generate_plots.py <results.csv>")
        print("\nExample: python generate_plots.py batch_analysis_results.csv")
        sys.exit(1)
    
    input_file = Path(sys.argv[1])
    
    if not input_file.exists():
        print(f"Error: File not found: {input_file}")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path("plots")
    output_dir.mkdir(exist_ok=True)
    
    print(f"Pyrene Dimer Analyzer - Plot Generation")
    print("=" * 50)
    print(f"Input file: {input_file}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Load results
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error loading file: {e}")
        sys.exit(1)
    
    print(f"Loaded {len(df)} rows from {input_file}")
    print()
    
    # Generate plots
    plots_generated = []
    
    # 1. Angle vs Energy (if energy data available)
    if 'rel_energy_kcal_mol' in df.columns:
        print("Generating angle vs energy plot...")
        plot_angle_vs_energy(
            df,
            output=output_dir / "angle_vs_energy.png",
            color_by='molecule' if 'molecule' in df.columns else None
        )
        plots_generated.append("angle_vs_energy.png")
    
    # 2. Distance vs Overlap
    print("Generating distance vs overlap plot...")
    plot_distance_vs_overlap(
        df,
        output=output_dir / "distance_vs_overlap.png",
        color_by='molecule' if 'molecule' in df.columns else None
    )
    plots_generated.append("distance_vs_overlap.png")
    
    # 3. Conformer Distributions
    print("Generating distribution plots...")
    plot_conformer_distribution(
        df,
        output=output_dir / "distributions.png"
    )
    plots_generated.append("distributions.png")
    
    # 4. Variant Comparison (if multiple molecules)
    if 'molecule' in df.columns and df['molecule'].nunique() > 1:
        print("Generating variant comparison plot...")
        plot_variant_comparison(
            df,
            output=output_dir / "variant_comparison.png"
        )
        plots_generated.append("variant_comparison.png")
    
    # 5. Summary Figure
    print("Generating summary figure...")
    create_summary_figure(
        df,
        output=output_dir / "summary.png"
    )
    plots_generated.append("summary.png")
    
    # Print summary
    print("\n" + "=" * 50)
    print("PLOT GENERATION COMPLETE")
    print("=" * 50)
    print(f"\nGenerated {len(plots_generated)} plot(s):")
    for plot in plots_generated:
        print(f"  - {output_dir / plot}")


if __name__ == "__main__":
    main()
