#!/usr/bin/env python3
"""
Analyze Single Molecule
=======================

This example demonstrates how to analyze a single SDF file containing
pyrene dimer conformers using the pyrene-dimer-analyzer package.

Usage:
    python analyze_single_molecule.py <input_file.sdf>
"""

import sys
from pathlib import Path

from pyrene_analyzer import PyreneDimerAnalyzer
from pyrene_analyzer.io import export_to_csv


def main():
    # Check command line arguments
    if len(sys.argv) < 2:
        print("Usage: python analyze_single_molecule.py <input_file.sdf>")
        print("\nExample: python analyze_single_molecule.py conformers.sdf")
        sys.exit(1)
    
    input_file = Path(sys.argv[1])
    
    if not input_file.exists():
        print(f"Error: File not found: {input_file}")
        sys.exit(1)
    
    # Define output file
    output_file = input_file.stem + "_analysis.csv"
    
    print(f"Pyrene Dimer Analyzer - Single Molecule Analysis")
    print("=" * 50)
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    print()
    
    # Initialize the analyzer
    analyzer = PyreneDimerAnalyzer(verbose=True)
    
    # Run analysis
    try:
        results = analyzer.analyze_file(input_file)
        
        if results.empty:
            print("No results generated. Check your input file.")
            sys.exit(1)
        
        # Add classification
        results = analyzer.add_classification(results)
        
        # Export results
        export_to_csv(results, output_file)
        
        # Print summary
        print("\n" + "=" * 50)
        print("ANALYSIS SUMMARY")
        print("=" * 50)
        print(f"Total conformers analyzed: {len(results)}")
        print(f"\nPlane angle (θ):")
        print(f"  Mean: {results['plane_angle_deg'].mean():.1f}°")
        print(f"  Min:  {results['plane_angle_deg'].min():.1f}°")
        print(f"  Max:  {results['plane_angle_deg'].max():.1f}°")
        print(f"\nInter-plane distance (d):")
        print(f"  Mean: {results['interplane_distance_A'].mean():.2f} Å")
        print(f"  Min:  {results['interplane_distance_A'].min():.2f} Å")
        print(f"  Max:  {results['interplane_distance_A'].max():.2f} Å")
        print(f"\nπ-overlap:")
        print(f"  Mean: {results['pi_overlap_pct'].mean():.1f}%")
        print(f"  Min:  {results['pi_overlap_pct'].min():.1f}%")
        print(f"  Max:  {results['pi_overlap_pct'].max():.1f}%")
        
        # Classification summary
        print(f"\nClassification:")
        for cls, count in results['classification'].value_counts().items():
            pct = 100 * count / len(results)
            print(f"  {cls}: {count} ({pct:.1f}%)")
        
        print(f"\nResults saved to: {output_file}")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
