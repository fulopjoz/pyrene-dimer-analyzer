#!/usr/bin/env python3
"""
Batch Analysis
==============

This example demonstrates how to analyze multiple SDF files in batch mode
with parallel processing support.

Usage:
    python batch_analysis.py <file1.sdf> <file2.sdf> ...
    python batch_analysis.py *.sdf
"""

import sys
from pathlib import Path

from pyrene_analyzer import PyreneDimerAnalyzer
from pyrene_analyzer.io import export_results


def main():
    # Check command line arguments
    if len(sys.argv) < 2:
        print("Usage: python batch_analysis.py <file1.sdf> <file2.sdf> ...")
        print("\nExample: python batch_analysis.py Et.sdf iPr.sdf cHex.sdf tBu.sdf")
        sys.exit(1)
    
    input_files = [Path(f) for f in sys.argv[1:]]
    
    # Validate files
    valid_files = []
    for f in input_files:
        if f.exists():
            valid_files.append(f)
        else:
            print(f"Warning: File not found, skipping: {f}")
    
    if not valid_files:
        print("Error: No valid input files found.")
        sys.exit(1)
    
    print(f"Pyrene Dimer Analyzer - Batch Analysis")
    print("=" * 50)
    print(f"Processing {len(valid_files)} file(s):")
    for f in valid_files:
        print(f"  - {f}")
    print()
    
    # Initialize the analyzer
    analyzer = PyreneDimerAnalyzer(verbose=True)
    
    # Run batch analysis with parallel processing
    try:
        results = analyzer.batch_analyze(
            [str(f) for f in valid_files],
            n_jobs=-1,  # Use all available CPUs
            show_progress=True
        )
        
        if results.empty:
            print("No results generated. Check your input files.")
            sys.exit(1)
        
        # Add classification
        results = analyzer.add_classification(results)
        
        # Export results to multiple formats
        output_base = "batch_analysis_results"
        export_results(results, output_base, formats=['csv', 'json', 'excel'])
        
        # Print summary by molecule
        print("\n" + "=" * 50)
        print("BATCH ANALYSIS SUMMARY")
        print("=" * 50)
        
        summary = results.groupby('molecule').agg({
            'conformer_id': 'count',
            'plane_angle_deg': ['mean', 'std'],
            'interplane_distance_A': ['mean', 'std'],
            'pi_overlap_pct': ['mean', 'std'],
        })
        
        print("\nResults by molecule:")
        print(summary.to_string())
        
        # Classification summary
        print("\n\nClassification summary:")
        class_summary = results.groupby(['molecule', 'classification']).size().unstack(fill_value=0)
        print(class_summary.to_string())
        
        print(f"\nResults saved to:")
        print(f"  - {output_base}.csv")
        print(f"  - {output_base}.json")
        print(f"  - {output_base}.xlsx")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
