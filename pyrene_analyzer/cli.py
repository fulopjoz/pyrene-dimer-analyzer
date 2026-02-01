"""
Command-Line Interface Module
=============================

Click-based CLI for pyrene dimer conformational analysis.

Usage:
    pyrene-analyze input.sdf -o results.csv
    pyrene-analyze Et.sdf iPr.sdf -o all_results.csv --plot
    pyrene-analyze input.sdf -o results --format csv,json,excel
"""

import sys
from pathlib import Path
from typing import Optional, Tuple

import click
import pandas as pd

from pyrene_analyzer import __version__
from pyrene_analyzer.core import PyreneDimerAnalyzer
from pyrene_analyzer.io import export_to_csv, export_to_json, export_to_excel
from pyrene_analyzer.visualization import (
    plot_angle_vs_energy,
    plot_distance_vs_overlap,
    plot_conformer_distribution,
    create_summary_figure,
)


@click.group(invoke_without_command=True)
@click.option('--version', '-V', is_flag=True, help='Show version and exit.')
@click.pass_context
def cli(ctx: click.Context, version: bool) -> None:
    """
    Pyrene Dimer Analyzer - Automated geometric analysis of pyrene dimer conformers.
    
    Analyzes pyrene dimer molecular geometries including plane-plane angles,
    inter-plane distances, and π-π overlap calculations.
    
    \b
    Examples:
        pyrene-analyze analyze input.sdf -o results.csv
        pyrene-analyze analyze Et.sdf iPr.sdf -o all_results.csv --plot
        pyrene-analyze info
    """
    if version:
        click.echo(f"pyrene-dimer-analyzer version {__version__}")
        ctx.exit()
    
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@cli.command()
@click.argument('input_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('-o', '--output', required=True, type=click.Path(),
              help='Output file path (extension determines format).')
@click.option('-f', '--format', 'formats', default='csv',
              help='Output format(s): csv, json, excel (comma-separated).')
@click.option('-j', '--jobs', default=1, type=int,
              help='Number of parallel jobs (-1 for all CPUs).')
@click.option('-p', '--plot', is_flag=True,
              help='Generate visualization plots.')
@click.option('--plot-dir', type=click.Path(),
              help='Directory for plot output (default: same as output).')
@click.option('-v', '--verbose', is_flag=True,
              help='Enable verbose output.')
@click.option('-q', '--quiet', is_flag=True,
              help='Suppress all output except errors.')
@click.option('--classify', is_flag=True,
              help='Add excimer/monomer classification.')
@click.option('--no-shapely', is_flag=True,
              help='Disable Shapely for π-overlap (use approximation).')
def analyze(
    input_files: Tuple[str, ...],
    output: str,
    formats: str,
    jobs: int,
    plot: bool,
    plot_dir: Optional[str],
    verbose: bool,
    quiet: bool,
    classify: bool,
    no_shapely: bool
) -> None:
    """
    Analyze pyrene dimer conformers from molecular structure files.
    
    \b
    INPUT_FILES: One or more SDF, MOL2, or PDB files containing pyrene dimers.
    
    \b
    Examples:
        pyrene-analyze analyze conformers.sdf -o results.csv
        pyrene-analyze analyze Et.sdf iPr.sdf cHex.sdf -o all_results.csv
        pyrene-analyze analyze input.sdf -o results.csv --plot --verbose
    """
    # Set verbosity
    if quiet:
        verbose = False
    
    # Initialize analyzer
    analyzer = PyreneDimerAnalyzer(
        verbose=verbose,
        use_shapely=not no_shapely
    )
    
    # Process files
    if not quiet:
        click.echo(f"Processing {len(input_files)} file(s)...")
    
    try:
        if len(input_files) == 1:
            results = analyzer.analyze_file(input_files[0])
        else:
            results = analyzer.batch_analyze(list(input_files), n_jobs=jobs)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    
    if results.empty:
        click.echo("Error: No results generated.", err=True)
        sys.exit(1)
    
    # Add classification if requested
    if classify:
        results = analyzer.add_classification(results)
    
    # Export results
    output_path = Path(output)
    format_list = [f.strip().lower() for f in formats.split(',')]
    
    for fmt in format_list:
        if fmt == 'csv':
            csv_path = output_path.with_suffix('.csv')
            export_to_csv(results, csv_path)
            if not quiet:
                click.echo(f"Saved: {csv_path}")
        elif fmt == 'json':
            json_path = output_path.with_suffix('.json')
            export_to_json(results, json_path)
            if not quiet:
                click.echo(f"Saved: {json_path}")
        elif fmt == 'excel':
            excel_path = output_path.with_suffix('.xlsx')
            export_to_excel(results, excel_path)
            if not quiet:
                click.echo(f"Saved: {excel_path}")
        else:
            click.echo(f"Warning: Unknown format '{fmt}'", err=True)
    
    # Generate plots if requested
    if plot:
        if plot_dir:
            plot_path = Path(plot_dir)
        else:
            plot_path = output_path.parent / 'plots'
        
        plot_path.mkdir(parents=True, exist_ok=True)
        
        try:
            plot_angle_vs_energy(
                results,
                output=plot_path / 'angle_vs_energy.png',
                color_by='molecule' if 'molecule' in results.columns else None
            )
            if not quiet:
                click.echo(f"Saved: {plot_path / 'angle_vs_energy.png'}")
            
            plot_distance_vs_overlap(
                results,
                output=plot_path / 'distance_vs_overlap.png',
                color_by='molecule' if 'molecule' in results.columns else None
            )
            if not quiet:
                click.echo(f"Saved: {plot_path / 'distance_vs_overlap.png'}")
            
            plot_conformer_distribution(
                results,
                output=plot_path / 'distributions.png'
            )
            if not quiet:
                click.echo(f"Saved: {plot_path / 'distributions.png'}")
            
            create_summary_figure(
                results,
                output=plot_path / 'summary.png'
            )
            if not quiet:
                click.echo(f"Saved: {plot_path / 'summary.png'}")
            
        except Exception as e:
            click.echo(f"Warning: Could not generate some plots: {e}", err=True)
    
    # Print summary
    if not quiet:
        click.echo(f"\n{'='*60}")
        click.echo("ANALYSIS COMPLETE")
        click.echo(f"{'='*60}")
        click.echo(f"Total conformers analyzed: {len(results)}")
        click.echo(f"Plane angle range: {results['plane_angle_deg'].min():.1f} - "
                   f"{results['plane_angle_deg'].max():.1f}°")
        click.echo(f"Distance range: {results['interplane_distance_A'].min():.2f} - "
                   f"{results['interplane_distance_A'].max():.2f} Å")
        click.echo(f"π-overlap range: {results['pi_overlap_pct'].min():.1f} - "
                   f"{results['pi_overlap_pct'].max():.1f}%")
        
        if classify and 'classification' in results.columns:
            click.echo("\nClassification summary:")
            for cls, count in results['classification'].value_counts().items():
                pct = 100 * count / len(results)
                click.echo(f"  {cls}: {count} ({pct:.1f}%)")


@cli.command()
def info() -> None:
    """
    Display information about the package and its capabilities.
    """
    click.echo(f"""
Pyrene Dimer Analyzer v{__version__}
{'='*50}

A Python package for automated geometric analysis of pyrene dimer conformers.

CAPABILITIES:
  • Automatic pyrene ring detection using SMARTS patterns
  • Plane-plane angle (θ) calculation using SVD
  • Inter-plane distance (d) measurement
  • π-π overlap estimation using Shapely polygon intersection
  • Bridge dihedral angle calculation
  • Batch processing with parallel execution
  • Publication-quality visualization

SUPPORTED FORMATS:
  Input:  SDF, MOL2, PDB
  Output: CSV, JSON, Excel

EXCIMER FORMATION CRITERIA:
  Strong excimer: θ < 20°, d = 3.3-3.7 Å, overlap > 70%
  Weak excimer:   θ = 20-60°, d < 4.5 Å, overlap > 30%
  Monomer:        θ > 60° or d > 4.5 Å or overlap < 30%

REFERENCES:
  • Birks, J.B. (1970). Photophysics of Aromatic Molecules
  • Stevens, B. (1968). Proc. Royal Soc. A, 305, 55-70
  • Poisson, L. et al. (2017). Phys. Chem. Chem. Phys., 19, 23492
  • Ge, Y. et al. (2020). J. Mater. Chem. C, 8, 10223

USAGE:
  pyrene-analyze analyze input.sdf -o results.csv
  pyrene-analyze analyze *.sdf -o results.csv --plot --verbose

For more information, visit:
  https://github.com/research-team/pyrene-dimer-analyzer
""")


@cli.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('-n', '--num', default=5, type=int,
              help='Number of conformers to display.')
def preview(input_file: str, num: int) -> None:
    """
    Preview analysis results for a file without full processing.
    
    Shows the first N conformers with their geometric properties.
    """
    from pyrene_analyzer.io import load_molecules
    
    click.echo(f"Loading: {input_file}")
    
    try:
        molecules = load_molecules(input_file)
    except Exception as e:
        click.echo(f"Error loading file: {e}", err=True)
        sys.exit(1)
    
    click.echo(f"Found {len(molecules)} molecule(s)")
    
    for mol, name in molecules:
        n_confs = mol.GetNumConformers()
        click.echo(f"\nMolecule: {name}")
        click.echo(f"  Atoms: {mol.GetNumAtoms()}")
        click.echo(f"  Conformers: {n_confs}")
        
        if n_confs > 0:
            analyzer = PyreneDimerAnalyzer(verbose=False)
            try:
                pyrene1, pyrene2 = analyzer.identify_pyrene_rings(mol)
                click.echo(f"  Pyrene 1 atoms: {len(pyrene1)}")
                click.echo(f"  Pyrene 2 atoms: {len(pyrene2)}")
                
                click.echo(f"\n  First {min(num, n_confs)} conformers:")
                click.echo(f"  {'ID':>4} {'θ (°)':>8} {'d (Å)':>8} {'Overlap':>8}")
                click.echo(f"  {'-'*32}")
                
                for i in range(min(num, n_confs)):
                    result = analyzer.analyze_conformer(mol, i, pyrene1, pyrene2)
                    click.echo(f"  {i:>4} {result['plane_angle_deg']:>8.1f} "
                               f"{result['interplane_distance_A']:>8.2f} "
                               f"{result['pi_overlap_pct']:>7.1f}%")
            except Exception as e:
                click.echo(f"  Error analyzing: {e}", err=True)


def main() -> None:
    """Main entry point for the CLI."""
    cli()


if __name__ == '__main__':
    main()
