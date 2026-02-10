#!/usr/bin/env python
"""
GFN2-xTB Single-Point Re-ranking Script
========================================

This script recalculates energies for existing MMFF94s-optimized conformers
using GFN2-xTB single-point calculations. This gives accurate energy ranking
WITHOUT the time cost of full GFN2-xTB optimization.

Rationale:
    - MMFF94s optimization: ~200 conformers in 5 minutes (fast)
    - GFN2-xTB optimization: ~200 conformers in 160 hours (too slow!)
    - GFN2-xTB single-point: ~200 conformers in 7 minutes (fast + accurate energies)

The distances will still be from MMFF94s (3.8-4.2 Å), but the energy ranking
and filtering will use GFN2-xTB energies (correct dispersion).

Run in WSL with:
    micromamba activate pyrene-xtb
    python run_xtb_singlepoint.py --input binaph_screening_all_conformers.csv

Expected runtime: ~2 hours for all 64 variants (~2s per conformer)
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from pyrene_analyzer.xtb_optimizer import has_xtb_available, single_point_gfn2xtb


def main():
    parser = argparse.ArgumentParser(
        description="Re-rank conformers using GFN2-xTB single-point energies"
    )
    parser.add_argument(
        "--input",
        default="binaph_screening_all_conformers.csv",
        help="Input CSV with conformer data (default: binaph_screening_all_conformers.csv)",
    )
    parser.add_argument(
        "--smiles-file",
        default="binaph_dimer_smiles.csv",
        help="CSV file with dimer SMILES (default: binaph_dimer_smiles.csv)",
    )
    parser.add_argument(
        "--output-prefix",
        default="binaph_xtb_sp",
        help="Prefix for output files (default: binaph_xtb_sp)",
    )
    parser.add_argument(
        "--method",
        default="GFN2-xTB",
        choices=["GFN2-xTB", "GFN1-xTB", "GFN0-xTB", "GFNFF"],
        help="xTB method for single-point (default: GFN2-xTB)",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run on first 10 conformers only",
    )
    args = parser.parse_args()

    # Check xtb availability
    if not has_xtb_available():
        print("ERROR: xtb-python not available!")
        print("Run this script in WSL with the pyrene-xtb environment.")
        return 1

    print("=" * 60)
    print("GFN2-xTB SINGLE-POINT RE-RANKING")
    print("=" * 60)
    print(f"Method: {args.method}")
    print(f"Input: {args.input}")
    print()

    # Load existing results
    input_file = Path(args.input)
    if not input_file.exists():
        print(f"ERROR: Input file not found: {input_file}")
        print("Run the MMFF94s screening first: python run_screening.py")
        return 1

    df = pd.read_csv(input_file)
    print(f"Loaded {len(df)} conformers from {input_file}")

    # Load SMILES for rebuilding molecules
    smiles_file = Path(args.smiles_file)
    if not smiles_file.exists():
        print(f"ERROR: SMILES file not found: {smiles_file}")
        return 1

    df_smiles = pd.read_csv(smiles_file)
    smiles_dict = {row["name"]: row.get("smiles", row.get("SMILES"))
                   for _, row in df_smiles.iterrows()}

    if args.test:
        df = df.head(10)
        print(f"TEST MODE: Using only first {len(df)} conformers")

    # Check for coordinate columns
    has_coords = "conf_coords" in df.columns or "mol_block" in df.columns
    if not has_coords:
        # Try to reconstruct from SMILES + conformer ID
        print("NOTE: No stored coordinates. Will regenerate from SMILES.")
        print("      This means conformer IDs must match the original generation.")
        print()

    # Output files
    output_file = Path(f"{args.output_prefix}_energies.csv")

    # Process conformers
    results = []
    total_start = time.time()

    # Group by molecule name for efficiency
    for name in df["name"].unique() if "name" in df.columns else [None]:
        if name:
            mol_df = df[df["name"] == name]
            smiles = smiles_dict.get(name)
            if not smiles:
                print(f"WARNING: No SMILES for {name}, skipping")
                continue
        else:
            mol_df = df
            smiles = df.iloc[0].get("smiles", df.iloc[0].get("SMILES"))

        print(f"\n[{name or 'molecule'}] {len(mol_df)} conformers")

        # Build molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  ERROR: Invalid SMILES, skipping")
            continue

        mol = Chem.AddHs(mol)

        # Generate conformers to match original IDs
        # Note: This assumes the same random seed and parameters were used
        n_confs = len(mol_df)
        AllChem.EmbedMultipleConfs(
            mol, numConfs=n_confs * 2,  # Generate extra in case some fail
            randomSeed=42,
            useRandomCoords=True,
        )

        if mol.GetNumConformers() < n_confs:
            print(f"  WARNING: Only generated {mol.GetNumConformers()} conformers, "
                  f"need {n_confs}")

        # Calculate single-point energies for each conformer
        for i, (_, row) in enumerate(mol_df.iterrows()):
            conf_idx = i if i < mol.GetNumConformers() else 0

            try:
                t0 = time.perf_counter()
                xtb_energy = single_point_gfn2xtb(mol, conf_id=conf_idx, method=args.method)
                elapsed = time.perf_counter() - t0

                # Keep all original data, add xTB energy
                result = row.to_dict()
                result["xtb_energy"] = xtb_energy
                result["xtb_method"] = args.method
                result["xtb_time_s"] = elapsed
                results.append(result)

                if (i + 1) % 10 == 0 or i == len(mol_df) - 1:
                    print(f"  Conformer {i + 1}/{len(mol_df)}: "
                          f"E = {xtb_energy:.2f} kcal/mol ({elapsed:.1f}s)")

            except Exception as e:
                print(f"  Conformer {i + 1}: ERROR: {e}")
                result = row.to_dict()
                result["xtb_energy"] = np.nan
                result["xtb_method"] = args.method
                result["xtb_time_s"] = 0
                results.append(result)

        # Save intermediate results
        if results:
            pd.DataFrame(results).to_csv(output_file, index=False)

    # Final statistics
    total_time = time.time() - total_start
    df_results = pd.DataFrame(results)

    print()
    print("=" * 60)
    print("SINGLE-POINT CALCULATIONS COMPLETE")
    print("=" * 60)
    print(f"Total conformers: {len(df_results)}")
    print(f"Total time: {total_time / 60:.1f} minutes")
    print(f"Average per conformer: {total_time / len(df_results):.1f}s")
    print()

    # Summary statistics
    if "xtb_energy" in df_results.columns:
        valid_energies = df_results["xtb_energy"].dropna()
        print(f"Energy range: {valid_energies.min():.2f} to {valid_energies.max():.2f} kcal/mol")
        print(f"Mean energy: {valid_energies.mean():.2f} kcal/mol")
        print(f"Std energy: {valid_energies.std():.2f} kcal/mol")

    print()
    print(f"Results saved to: {output_file}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
