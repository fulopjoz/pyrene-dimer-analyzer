#!/usr/bin/env python
"""
Analyze MOE Conformer Search Results
=====================================

Maps MOE conformers to molecule names using the mseq property,
runs geometric analysis, classifies conformers, and generates
per-molecule summary statistics.

Input:
    - moe_conformers/cnph_th_cf3_3d_conformers.sdf (MOE V3000 output)
    - binaph_dimer_smiles.csv (molecule name/metadata mapping)

Output:
    - moe_screening_all_conformers.csv (per-conformer analysis)
    - moe_screening_summary.csv (per-molecule summary)
"""

import argparse
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

from pyrene_analyzer.core import AromaticDimerAnalyzer


def load_molecule_map(smiles_csv: str) -> dict:
    """Load mseq -> molecule metadata mapping from binaph_dimer_smiles.csv.

    mseq is 1-indexed and maps to CSV row order.
    """
    df = pd.read_csv(smiles_csv)
    mol_map = {}
    for i, row in df.iterrows():
        mseq = i + 1  # 1-indexed
        mol_map[mseq] = {
            "name": row["name"],
            "r_group": row["r_group"],
            "screen_group": row["screen_group"],
        }
    return mol_map


def main():
    parser = argparse.ArgumentParser(
        description="Analyze MOE conformer search results"
    )
    parser.add_argument(
        "--sdf",
        default="moe_conformers/cnph_th_cf3_3d_conformers.sdf",
        help="MOE conformer SDF file",
    )
    parser.add_argument(
        "--smiles-csv",
        default="binaph_dimer_smiles.csv",
        help="Molecule name mapping CSV",
    )
    parser.add_argument(
        "--output-prefix",
        default="moe_screening",
        help="Prefix for output files",
    )
    args = parser.parse_args()

    sdf_path = Path(args.sdf)
    smiles_csv = Path(args.smiles_csv)

    if not sdf_path.exists():
        print(f"ERROR: SDF file not found: {sdf_path}")
        return 1
    if not smiles_csv.exists():
        print(f"ERROR: SMILES CSV not found: {smiles_csv}")
        return 1

    # Load molecule name mapping
    mol_map = load_molecule_map(str(smiles_csv))
    print(f"Loaded {len(mol_map)} molecule mappings from {smiles_csv}")

    # Initialize analyzer for binaphthalene
    analyzer = AromaticDimerAnalyzer(aromatic_system="binaphthalene")

    # Read SDF entries one by one
    print(f"\nReading SDF: {sdf_path}")
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)

    all_results = []
    errors = []
    ring_cache = {}  # Cache aromatic ring detection per mseq
    start_time = time.time()

    for entry_idx, mol in enumerate(supplier):
        if mol is None:
            errors.append({"entry": entry_idx, "error": "RDKit could not parse"})
            continue

        # Extract mseq
        if mol.HasProp("mseq"):
            mseq = int(mol.GetProp("mseq"))
        else:
            errors.append({"entry": entry_idx, "error": "No mseq property"})
            continue

        # Get molecule metadata
        meta = mol_map.get(mseq)
        if meta is None:
            errors.append({"entry": entry_idx, "error": f"Unknown mseq={mseq}"})
            continue

        mol_name = meta["name"]

        # Extract MOE energy properties
        energy = None
        dE = None
        for prop in ["E", "energy", "Energy"]:
            if mol.HasProp(prop):
                try:
                    energy = float(mol.GetProp(prop))
                    break
                except ValueError:
                    pass
        if mol.HasProp("dE"):
            try:
                dE = float(mol.GetProp("dE"))
            except ValueError:
                pass

        # Detect aromatic rings (cache per mseq since all conformers of same
        # molecule have same topology)
        try:
            if mseq not in ring_cache:
                pyrene1, pyrene2 = analyzer.identify_aromatic_rings(mol)
                ring_cache[mseq] = (pyrene1, pyrene2)
            else:
                pyrene1, pyrene2 = ring_cache[mseq]
        except ValueError as e:
            errors.append({
                "entry": entry_idx, "mseq": mseq,
                "name": mol_name, "error": str(e)
            })
            continue

        # Analyze geometry (each MOE entry has exactly 1 conformer at id=0)
        try:
            result = analyzer.analyze_conformer(mol, 0, pyrene1, pyrene2)
        except Exception as e:
            errors.append({
                "entry": entry_idx, "mseq": mseq,
                "name": mol_name, "error": str(e)
            })
            continue

        # Classify
        classification = analyzer.classify_conformer(
            result["plane_angle_deg"],
            result["interplane_distance_A"],
            result["pi_overlap_pct"],
        )

        # Build result record
        record = {
            "name": mol_name,
            "r_group": meta["r_group"],
            "screen_group": meta["screen_group"],
            "mseq": mseq,
            "conformer_id": result["conformer_id"],
            "plane_angle_deg": result["plane_angle_deg"],
            "interplane_distance_A": result["interplane_distance_A"],
            "pi_overlap_pct": result["pi_overlap_pct"],
            "centroid_distance_A": result["centroid_distance_A"],
            "slip_stack_A": result["slip_stack_A"],
            "bridge_dihedral_L_deg": result.get("bridge_dihedral_L_deg"),
            "bridge_dihedral_R_deg": result.get("bridge_dihedral_R_deg"),
            "energy_kcal_mol": energy,
            "dE_kcal_mol": dE,
            "geometry_warnings": result.get("geometry_warnings"),
            "classification": classification,
        }
        all_results.append(record)

        # Progress
        if (entry_idx + 1) % 500 == 0:
            elapsed = time.time() - start_time
            print(f"  Processed {entry_idx + 1} conformers "
                  f"({elapsed:.1f}s, {len(errors)} errors)")

    elapsed = time.time() - start_time
    print(f"\nProcessed {len(all_results)} conformers in {elapsed:.1f}s "
          f"({len(errors)} errors)")

    if errors:
        print(f"\nFirst 5 errors:")
        for e in errors[:5]:
            print(f"  {e}")

    if not all_results:
        print("ERROR: No results to save!")
        return 1

    # Create all-conformers DataFrame
    df_all = pd.DataFrame(all_results)

    # Add relative energy per molecule
    for name in df_all["name"].unique():
        mask = df_all["name"] == name
        energies = df_all.loc[mask, "energy_kcal_mol"]
        if energies.notna().any():
            df_all.loc[mask, "rel_energy_kcal_mol"] = energies - energies.min()

    # Save all conformers
    all_conf_file = f"{args.output_prefix}_all_conformers.csv"
    df_all.to_csv(all_conf_file, index=False)
    print(f"\nSaved {len(df_all)} conformers to {all_conf_file}")

    # Generate per-molecule summary
    summaries = []
    for name in df_all["name"].unique():
        df_mol = df_all[df_all["name"] == name]
        meta = {"name": name}
        meta["r_group"] = df_mol["r_group"].iloc[0]
        meta["screen_group"] = df_mol["screen_group"].iloc[0]
        meta["mseq"] = df_mol["mseq"].iloc[0]
        meta["n_conformers"] = len(df_mol)

        # Classification counts
        n_strong = (df_mol["classification"] == "strong_excimer").sum()
        n_weak = (df_mol["classification"] == "weak_excimer").sum()
        n_monomer = (df_mol["classification"] == "monomer").sum()
        n_excimer = n_strong + n_weak
        meta["n_strong_excimer"] = n_strong
        meta["n_weak_excimer"] = n_weak
        meta["n_monomer"] = n_monomer
        meta["n_excimer"] = n_excimer
        meta["excimer_fraction"] = n_excimer / len(df_mol) if len(df_mol) > 0 else 0

        # Geometric stats
        meta["mean_angle_deg"] = df_mol["plane_angle_deg"].mean()
        meta["mean_distance_A"] = df_mol["interplane_distance_A"].mean()
        meta["mean_overlap_pct"] = df_mol["pi_overlap_pct"].mean()
        meta["min_distance_A"] = df_mol["interplane_distance_A"].min()
        meta["max_overlap_pct"] = df_mol["pi_overlap_pct"].max()

        # Excimer-only stats (if any)
        df_exc = df_mol[df_mol["classification"].str.contains("excimer")]
        if len(df_exc) > 0:
            meta["excimer_mean_angle_deg"] = df_exc["plane_angle_deg"].mean()
            meta["excimer_mean_distance_A"] = df_exc["interplane_distance_A"].mean()
            meta["excimer_mean_overlap_pct"] = df_exc["pi_overlap_pct"].mean()

        # Energy stats
        if df_mol["energy_kcal_mol"].notna().any():
            meta["min_energy_kcal_mol"] = df_mol["energy_kcal_mol"].min()
            meta["energy_range_kcal_mol"] = (
                df_mol["energy_kcal_mol"].max() - df_mol["energy_kcal_mol"].min()
            )

        summaries.append(meta)

    df_summary = pd.DataFrame(summaries)
    summary_file = f"{args.output_prefix}_summary.csv"
    df_summary.to_csv(summary_file, index=False)
    print(f"Saved {len(df_summary)} molecule summaries to {summary_file}")

    # Print overview
    print(f"\n{'='*60}")
    print("MOE CONFORMER ANALYSIS SUMMARY")
    print(f"{'='*60}")
    print(f"Total molecules: {len(df_summary)}")
    print(f"Total conformers: {len(df_all)}")

    # Per R-group summary
    print(f"\nPer R-group:")
    for rg in df_summary["r_group"].unique():
        rg_data = df_summary[df_summary["r_group"] == rg]
        total_confs = rg_data["n_conformers"].sum()
        total_exc = rg_data["n_excimer"].sum()
        n_with_exc = (rg_data["n_excimer"] > 0).sum()
        conf_range = f"{rg_data['n_conformers'].min()}-{rg_data['n_conformers'].max()}"
        print(f"  {rg:12s}: {len(rg_data)} mols, {total_confs:4d} confs "
              f"(range {conf_range}), {total_exc} excimer, "
              f"{n_with_exc}/{len(rg_data)} variants with excimer")

    # Overall classification
    print(f"\nClassification:")
    for cls in ["strong_excimer", "weak_excimer", "monomer"]:
        n = (df_all["classification"] == cls).sum()
        pct = n / len(df_all) * 100
        print(f"  {cls:16s}: {n:5d} ({pct:.1f}%)")

    # Top excimer molecules
    top = df_summary.nlargest(10, "excimer_fraction")
    if top["excimer_fraction"].max() > 0:
        print(f"\nTop excimer-forming molecules:")
        for _, row in top.iterrows():
            if row["excimer_fraction"] > 0:
                print(f"  {row['name']:20s}: {row['excimer_fraction']:.1%} "
                      f"({row['n_excimer']}/{row['n_conformers']} conformers)")

    # Conformer count distribution
    print(f"\nConformer count per molecule:")
    print(f"  Min: {df_summary['n_conformers'].min()}")
    print(f"  Max: {df_summary['n_conformers'].max()}")
    print(f"  Mean: {df_summary['n_conformers'].mean():.1f}")
    print(f"  Median: {df_summary['n_conformers'].median():.0f}")

    # Molecules with very few conformers
    sparse = df_summary[df_summary["n_conformers"] <= 5]
    if len(sparse) > 0:
        print(f"\n  WARNING: {len(sparse)} molecules with <= 5 conformers:")
        for _, row in sparse.iterrows():
            print(f"    {row['name']}: {row['n_conformers']} conformers")

    return 0


if __name__ == "__main__":
    exit(main())
