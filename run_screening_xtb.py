#!/usr/bin/env python
"""
Binaphthalene Dimer Screening with GFN2-xTB
===========================================

This script runs the full 64-variant binaphthalene dimer screening
using GFN2-xTB optimization instead of MMFF94s.

GFN2-xTB correctly reproduces aromatic pi-pi stacking distances
(3.3-3.6 A vs MMFF94s overestimate of 3.8-4.2 A).

Run in WSL with:
    micromamba activate pyrene-xtb
    python run_screening_xtb.py --num-confs 50 --output-prefix binaph_xtb

Expected runtime: 64 variants x 50 conformers x ~15s/conf = ~13 hours
For quick test: python run_screening_xtb.py --num-confs 10 --test
"""

import argparse
import time
from datetime import datetime
from pathlib import Path

import pandas as pd

from pyrene_analyzer.screening import analyze_from_smiles
from pyrene_analyzer.xtb_optimizer import has_xtb_available


def main():
    parser = argparse.ArgumentParser(
        description="Screen binaphthalene dimers with GFN2-xTB optimization"
    )
    parser.add_argument(
        "--smiles-file",
        default="binaph_dimer_smiles.csv",
        help="CSV file with dimer SMILES (default: binaph_dimer_smiles.csv)",
    )
    parser.add_argument(
        "--output-prefix",
        default="binaph_xtb_screening",
        help="Prefix for output files (default: binaph_xtb_screening)",
    )
    parser.add_argument(
        "--num-confs",
        type=int,
        default=50,
        help="Number of conformers per molecule (default: 50)",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run quick test on first 3 molecules only",
    )
    parser.add_argument(
        "--resume",
        type=int,
        default=0,
        help="Resume from molecule index (default: 0 = start fresh)",
    )
    parser.add_argument(
        "--optimizer",
        default="GFN2-xTB",
        choices=["GFN2-xTB", "MMFF94s", "none"],
        help="Optimizer to use (default: GFN2-xTB)",
    )
    parser.add_argument(
        "--priority",
        action="store_true",
        help="Only run priority 10 molecules (EtynPyr_Me, etc.)",
    )
    args = parser.parse_args()

    # Check xtb availability
    if args.optimizer == "GFN2-xTB" and not has_xtb_available():
        print("ERROR: GFN2-xTB requested but xtb-python not available!")
        print("Either install xtb-python or use --optimizer MMFF94s")
        return 1

    print("=" * 60)
    print("BINAPHTHALENE DIMER SCREENING")
    print("=" * 60)
    print(f"Optimizer: {args.optimizer}")
    print(f"Conformers per molecule: {args.num_confs}")
    print(f"Test mode: {args.test}")
    print(f"Resume from: {args.resume}")
    print()

    # Load SMILES
    smiles_file = Path(args.smiles_file)
    if not smiles_file.exists():
        print(f"ERROR: SMILES file not found: {smiles_file}")
        return 1

    df_smiles = pd.read_csv(smiles_file)
    print(f"Loaded {len(df_smiles)} molecules from {smiles_file}")

    if args.test:
        df_smiles = df_smiles.head(3)
        print(f"TEST MODE: Using only first {len(df_smiles)} molecules")

    if args.priority:
        priority_names = [
            "EtynPyr_Me", "EtynPyr_Et", "EtynPyr_H", "EtynPyr_iPr",
            "CNPh_Th_tBu", "CNPh_Th_Me", "Pyr_Me", "Pyr_H",
            "DCV_Th_Me", "DCV_Th_H",
        ]
        df_smiles = df_smiles[df_smiles["name"].isin(priority_names)]
        print(f"PRIORITY MODE: Running {len(df_smiles)} priority molecules")

    # Resume logic
    if args.resume > 0:
        df_smiles = df_smiles.iloc[args.resume:]
        print(f"RESUMING: Starting from molecule index {args.resume}")

    # Output files
    all_results_file = Path(f"{args.output_prefix}_all_conformers.csv")
    summary_file = Path(f"{args.output_prefix}_summary.csv")

    # Track results
    all_conformers = []
    summaries = []

    total_start = time.time()

    for mol_num, (idx, row) in enumerate(df_smiles.iterrows(), 1):
        name = str(row.get("name", f"mol_{idx}"))
        smiles = row.get("smiles", row.get("SMILES", None))

        if not smiles:
            print(f"WARNING: No SMILES for {name}, skipping")
            continue

        print(f"\n[{mol_num}/{len(df_smiles)}] {name}")
        print(f"  SMILES: {smiles[:60]}...")

        mol_start = time.time()

        try:
            results, summary = analyze_from_smiles(
                smiles=smiles,
                aromatic_system="binaphthalene",
                num_confs=args.num_confs,
                optimizer=args.optimizer,
                verbose=True,
            )

            mol_time = time.time() - mol_start

            # Debug: verify return types
            print(f"  DEBUG: results type={type(results).__name__}, "
                  f"summary type={type(summary).__name__}, "
                  f"len(results)={len(results)}")

            # Safely convert results to records and add metadata
            r_group = str(row.get("r_group", ""))
            screen_group = str(row.get("screen_group", ""))

            records = results.to_dict("records")
            for rec in records:
                rec["name"] = name
                rec["r_group"] = r_group
                rec["screen_group"] = screen_group
            all_conformers.extend(records)

            # Safely build summary entry
            summary_entry = dict(summary) if isinstance(summary, dict) else {}
            summary_entry["name"] = name
            summary_entry["r_group"] = r_group
            summary_entry["screen_group"] = screen_group
            summary_entry["time_seconds"] = mol_time
            summaries.append(summary_entry)

            n_confs = summary_entry.get("n_conformers", len(records))
            n_exc = summary_entry.get("n_excimer", 0)
            exc_frac = summary_entry.get("excimer_fraction", 0.0)
            print(f"  Result: {n_confs} conformers, "
                  f"{n_exc} excimer ({exc_frac:.1%})")
            print(f"  Time: {mol_time:.1f}s")

            # Save intermediate results every molecule (not every 5)
            _save_results(all_conformers, summaries, all_results_file, summary_file)
            print(f"  [Saved results: {len(summaries)} molecules, {len(all_conformers)} conformers]")

        except Exception as e:
            import traceback
            print(f"  ERROR: {e}")
            traceback.print_exc()
            continue

    # Final save
    _save_results(all_conformers, summaries, all_results_file, summary_file)

    total_time = time.time() - total_start
    print()
    print("=" * 60)
    print("SCREENING COMPLETE")
    print("=" * 60)
    print(f"Total molecules: {len(summaries)}")
    print(f"Total conformers: {len(all_conformers)}")
    print(f"Total time: {total_time / 3600:.1f} hours")
    print()
    print(f"Results saved to:")
    print(f"  {all_results_file}")
    print(f"  {summary_file}")

    return 0


def _save_results(conformers, summaries, conf_file, summary_file):
    """Save intermediate results to CSV."""
    if conformers:
        pd.DataFrame(conformers).to_csv(conf_file, index=False)
    if summaries:
        pd.DataFrame(summaries).to_csv(summary_file, index=False)


if __name__ == "__main__":
    exit(main())
