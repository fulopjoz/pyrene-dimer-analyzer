"""Run conformer generation and geometric analysis for all R-group dimer variants.

Uses the 1,1'-binaphthalene dimer SMILES from binaph_dimer_smiles.csv.
Generates conformers, runs geometric analysis, and computes ensemble features.
"""

import os
import sys
import time
import traceback

import pandas as pd
import numpy as np

from pyrene_analyzer.screening import analyze_from_smiles
from pyrene_analyzer.ensemble import compute_ensemble_features


def main():
    # Load dimer SMILES
    df = pd.read_csv("binaph_dimer_smiles.csv")
    print(f"Loaded {len(df)} dimer variants from binaph_dimer_smiles.csv")
    print("=" * 70)

    # Parameters
    NUM_CONFS = 50  # conformers per variant (moderate for screening)
    ENERGY_WINDOW = 10.0  # kcal/mol

    all_conformer_results = []
    summaries = []
    failures = []

    for idx, row in df.iterrows():
        name = row["name"]
        smiles = row["smiles"]
        n_atoms = row["n_atoms"]

        print(f"\n[{idx+1}/{len(df)}] {name} (atoms={n_atoms}, MW={row['mw']})")
        print(f"  SMILES: {smiles[:80]}{'...' if len(smiles) > 80 else ''}")

        t0 = time.perf_counter()
        try:
            results_df, summary = analyze_from_smiles(
                smiles=smiles,
                aromatic_system="binaphthalene",
                num_confs=NUM_CONFS,
                energy_window=ENERGY_WINDOW,
                use_biased=True,
                random_seed=42,
                verbose=False,
            )

            elapsed = time.perf_counter() - t0

            if results_df.empty:
                print(f"  WARNING: No conformers survived for {name} ({elapsed:.1f}s)")
                failures.append(name)
                continue

            # Add substituent name
            results_df["substituent"] = name

            n_conf = len(results_df)
            exc_frac = summary.get("excimer_fraction", 0)
            mean_angle = summary.get("mean_angle", float("nan"))
            mean_dist = summary.get("mean_distance", float("nan"))
            mean_overlap = summary.get("mean_overlap", float("nan"))

            print(f"  {n_conf} conformers, excimer={exc_frac:.1%}, "
                  f"angle={mean_angle:.1f} deg, dist={mean_dist:.2f} A, "
                  f"overlap={mean_overlap:.1f}% ({elapsed:.1f}s)")

            all_conformer_results.append(results_df)
            summary["substituent"] = name
            summary["n_atoms"] = n_atoms
            summary["mw"] = row["mw"]
            summaries.append(summary)

        except Exception as e:
            elapsed = time.perf_counter() - t0
            print(f"  ERROR: {e} ({elapsed:.1f}s)")
            traceback.print_exc()
            failures.append(name)

    # Combine all conformer results
    print("\n" + "=" * 70)
    print("COMBINING RESULTS")
    print("=" * 70)

    if not all_conformer_results:
        print("No results to combine!")
        return

    all_results = pd.concat(all_conformer_results, ignore_index=True)
    all_results.to_csv("binaph_screening_all_conformers.csv", index=False)
    print(f"Saved {len(all_results)} total conformers to binaph_screening_all_conformers.csv")

    # Compute ensemble features
    print("\nComputing ensemble features...")
    features = compute_ensemble_features(
        all_results,
        group_col="substituent",
        energy_col="energy_kcal_mol",
    )
    features.to_csv("binaph_screening_ensemble_features.csv")
    print(f"Saved ensemble features ({len(features)} variants, "
          f"{len(features.columns)} features) to binaph_screening_ensemble_features.csv")

    # Save simple summary
    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv("binaph_screening_summary.csv", index=False)
    print(f"Saved summary to binaph_screening_summary.csv")

    # Print ranking by excimer fraction
    print("\n" + "=" * 70)
    print("SUBSTITUENT RANKING BY EXCIMER FRACTION")
    print("=" * 70)

    # Use ensemble features for ranking
    rank_cols = []
    for col in ["frac_any_excimer", "boltz_excimer_fraction",
                 "plane_angle_deg_boltz", "interplane_distance_A_boltz",
                 "pi_overlap_pct_boltz", "n_conformers"]:
        if col in features.columns:
            rank_cols.append(col)

    if rank_cols:
        ranking = features[rank_cols].copy()
        ranking = ranking.sort_values("frac_any_excimer", ascending=False)
        print(ranking.to_string())

    # Print simple summary ranking
    print("\n" + "=" * 70)
    print("SIMPLE SUMMARY RANKING")
    print("=" * 70)
    print(f"\n{'Rank':<5} {'Name':<20} {'Exc%':>6} {'Angle':>7} {'Dist':>6} "
          f"{'Overlap':>8} {'nConf':>6}")
    print("-" * 63)

    summary_df_sorted = summary_df.sort_values("excimer_fraction", ascending=False)
    for rank, (_, row) in enumerate(summary_df_sorted.iterrows(), 1):
        print(f"{rank:<5} {row['substituent']:<20} "
              f"{row.get('excimer_fraction', 0):>5.1%} "
              f"{row.get('mean_angle', float('nan')):>6.1f} "
              f"{row.get('mean_distance', float('nan')):>6.2f} "
              f"{row.get('mean_overlap', float('nan')):>7.1f} "
              f"{row.get('n_conformers', 0):>6}")

    if failures:
        print(f"\nFailed variants: {', '.join(failures)}")

    print(f"\nTotal: {len(summaries)} successful, {len(failures)} failed")


if __name__ == "__main__":
    main()
