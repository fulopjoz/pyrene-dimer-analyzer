"""
Re-optimize MOE conformers with MACE-OFF23.

Takes 3,347 MOE conformers (Amber:EHT force field) and re-optimizes
geometries using MACE-OFF23 (DFT-quality dispersion). This corrects the
systematic distance overestimation of classical force fields for pi-stacking.

Scientific basis:
    MOE/Amber:EHT overestimates pi-stacking distances by ~1-2 A due to
    incomplete London dispersion treatment. MACE-OFF23 (trained on
    wB97M-D3BJ/def2-TZVPPD) reproduces CCSD(T)/CBS benchmark distances
    (3.43 A for pyrene dimer, PMC11476719 2024).

    Expected effect: conformers at 4.5-5.5 A (force field artifact) should
    compress to 3.3-4.2 A after MACE-OFF23 optimization.

Usage:
    python reoptimize_moe_conformers.py --test
    python reoptimize_moe_conformers.py --model small --device auto
    python reoptimize_moe_conformers.py --resume --output-prefix moe_mace_reopt

References:
    - Batatia et al. (2024) arXiv:2401.00096 (MACE-OFF23)
    - PMC11476719 (2024): DLPNO-CCSD(T)/CBS pyrene dimer d=3.43 A
    - Ge et al. (2020) J. Mater. Chem. C 8, 10223 (4.5 A threshold)
"""

import argparse
import sys
import time
import traceback
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from pyrene_analyzer.core import AromaticDimerAnalyzer
from pyrene_analyzer.mace_optimizer import (
    has_mace_available,
    get_mace_import_error,
    optimize_with_mace,
)


def load_molecule_map(smiles_csv: str) -> dict:
    """Load molecule name mapping from binaph_dimer_smiles.csv.

    Returns dict mapping mseq (1-indexed) to molecule info.
    """
    df = pd.read_csv(smiles_csv)
    mol_map = {}
    for i, row in df.iterrows():
        mseq = i + 1  # 1-indexed to match MOE mseq property
        mol_map[mseq] = {
            "name": row["name"],
            "r_group": row["r_group"],
            "screen_group": row["screen_group"],
        }
    return mol_map


def process_conformer(
    mol: Chem.Mol,
    conf_id: int,
    analyzer: AromaticDimerAnalyzer,
    ring1_atoms: list,
    ring2_atoms: list,
    model: str = "small",
    device: str = "cpu",
    mace_dtype: str = "float64",
    max_steps: int = 500,
    fmax: float = 0.005,
) -> dict:
    """Optimize a single conformer with MACE-OFF23 and analyze geometry.

    Returns dict with both pre- and post-optimization geometry.
    """
    # Pre-optimization geometry
    pre_result = analyzer.analyze_conformer(mol, conf_id, ring1_atoms, ring2_atoms)

    # Optimize with MACE-OFF23
    opt_mol, energy = optimize_with_mace(
        mol,
        conf_id=conf_id,
        model=model,
        device=device,
        default_dtype=mace_dtype,
        max_steps=max_steps,
        fmax=fmax,
        verbose=False,
    )

    # Post-optimization geometry
    post_result = analyzer.analyze_conformer(opt_mol, conf_id, ring1_atoms, ring2_atoms)

    return {
        "mace_energy_kcal": energy,
        # Pre-optimization (MOE) geometry
        "moe_plane_angle_deg": pre_result["plane_angle_deg"],
        "moe_interplane_distance_A": pre_result["interplane_distance_A"],
        "moe_pi_overlap_pct": pre_result["pi_overlap_pct"],
        # Post-optimization (MACE) geometry
        "plane_angle_deg": post_result["plane_angle_deg"],
        "interplane_distance_A": post_result["interplane_distance_A"],
        "pi_overlap_pct": post_result["pi_overlap_pct"],
        "centroid_distance_A": post_result.get("centroid_distance_A", np.nan),
        "slip_stack_displacement_A": post_result.get("slip_stack_displacement_A", np.nan),
        "bridge_dihedral_deg": post_result.get("bridge_dihedral_deg", np.nan),
        "geometry_warnings": post_result.get("geometry_warnings", ""),
    }


def resolve_device(device_arg: str) -> str:
    """Resolve requested device to an available backend."""
    requested = device_arg.lower()
    if requested not in {"auto", "cpu", "cuda"}:
        raise ValueError(f"Unsupported device '{device_arg}'. Use auto/cpu/cuda.")

    if requested == "cpu":
        return "cpu"

    try:
        import torch
    except ImportError:
        if requested == "cuda":
            raise RuntimeError(
                "CUDA requested but torch is not importable in this environment."
            )
        return "cpu"

    if requested == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"

    if requested == "cuda" and not torch.cuda.is_available():
        raise RuntimeError(
            "CUDA requested but no CUDA device is available (torch.cuda.is_available() is False)."
        )
    return requested


def main():
    parser = argparse.ArgumentParser(
        description="Re-optimize MOE conformers with MACE-OFF23"
    )
    parser.add_argument(
        "--sdf",
        default="moe_conformers/cnph_th_cf3_3d_conformers.sdf",
        help="SDF file with MOE conformers (default: moe_conformers/...)",
    )
    parser.add_argument(
        "--smiles-csv",
        default="binaph_dimer_smiles.csv",
        help="CSV with molecule names and mseq mapping",
    )
    parser.add_argument(
        "--output-prefix",
        default="moe_mace_reopt",
        help="Output file prefix (default: moe_mace_reopt)",
    )
    parser.add_argument(
        "--model",
        default="small",
        choices=["small", "medium", "large"],
        help="MACE model size (default: small, best throughput)",
    )
    parser.add_argument(
        "--device",
        default="auto",
        choices=["auto", "cpu", "cuda"],
        help="Device for MACE (auto|cpu|cuda, default: auto)",
    )
    parser.add_argument(
        "--mace-dtype",
        default="float64",
        choices=["float64", "float32"],
        help="MACE dtype (float64=more accurate, float32=faster GPU)",
    )
    parser.add_argument(
        "--max-steps",
        default=500,
        type=int,
        help="Max optimization steps per conformer (default: 500)",
    )
    parser.add_argument(
        "--fmax",
        default=0.005,
        type=float,
        help="Force convergence threshold (eV/A, default: 0.005)",
    )
    parser.add_argument(
        "--save-every",
        default=25,
        type=int,
        help="Write intermediate CSV every N processed conformers (default: 25)",
    )
    parser.add_argument(
        "--start-index",
        default=0,
        type=int,
        help="Global conformer start index for chunked runs (default: 0)",
    )
    parser.add_argument(
        "--end-index",
        default=None,
        type=int,
        help="Global conformer end index (exclusive) for chunked runs",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Test mode: process only first 3 conformers",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from existing output file (skip already processed)",
    )
    args = parser.parse_args()

    if args.save_every < 1:
        parser.error("--save-every must be >= 1")
    if args.start_index < 0:
        parser.error("--start-index must be >= 0")
    if args.end_index is not None and args.end_index <= args.start_index:
        parser.error("--end-index must be greater than --start-index")

    # Check MACE availability
    if not has_mace_available():
        print(f"ERROR: MACE-OFF23 not available: {get_mace_import_error()}")
        print("Install with: pip install mace-torch ase torch")
        sys.exit(1)

    try:
        resolved_device = resolve_device(args.device)
    except Exception as e:
        print(f"ERROR: could not resolve device: {e}")
        sys.exit(1)

    print(
        f"Runtime config: model={args.model}, device={resolved_device}, "
        f"dtype={args.mace_dtype}, max_steps={args.max_steps}, fmax={args.fmax}"
    )

    # Load molecule name mapping
    mol_map = load_molecule_map(args.smiles_csv)
    print(f"Loaded {len(mol_map)} molecule definitions from {args.smiles_csv}")

    # Load existing results for resume
    output_csv = f"{args.output_prefix}_all_conformers.csv"
    processed_keys = set()
    existing_records = []
    if args.resume and Path(output_csv).exists():
        existing_df = pd.read_csv(output_csv)
        for _, row in existing_df.iterrows():
            key = (row["name"], row.get("conformer_id", row.get("conf_idx", -1)))
            processed_keys.add(key)
        existing_records = existing_df.to_dict("records")
        print(f"Resuming: {len(processed_keys)} conformers already processed")

    # Load SDF
    print(f"Loading conformers from {args.sdf}...")
    supplier = Chem.SDMolSupplier(args.sdf, removeHs=False, sanitize=True)
    if supplier is None:
        # Try V3000 format
        mols = []
        with open(args.sdf, "r") as f:
            content = f.read()
        blocks = content.split("$$$$")
        for block in blocks:
            block = block.strip()
            if not block:
                continue
            mol = Chem.MolFromMolBlock(block, removeHs=False, sanitize=True)
            if mol is not None:
                mols.append(mol)
        print(f"Loaded {len(mols)} conformers (V3000 format)")
    else:
        mols = [m for m in supplier if m is not None]
        print(f"Loaded {len(mols)} conformers")

    indexed_mols = list(enumerate(mols))
    start = args.start_index
    end: Optional[int] = args.end_index
    indexed_mols = indexed_mols[start:end]

    if args.test:
        indexed_mols = indexed_mols[:3]
        print(f"TEST MODE: processing only {len(indexed_mols)} conformers")

    print(
        f"Processing window: start={start}, "
        f"end={'EOF' if end is None else end}, selected={len(indexed_mols)}"
    )

    # Initialize analyzer
    analyzer = AromaticDimerAnalyzer(
        aromatic_system="binaphthalene", verbose=False
    )

    # Process conformers
    all_records = list(existing_records)
    n_total = len(indexed_mols)
    n_errors = 0
    n_processed = 0
    t_start = time.perf_counter()

    # Cache ring detection per molecule (mseq)
    ring_cache = {}

    for local_idx, (global_idx, mol) in enumerate(indexed_mols):
        # Get molecule identity
        try:
            mseq = int(mol.GetProp("mseq")) if mol.HasProp("mseq") else 0
        except (ValueError, KeyError):
            mseq = 0

        mol_info = mol_map.get(mseq, {"name": f"mol_{mseq}", "r_group": "?", "screen_group": "?"})
        name = mol_info["name"]

        # Skip if already processed (resume mode)
        if (name, global_idx) in processed_keys:
            continue

        # Get MOE energy if available
        try:
            moe_energy = float(mol.GetProp("E")) if mol.HasProp("E") else np.nan
        except (ValueError, KeyError):
            moe_energy = np.nan

        # Detect rings (cache per mseq)
        if mseq not in ring_cache:
            try:
                ring1, ring2 = analyzer.identify_aromatic_rings(mol)
                ring_cache[mseq] = (ring1, ring2)
            except Exception as e:
                print(
                    f"  [{local_idx+1}/{n_total}] {name} (idx={global_idx}): "
                    f"ring detection failed: {e}"
                )
                n_errors += 1
                continue
        ring1, ring2 = ring_cache[mseq]

        # Process conformer
        t_conf = time.perf_counter()
        try:
            result = process_conformer(
                mol,
                0,
                analyzer,
                ring1,
                ring2,
                model=args.model,
                device=resolved_device,
                mace_dtype=args.mace_dtype,
                max_steps=args.max_steps,
                fmax=args.fmax,
            )

            # Classify
            classified = analyzer.classify_conformer(
                result["plane_angle_deg"],
                result["interplane_distance_A"],
                result["pi_overlap_pct"],
            )

            # Dark excimer check (Dai et al. 2024 Molecules 29, 507)
            is_dark = (
                classified in ("strong_excimer", "weak_excimer")
                and result["plane_angle_deg"] < 5.0
            )

            record = {
                "name": name,
                "r_group": mol_info["r_group"],
                "screen_group": mol_info["screen_group"],
                "mseq": mseq,
                "conformer_id": global_idx,
                "moe_energy_kcal": moe_energy,
                **result,
                "classification": classified,
                "dark_excimer_warning": is_dark,
                # Distance correction
                "distance_delta_A": (
                    result["moe_interplane_distance_A"]
                    - result["interplane_distance_A"]
                ),
            }

            all_records.append(record)
            n_processed += 1

            elapsed_conf = time.perf_counter() - t_conf
            elapsed_total = time.perf_counter() - t_start
            remaining = (elapsed_total / max(n_processed, 1)) * (
                n_total - local_idx - 1
            )

            if n_processed % 10 == 1 or args.test:
                print(
                    f"  [{local_idx+1}/{n_total}] {name} (idx={global_idx}): "
                    f"d={result['interplane_distance_A']:.2f}A "
                    f"(MOE: {result['moe_interplane_distance_A']:.2f}A, "
                    f"delta={record['distance_delta_A']:+.2f}A) "
                    f"class={classified} "
                    f"({elapsed_conf:.1f}s, ETA {remaining/60:.0f}min)"
                )

        except Exception as e:
            print(f"  [{local_idx+1}/{n_total}] {name} (idx={global_idx}): FAILED - {e}")
            traceback.print_exc()
            n_errors += 1
            continue

        if n_processed % args.save_every == 0:
            df = pd.DataFrame(all_records)
            df.to_csv(output_csv, index=False, encoding="utf-8-sig")
            print(f"  [checkpoint] wrote {len(df)} rows to {output_csv}")

    # Final save
    df = pd.DataFrame(all_records)
    df.to_csv(output_csv, index=False, encoding="utf-8-sig")

    # Generate summary
    elapsed_total = time.perf_counter() - t_start
    print(f"\n{'='*60}")
    print("RE-OPTIMIZATION COMPLETE")
    print(f"{'='*60}")
    print(f"Processed: {n_processed}/{n_total} conformers")
    print(f"Errors: {n_errors}")
    print(f"Time: {elapsed_total:.0f}s ({elapsed_total/60:.1f} min)")
    print(f"Saved: {output_csv}")

    if not df.empty:
        # Per-molecule summary
        summary_records = []
        for mol_name, group in df.groupby("name"):
            n_conf = len(group)
            n_excimer = sum(group["classification"].isin(["strong_excimer", "weak_excimer"]))
            frac = n_excimer / n_conf if n_conf > 0 else 0

            summary_records.append({
                "name": mol_name,
                "r_group": group["r_group"].iloc[0],
                "screen_group": group["screen_group"].iloc[0],
                "n_conformers": n_conf,
                "n_excimer": n_excimer,
                "excimer_fraction": frac,
                "mean_distance_moe": group["moe_interplane_distance_A"].mean(),
                "mean_distance_mace": group["interplane_distance_A"].mean(),
                "mean_distance_delta": group["distance_delta_A"].mean(),
                "mean_angle": group["plane_angle_deg"].mean(),
                "mean_overlap": group["pi_overlap_pct"].mean(),
                "best_overlap": group["pi_overlap_pct"].max(),
                "mean_energy": group["mace_energy_kcal"].mean(),
            })

        summary_df = pd.DataFrame(summary_records)
        summary_df = summary_df.sort_values("excimer_fraction", ascending=False)
        summary_csv = f"{args.output_prefix}_summary.csv"
        summary_df.to_csv(summary_csv, index=False, encoding="utf-8-sig")
        print(f"Saved: {summary_csv}")

        # Print top results
        print(f"\nTop 10 molecules by excimer fraction:")
        print(f"{'Name':<20} {'Frac':>6} {'n_exc':>5} {'d_MOE':>6} {'d_MACE':>7} {'delta':>6}")
        print("-" * 56)
        for _, row in summary_df.head(10).iterrows():
            print(
                f"{row['name']:<20} {row['excimer_fraction']:>5.1%} "
                f"{row['n_excimer']:>5d} {row['mean_distance_moe']:>6.2f} "
                f"{row['mean_distance_mace']:>7.2f} {row['mean_distance_delta']:>+6.2f}"
            )

    # Boltzmann temperature sweep (addresses energy weighting accuracy)
    if not df.empty and "mace_energy_kcal" in df.columns:
        try:
            from pyrene_analyzer.ensemble import compute_boltzmann_weighted_features

            print("\nBoltzmann temperature sweep (MACE energies):")
            for T in [200, 298, 350, 400]:
                try:
                    boltz = compute_boltzmann_weighted_features(
                        df, group_col="name",
                        energy_col="mace_energy_kcal",
                        temperature_K=T,
                    )
                    boltz_csv = f"{args.output_prefix}_boltzmann_{T}K.csv"
                    boltz.to_csv(boltz_csv, encoding="utf-8-sig")
                    # Quick summary: how many molecules have boltz_excimer_fraction > 0?
                    if "boltz_excimer_fraction" in boltz.columns:
                        n_active = (boltz["boltz_excimer_fraction"] > 0).sum()
                        print(f"  {T}K: {n_active}/{len(boltz)} molecules with excimer population")
                except Exception as e:
                    print(f"  {T}K: FAILED - {e}")
        except ImportError:
            print("  (skipping: ensemble module not available)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
