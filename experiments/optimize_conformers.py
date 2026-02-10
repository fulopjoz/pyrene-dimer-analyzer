"""Parameter sweep for conformer generation optimization.

Tests different ETKDGv3 parameters, constraint levels, and multi-seed strategies
on representative binaphthalene dimer variants to maximize conformer count and
diversity while maintaining geometric accuracy.

Representative variants selected for diversity:
  - EtynPyr_Me:   Top excimer performer (42.9%), 134 atoms, MW 1737
  - DCV_Th_H:     Baseline DCV_Th (0% excimer), 104 atoms, MW 1445
  - CNPh_Th_tBu:  Bulky + excimer anomaly (33.3%), 136 atoms, MW 1861

Usage:
    python experiments/optimize_conformers.py
    python experiments/optimize_conformers.py --variant EtynPyr_Me
    python experiments/optimize_conformers.py --quick  # fast subset only
"""

import argparse
import itertools
import sys
import time
from pathlib import Path

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from rdkit.DistanceGeometry import DoTriangleSmoothing

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pyrene_analyzer.screening import (
    prepare_molecule,
    generate_conformers,
    generate_conformers_biased,
    filter_by_energy,
)
from pyrene_analyzer.core import AromaticDimerAnalyzer


# Representative test variants
TEST_VARIANTS = {
    "EtynPyr_Me": {
        "reason": "Top excimer (42.9%), ethynyl-pyrene + methyl",
        "atoms": 134,
    },
    "DCV_Th_H": {
        "reason": "Baseline DCV_Th (0% excimer), smallest molecule",
        "atoms": 104,
    },
    "CNPh_Th_tBu": {
        "reason": "Bulky + excimer anomaly (33.3%)",
        "atoms": 136,
    },
}


def load_smiles(variant_name: str) -> str:
    """Load SMILES for a variant from binaph_dimer_smiles.csv."""
    csv_path = Path(__file__).resolve().parent.parent / "binaph_dimer_smiles.csv"
    df = pd.read_csv(csv_path)
    row = df[df["name"] == variant_name]
    if row.empty:
        raise ValueError(f"Variant '{variant_name}' not found in {csv_path}")
    return row["smiles"].values[0]


def run_standard_sweep(smiles: str, variant_name: str) -> pd.DataFrame:
    """Sweep standard ETKDGv3 parameters (no biased generation).

    Parameters tested:
      - num_confs: 50, 100, 200, 500
      - pruneRmsThresh: 0.25, 0.5, 0.75, 1.0
      - maxAttempts: default (0), 200, 500, 1000
    """
    mol = Chem.MolFromSmiles(smiles)
    results = []

    num_confs_values = [50, 100, 200, 500]
    prune_rms_values = [0.25, 0.5, 0.75, 1.0]
    max_attempts_values = [0, 200, 500]  # 0 = RDKit default

    # Full sweep would be 48 combinations - do it selectively
    # Fix prune_rms=0.5 and sweep num_confs x maxAttempts
    for num_confs in num_confs_values:
        for max_attempts in max_attempts_values:
            prepped = prepare_molecule(Chem.RWMol(mol), verbose=False)
            if prepped is None:
                continue

            # Ensure Hs
            if not any(a.GetAtomicNum() == 1 for a in prepped.GetAtoms()):
                prepped = Chem.AddHs(prepped)
            prepped.RemoveAllConformers()

            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.numThreads = 0
            params.pruneRmsThresh = 0.5
            params.useRandomCoords = True
            params.enforceChirality = True
            if max_attempts > 0:
                params.maxAttempts = max_attempts

            t0 = time.perf_counter()
            cids = AllChem.EmbedMultipleConfs(
                prepped, numConfs=num_confs, params=params
            )
            embed_time = time.perf_counter() - t0
            n_embedded = len(cids)

            # Optimize
            t0 = time.perf_counter()
            if n_embedded > 0 and AllChem.MMFFHasAllMoleculeParams(prepped):
                try:
                    AllChem.MMFFOptimizeMoleculeConfs(
                        prepped, mmffVariant="MMFF94s", numThreads=0
                    )
                except Exception:
                    pass
            opt_time = time.perf_counter() - t0

            # Energy filter
            n_after_filter = 0
            if n_embedded > 0:
                filtered = filter_by_energy(prepped, 10.0)
                n_after_filter = filtered.GetNumConformers()

            results.append({
                "variant": variant_name,
                "method": "standard",
                "num_confs": num_confs,
                "prune_rms": 0.5,
                "max_attempts": max_attempts,
                "seed": 42,
                "n_embedded": n_embedded,
                "n_after_filter": n_after_filter,
                "embed_time_s": round(embed_time, 1),
                "opt_time_s": round(opt_time, 1),
                "total_time_s": round(embed_time + opt_time, 1),
            })
            print(f"  std num={num_confs} maxAtt={max_attempts}: "
                  f"embedded={n_embedded}, filtered={n_after_filter}, "
                  f"time={embed_time + opt_time:.1f}s")

    # Fix num_confs=200 and sweep prune_rms
    for prune_rms in prune_rms_values:
        if prune_rms == 0.5:
            continue  # already tested above

        prepped = prepare_molecule(Chem.RWMol(mol), verbose=False)
        if prepped is None:
            continue
        if not any(a.GetAtomicNum() == 1 for a in prepped.GetAtoms()):
            prepped = Chem.AddHs(prepped)
        prepped.RemoveAllConformers()

        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.numThreads = 0
        params.pruneRmsThresh = prune_rms
        params.useRandomCoords = True
        params.enforceChirality = True

        t0 = time.perf_counter()
        cids = AllChem.EmbedMultipleConfs(prepped, numConfs=200, params=params)
        embed_time = time.perf_counter() - t0
        n_embedded = len(cids)

        t0 = time.perf_counter()
        if n_embedded > 0 and AllChem.MMFFHasAllMoleculeParams(prepped):
            try:
                AllChem.MMFFOptimizeMoleculeConfs(
                    prepped, mmffVariant="MMFF94s", numThreads=0
                )
            except Exception:
                pass
        opt_time = time.perf_counter() - t0

        n_after_filter = 0
        if n_embedded > 0:
            filtered = filter_by_energy(prepped, 10.0)
            n_after_filter = filtered.GetNumConformers()

        results.append({
            "variant": variant_name,
            "method": "standard",
            "num_confs": 200,
            "prune_rms": prune_rms,
            "max_attempts": 0,
            "seed": 42,
            "n_embedded": n_embedded,
            "n_after_filter": n_after_filter,
            "embed_time_s": round(embed_time, 1),
            "opt_time_s": round(opt_time, 1),
            "total_time_s": round(embed_time + opt_time, 1),
        })
        print(f"  std prune={prune_rms}: "
              f"embedded={n_embedded}, filtered={n_after_filter}, "
              f"time={embed_time + opt_time:.1f}s")

    return pd.DataFrame(results)


def run_biased_sweep(smiles: str, variant_name: str) -> pd.DataFrame:
    """Sweep biased conformer generation parameters.

    Tests:
      - Different constraint level sets
      - Different num_confs values
    """
    mol = Chem.MolFromSmiles(smiles)
    results = []

    constraint_sets = {
        "default": [(3.0, 4.5), (4.0, 6.0), (5.0, 8.0), None],
        "wider": [(3.0, 5.0), (4.0, 7.0), (5.0, 9.0), None],
        "tighter": [(3.0, 4.0), (3.5, 5.0), (4.0, 6.0), None],
        "5_levels": [(3.0, 4.0), (3.5, 5.0), (4.0, 6.0), (5.0, 8.0), None],
        "binaph_wide": [(3.0, 5.5), (4.0, 7.5), (5.5, 10.0), None],
    }

    num_confs_values = [50, 100, 200]

    for cs_name, constraints in constraint_sets.items():
        for num_confs in num_confs_values:
            prepped = prepare_molecule(Chem.RWMol(mol), verbose=False)
            if prepped is None:
                continue

            t0 = time.perf_counter()
            try:
                conf_mol = generate_conformers_biased(
                    prepped,
                    num_confs=num_confs,
                    aromatic_system="binaphthalene",
                    constraint_levels=constraints,
                    prune_rms=0.5,
                    random_seed=42,
                    verbose=False,
                )
                n_generated = conf_mol.GetNumConformers()
                conf_mol = filter_by_energy(conf_mol, 10.0)
                n_after_filter = conf_mol.GetNumConformers()
            except Exception as e:
                print(f"  biased {cs_name} num={num_confs}: FAILED ({e})")
                n_generated = 0
                n_after_filter = 0

            elapsed = time.perf_counter() - t0

            # Analyze excimer fraction if we have conformers
            excimer_frac = 0.0
            if n_after_filter > 0:
                analyzer = AromaticDimerAnalyzer(
                    aromatic_system="binaphthalene", verbose=False
                )
                res_df = analyzer.analyze_molecule(
                    conf_mol, mol_name="test", show_progress=False
                )
                if not res_df.empty:
                    res_df = analyzer.add_classification(res_df)
                    n_exc = sum(
                        res_df["classification"].isin(
                            ["strong_excimer", "weak_excimer"]
                        )
                    )
                    excimer_frac = n_exc / len(res_df)

            results.append({
                "variant": variant_name,
                "method": f"biased_{cs_name}",
                "num_confs": num_confs,
                "constraint_set": cs_name,
                "n_levels": len(constraints),
                "n_generated": n_generated,
                "n_after_filter": n_after_filter,
                "excimer_fraction": round(excimer_frac, 3),
                "total_time_s": round(elapsed, 1),
            })
            print(f"  biased {cs_name} num={num_confs}: "
                  f"gen={n_generated}, filt={n_after_filter}, "
                  f"exc={excimer_frac:.1%}, time={elapsed:.1f}s")

    return pd.DataFrame(results)


def run_multiseed(smiles: str, variant_name: str) -> pd.DataFrame:
    """Test multi-seed strategy: multiple seeds + merge + deduplicate."""
    mol = Chem.MolFromSmiles(smiles)
    results = []

    seed_sets = {
        "1_seed": [42],
        "3_seeds": [42, 137, 512],
        "5_seeds": [42, 137, 256, 512, 1024],
        "10_seeds": [42, 137, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768],
    }

    analyzer = AromaticDimerAnalyzer(
        aromatic_system="binaphthalene", verbose=False
    )

    for ss_name, seeds in seed_sets.items():
        confs_per_seed = max(50 // len(seeds), 10)

        t0 = time.perf_counter()
        all_conf_mols = []

        for seed in seeds:
            prepped = prepare_molecule(Chem.RWMol(mol), verbose=False)
            if prepped is None:
                continue
            try:
                conf_mol = generate_conformers_biased(
                    prepped,
                    num_confs=confs_per_seed,
                    aromatic_system="binaphthalene",
                    prune_rms=0.5,
                    random_seed=seed,
                    verbose=False,
                )
                conf_mol = filter_by_energy(conf_mol, 10.0)
                all_conf_mols.append(conf_mol)
            except Exception:
                pass

        # Count total unique conformers (simple count, no RMSD dedup across seeds)
        total_confs = sum(m.GetNumConformers() for m in all_conf_mols)

        # Analyze each seed's conformers
        all_results = []
        for cm in all_conf_mols:
            if cm.GetNumConformers() > 0:
                res_df = analyzer.analyze_molecule(
                    cm, mol_name="test", show_progress=False
                )
                if not res_df.empty:
                    res_df = analyzer.add_classification(res_df)
                    all_results.append(res_df)

        excimer_frac = 0.0
        n_analyzed = 0
        if all_results:
            combined = pd.concat(all_results, ignore_index=True)
            n_analyzed = len(combined)
            n_exc = sum(
                combined["classification"].isin(
                    ["strong_excimer", "weak_excimer"]
                )
            )
            excimer_frac = n_exc / n_analyzed if n_analyzed > 0 else 0.0

        elapsed = time.perf_counter() - t0

        results.append({
            "variant": variant_name,
            "method": f"multiseed_{ss_name}",
            "n_seeds": len(seeds),
            "confs_per_seed": confs_per_seed,
            "total_conformers": total_confs,
            "n_analyzed": n_analyzed,
            "excimer_fraction": round(excimer_frac, 3),
            "total_time_s": round(elapsed, 1),
        })
        print(f"  {ss_name}: {len(seeds)} seeds x {confs_per_seed} = "
              f"{total_confs} confs, exc={excimer_frac:.1%}, time={elapsed:.1f}s")

    return pd.DataFrame(results)


def run_etkdg_options(smiles: str, variant_name: str) -> pd.DataFrame:
    """Test ETKDGv3 advanced options.

    Tests:
      - useBasicKnowledge
      - useMacrocycleTorsions
      - useSmallRingTorsions
      - ETKDG vs ETKDGv2 vs ETKDGv3
    """
    mol = Chem.MolFromSmiles(smiles)
    results = []

    configs = {
        "ETKDGv3_default": {"version": "v3"},
        "ETKDGv3_macrocycle": {"version": "v3", "useMacrocycleTorsions": True},
        "ETKDGv3_smallring": {"version": "v3", "useSmallRingTorsions": True},
        "ETKDGv3_both_torsions": {
            "version": "v3",
            "useMacrocycleTorsions": True,
            "useSmallRingTorsions": True,
        },
        "ETKDGv2": {"version": "v2"},
        "ETKDG": {"version": "v1"},
    }

    for config_name, config in configs.items():
        prepped = prepare_molecule(Chem.RWMol(mol), verbose=False)
        if prepped is None:
            continue
        if not any(a.GetAtomicNum() == 1 for a in prepped.GetAtoms()):
            prepped = Chem.AddHs(prepped)
        prepped.RemoveAllConformers()

        version = config.pop("version")
        if version == "v3":
            params = AllChem.ETKDGv3()
        elif version == "v2":
            params = AllChem.ETKDGv2()
        else:
            params = AllChem.ETKDG()
        config["version"] = version  # restore

        params.randomSeed = 42
        params.numThreads = 0
        params.pruneRmsThresh = 0.5
        params.useRandomCoords = True
        params.enforceChirality = True

        # Apply extra options
        if config.get("useMacrocycleTorsions"):
            try:
                params.useMacrocycleTorsions = True
            except AttributeError:
                pass
        if config.get("useSmallRingTorsions"):
            try:
                params.useSmallRingTorsions = True
            except AttributeError:
                pass

        t0 = time.perf_counter()
        cids = AllChem.EmbedMultipleConfs(prepped, numConfs=100, params=params)
        embed_time = time.perf_counter() - t0
        n_embedded = len(cids)

        # Optimize
        t0 = time.perf_counter()
        if n_embedded > 0 and AllChem.MMFFHasAllMoleculeParams(prepped):
            try:
                AllChem.MMFFOptimizeMoleculeConfs(
                    prepped, mmffVariant="MMFF94s", numThreads=0
                )
            except Exception:
                pass
        opt_time = time.perf_counter() - t0

        n_after_filter = 0
        if n_embedded > 0:
            filtered = filter_by_energy(prepped, 10.0)
            n_after_filter = filtered.GetNumConformers()

        results.append({
            "variant": variant_name,
            "method": config_name,
            "version": version,
            "n_embedded": n_embedded,
            "n_after_filter": n_after_filter,
            "embed_time_s": round(embed_time, 1),
            "opt_time_s": round(opt_time, 1),
            "total_time_s": round(embed_time + opt_time, 1),
        })
        print(f"  {config_name}: embedded={n_embedded}, "
              f"filtered={n_after_filter}, time={embed_time + opt_time:.1f}s")

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description="Parameter sweep for conformer generation optimization"
    )
    parser.add_argument(
        "--variant",
        type=str,
        default=None,
        help="Test a single variant (default: all 3 representatives)",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run quick subset (ETKDG options + multiseed only)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="experiments/results",
        help="Output directory for results CSVs",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine which variants to test
    if args.variant:
        if args.variant not in TEST_VARIANTS:
            print(f"Unknown variant '{args.variant}'. "
                  f"Available: {', '.join(TEST_VARIANTS.keys())}")
            sys.exit(1)
        variants = {args.variant: TEST_VARIANTS[args.variant]}
    else:
        variants = TEST_VARIANTS

    print("=" * 70)
    print("CONFORMER GENERATION PARAMETER SWEEP")
    print("=" * 70)
    print(f"Variants: {', '.join(variants.keys())}")
    print(f"Quick mode: {args.quick}")
    print(f"Output: {output_dir}")
    print()

    all_results = []

    for variant_name, info in variants.items():
        print(f"\n{'=' * 60}")
        print(f"VARIANT: {variant_name}")
        print(f"  {info['reason']}")
        print(f"  {info['atoms']} atoms")
        print("=" * 60)

        smiles = load_smiles(variant_name)

        # 1. ETKDG version/options sweep (fast)
        print("\n--- ETKDG Version & Options ---")
        etkdg_df = run_etkdg_options(smiles, variant_name)
        all_results.append(etkdg_df)

        # 2. Multi-seed strategy (moderate)
        print("\n--- Multi-Seed Strategy ---")
        multiseed_df = run_multiseed(smiles, variant_name)
        all_results.append(multiseed_df)

        if not args.quick:
            # 3. Standard parameter sweep (slow)
            print("\n--- Standard Parameter Sweep ---")
            standard_df = run_standard_sweep(smiles, variant_name)
            all_results.append(standard_df)

            # 4. Biased constraint sweep (slow)
            print("\n--- Biased Constraint Sweep ---")
            biased_df = run_biased_sweep(smiles, variant_name)
            all_results.append(biased_df)

    # Save results
    print("\n" + "=" * 70)
    print("SAVING RESULTS")
    print("=" * 70)

    for df in all_results:
        if not df.empty:
            # Determine filename from method column
            method = df["method"].iloc[0].split("_")[0] if "method" in df.columns else "unknown"
            variant = df["variant"].iloc[0] if "variant" in df.columns else "unknown"

    # Combine all results
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        outpath = output_dir / "conformer_optimization_sweep.csv"
        combined.to_csv(outpath, index=False)
        print(f"Saved {len(combined)} results to {outpath}")

        # Print summary of best configurations
        print("\n" + "=" * 70)
        print("BEST CONFIGURATIONS BY CONFORMER COUNT")
        print("=" * 70)

        for variant_name in variants:
            v_df = combined[combined["variant"] == variant_name]
            if v_df.empty:
                continue
            print(f"\n{variant_name}:")

            # Find column with conformer count
            count_col = None
            for col in ["n_after_filter", "n_analyzed", "total_conformers"]:
                if col in v_df.columns:
                    count_col = col
                    break

            if count_col:
                top5 = v_df.nlargest(5, count_col)
                for _, row in top5.iterrows():
                    method = row.get("method", "?")
                    n = row.get(count_col, 0)
                    t = row.get("total_time_s", 0)
                    exc = row.get("excimer_fraction", "N/A")
                    print(f"  {method:<30s}: {n:>4} conformers, "
                          f"exc={exc}, time={t:.0f}s")


if __name__ == "__main__":
    main()
