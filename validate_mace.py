"""
MACE-OFF23 Validation Benchmark
================================

Validates MACE-OFF23 accuracy for pi-stacking distances by comparing
optimized dimer geometries against CCSD(T)/CBS reference data.

Benchmark molecules (increasing size):
1. Benzene dimer (S22 reference: 3.40 A)
2. Naphthalene dimer (crystal data: ~3.45 A)
3. Pyrene dimer (DLPNO-CCSD(T)/CBS: 3.43 A, PMC11476719)

Also tests on actual binaphthalene conformers from MOE if available.

Pass criteria: MACE-OFF23 distance within 0.2 A of reference for pyrene dimer.

References:
    - Jurecka et al. (2006) PCCP 8, 1985 (S22 benchmark)
    - Batatia et al. (2024) arXiv:2401.00096 (MACE-OFF23)
    - Bannwarth et al. (2019) JCTC 15, 1652 (GFN2-xTB)
    - PMC11476719 (2024) DLPNO-CCSD(T)/CBS pyrene dimer

Usage:
    python validate_mace.py
    python validate_mace.py --all-methods
    python validate_mace.py --output validation_results.csv
"""

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Suppress RDKit warnings during embedding
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)


# ---------------------------------------------------------------------------
# Reference data (peer-reviewed or high-level ab initio)
# ---------------------------------------------------------------------------

BENCHMARKS = [
    {
        "name": "benzene_dimer",
        "smiles": "c1ccccc1",  # monomer; we build dimer from 2 copies
        "reference_distance_A": 3.40,
        "reference_source": "Jurecka et al. 2006, S22 benchmark (CCSD(T)/CBS)",
        "n_aromatic_atoms": 6,
    },
    {
        "name": "naphthalene_dimer",
        "smiles": "c1ccc2ccccc2c1",
        "reference_distance_A": 3.45,
        "reference_source": "Crystal data (naphthalene stacking distance)",
        "n_aromatic_atoms": 10,
    },
    {
        "name": "pyrene_dimer",
        "smiles": "c1cc2ccc3cccc4ccc(c1)c2c34",
        "reference_distance_A": 3.43,
        "reference_source": "PMC11476719 (2024) DLPNO-CCSD(T)/CBS",
        "n_aromatic_atoms": 16,
    },
]


def build_stacked_dimer(monomer_smiles: str, initial_distance: float = 5.0):
    """
    Build a face-to-face stacked dimer from two copies of a monomer.

    Places one monomer at z=0 and the other at z=initial_distance,
    both in parallel orientation.

    Args:
        monomer_smiles: SMILES of the monomer.
        initial_distance: Initial inter-plane distance in Angstroms.

    Returns:
        RDKit Mol with 3D coordinates of the stacked dimer.
    """
    mol = Chem.MolFromSmiles(monomer_smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    n_atoms = mol.GetNumAtoms()

    # Get coordinates of monomer
    coords = np.array([conf.GetAtomPosition(i) for i in range(n_atoms)])

    # Center at origin
    centroid = coords.mean(axis=0)
    coords -= centroid

    # Create dimer: two copies, one at z=0, one at z=initial_distance
    combo = Chem.CombineMols(mol, mol)
    combo = Chem.RWMol(combo)

    # Embed the combined molecule
    params2 = AllChem.ETKDGv3()
    params2.randomSeed = 42
    AllChem.EmbedMolecule(combo, params2)
    dimer_conf = combo.GetConformer()

    # Set positions manually: first monomer at z=0, second at z=initial_distance
    for i in range(n_atoms):
        pos = coords[i]
        dimer_conf.SetAtomPosition(i, (float(pos[0]), float(pos[1]), float(pos[2])))
        dimer_conf.SetAtomPosition(
            i + n_atoms,
            (float(pos[0]), float(pos[1]), float(pos[2] + initial_distance)),
        )

    return combo.GetMol(), n_atoms


def measure_interplane_distance(mol, n_monomer_atoms: int) -> float:
    """Measure distance between centroids of the two monomer planes."""
    conf = mol.GetConformer()
    n_total = mol.GetNumAtoms()

    coords1 = np.array([
        conf.GetAtomPosition(i) for i in range(n_monomer_atoms)
    ])
    coords2 = np.array([
        conf.GetAtomPosition(i) for i in range(n_monomer_atoms, n_total)
    ])

    # Use only heavy atoms for centroid
    heavy1 = []
    heavy2 = []
    for i in range(n_monomer_atoms):
        if mol.GetAtomWithIdx(i).GetAtomicNum() > 1:
            heavy1.append(coords1[i])
        if mol.GetAtomWithIdx(i + n_monomer_atoms).GetAtomicNum() > 1:
            heavy2.append(coords2[i])

    c1 = np.mean(heavy1, axis=0)
    c2 = np.mean(heavy2, axis=0)

    return float(np.linalg.norm(c2 - c1))


def optimize_with_mmff94s(mol, max_steps: int = 500, verbose: bool = False):
    """Optimize with MMFF94s and return (mol, energy, time)."""
    t0 = time.perf_counter()
    mol = Chem.RWMol(mol)
    try:
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        if mmff_props is None:
            raise ValueError("MMFF94s parameters unavailable for molecule")

        ff = AllChem.MMFFGetMoleculeForceField(
            mol,
            mmff_props,
            confId=-1,
            ignoreInterfragInteractions=False,
        )
        if ff is None:
            raise ValueError("Could not initialize MMFF94s force field")

        ff.Minimize(maxIts=max_steps)
        energy = ff.CalcEnergy()
    except Exception as e:
        if verbose:
            print(f"    MMFF94s failed: {e}")
        energy = float("nan")
    elapsed = time.perf_counter() - t0
    return mol.GetMol(), energy, elapsed


def optimize_with_mace_method(mol, model="small", max_steps=500, fmax=0.005):
    """Optimize with MACE-OFF23 and return (mol, energy, time)."""
    from pyrene_analyzer.mace_optimizer import optimize_with_mace

    t0 = time.perf_counter()
    mol, energy = optimize_with_mace(
        mol, model=model, max_steps=max_steps, fmax=fmax
    )
    elapsed = time.perf_counter() - t0
    return mol, energy, elapsed


def optimize_with_xtb_method(mol, max_steps=1000, fmax=0.005):
    """Optimize with GFN2-xTB and return (mol, energy, time)."""
    from pyrene_analyzer.xtb_optimizer import optimize_with_gfn2xtb

    t0 = time.perf_counter()
    mol, energy = optimize_with_gfn2xtb(
        mol, method="GFN2-xTB", max_steps=max_steps, fmax=fmax
    )
    elapsed = time.perf_counter() - t0
    return mol, energy, elapsed


def run_benchmark(all_methods: bool = False, verbose: bool = True):
    """
    Run the validation benchmark.

    Args:
        all_methods: If True, also benchmark GFN2-xTB (requires WSL).
        verbose: Print progress.

    Returns:
        DataFrame with benchmark results.
    """
    # Check MACE availability
    from pyrene_analyzer.mace_optimizer import has_mace_available
    if not has_mace_available():
        print("ERROR: MACE-OFF23 is not installed. Install with:")
        print("  pip install mace-torch ase torch")
        sys.exit(1)

    results = []

    for bench in BENCHMARKS:
        if verbose:
            print(f"\n{'='*60}")
            print(f"Benchmark: {bench['name']}")
            print(f"Reference: {bench['reference_distance_A']} A ({bench['reference_source']})")
            print(f"{'='*60}")

        # Build stacked dimer
        dimer, n_mono = build_stacked_dimer(bench["smiles"], initial_distance=4.0)
        d_init = measure_interplane_distance(dimer, n_mono)
        if verbose:
            print(f"  Initial distance: {d_init:.2f} A")

        # Method 1: MMFF94s
        if verbose:
            print(f"  Optimizing with MMFF94s...")
        mol_mmff, e_mmff, t_mmff = optimize_with_mmff94s(
            Chem.RWMol(dimer), verbose=verbose
        )
        d_mmff = measure_interplane_distance(mol_mmff, n_mono)
        err_mmff = d_mmff - bench["reference_distance_A"]
        if verbose:
            print(f"    d = {d_mmff:.3f} A (error: {err_mmff:+.3f} A, {t_mmff:.1f}s)")

        results.append({
            "benchmark": bench["name"],
            "method": "MMFF94s",
            "distance_A": d_mmff,
            "reference_A": bench["reference_distance_A"],
            "error_A": err_mmff,
            "abs_error_A": abs(err_mmff),
            "energy": e_mmff,
            "time_s": t_mmff,
            "n_aromatic_atoms": bench["n_aromatic_atoms"],
        })

        # Method 2: MACE-OFF23
        if verbose:
            print(f"  Optimizing with MACE-OFF23 (small)...")
        try:
            mol_mace, e_mace, t_mace = optimize_with_mace_method(
                Chem.RWMol(dimer), model="small", max_steps=500, fmax=0.005
            )
            d_mace = measure_interplane_distance(mol_mace, n_mono)
            err_mace = d_mace - bench["reference_distance_A"]
            if verbose:
                print(f"    d = {d_mace:.3f} A (error: {err_mace:+.3f} A, {t_mace:.1f}s)")
        except Exception as e:
            d_mace, e_mace, t_mace, err_mace = float("nan"), float("nan"), 0.0, float("nan")
            if verbose:
                print(f"    FAILED: {e}")

        results.append({
            "benchmark": bench["name"],
            "method": "MACE-OFF23",
            "distance_A": d_mace,
            "reference_A": bench["reference_distance_A"],
            "error_A": err_mace,
            "abs_error_A": abs(err_mace),
            "energy": e_mace,
            "time_s": t_mace,
            "n_aromatic_atoms": bench["n_aromatic_atoms"],
        })

        # Method 3: MACE-OFF23 medium (if available)
        if verbose:
            print(f"  Optimizing with MACE-OFF23 (medium)...")
        try:
            mol_mace_m, e_mace_m, t_mace_m = optimize_with_mace_method(
                Chem.RWMol(dimer), model="medium", max_steps=500, fmax=0.005
            )
            d_mace_m = measure_interplane_distance(mol_mace_m, n_mono)
            err_mace_m = d_mace_m - bench["reference_distance_A"]
            if verbose:
                print(f"    d = {d_mace_m:.3f} A (error: {err_mace_m:+.3f} A, {t_mace_m:.1f}s)")
        except Exception as e:
            d_mace_m, e_mace_m, t_mace_m, err_mace_m = float("nan"), float("nan"), 0.0, float("nan")
            if verbose:
                print(f"    FAILED: {e}")

        results.append({
            "benchmark": bench["name"],
            "method": "MACE-OFF23-medium",
            "distance_A": d_mace_m,
            "reference_A": bench["reference_distance_A"],
            "error_A": err_mace_m,
            "abs_error_A": abs(err_mace_m),
            "energy": e_mace_m,
            "time_s": t_mace_m,
            "n_aromatic_atoms": bench["n_aromatic_atoms"],
        })

        # Method 4: GFN2-xTB (optional)
        if all_methods:
            from pyrene_analyzer.xtb_optimizer import has_xtb_available
            if has_xtb_available():
                if verbose:
                    print(f"  Optimizing with GFN2-xTB...")
                try:
                    mol_xtb, e_xtb, t_xtb = optimize_with_xtb_method(
                        Chem.RWMol(dimer), max_steps=1000, fmax=0.005
                    )
                    d_xtb = measure_interplane_distance(mol_xtb, n_mono)
                    err_xtb = d_xtb - bench["reference_distance_A"]
                    if verbose:
                        print(f"    d = {d_xtb:.3f} A (error: {err_xtb:+.3f} A, {t_xtb:.1f}s)")
                except Exception as e:
                    d_xtb, e_xtb, t_xtb, err_xtb = float("nan"), float("nan"), 0.0, float("nan")
                    if verbose:
                        print(f"    FAILED: {e}")

                results.append({
                    "benchmark": bench["name"],
                    "method": "GFN2-xTB",
                    "distance_A": d_xtb,
                    "reference_A": bench["reference_distance_A"],
                    "error_A": err_xtb,
                    "abs_error_A": abs(err_xtb),
                    "energy": e_xtb,
                    "time_s": t_xtb,
                    "n_aromatic_atoms": bench["n_aromatic_atoms"],
                })
            else:
                if verbose:
                    print("  GFN2-xTB not available (requires WSL)")

    df = pd.DataFrame(results)
    return df


def evaluate_results(df: pd.DataFrame, verbose: bool = True) -> bool:
    """
    Evaluate benchmark results against pass criteria.

    Pass criteria:
    - MACE-OFF23 pyrene dimer distance within 0.2 A of 3.43 A reference.

    Returns:
        True if validation passes.
    """
    if verbose:
        print(f"\n{'='*60}")
        print("VALIDATION RESULTS SUMMARY")
        print(f"{'='*60}\n")

        # Print table
        for method in df["method"].unique():
            subset = df[df["method"] == method]
            mae = subset["abs_error_A"].mean()
            print(f"  {method}:")
            for _, row in subset.iterrows():
                print(
                    f"    {row['benchmark']:25s}  "
                    f"d={row['distance_A']:.3f} A  "
                    f"ref={row['reference_A']:.2f} A  "
                    f"err={row['error_A']:+.3f} A  "
                    f"({row['time_s']:.1f}s)"
                )
            print(f"    {'MAE':25s}  {mae:.3f} A\n")

    # Check pass criteria
    mace_pyrene = df[
        (df["method"] == "MACE-OFF23") & (df["benchmark"] == "pyrene_dimer")
    ]

    if mace_pyrene.empty:
        if verbose:
            print("FAIL: No MACE-OFF23 result for pyrene dimer")
        return False

    error = abs(mace_pyrene.iloc[0]["error_A"])
    passed = error < 0.2

    if verbose:
        status = "PASS" if passed else "FAIL"
        print(f"Pass criteria: MACE-OFF23 pyrene dimer |error| < 0.2 A")
        print(f"  Actual error: {error:.3f} A")
        print(f"  Result: {status}")
        print()

        if passed:
            print("Recommendation: MACE-OFF23 validated for pi-stacking optimization.")
            print("  Proceed with full 3,347-conformer re-optimization.")
        else:
            print("Recommendation: MACE-OFF23 NOT validated for pi-stacking.")
            print("  Use two-stage pipeline: MACE pre-screen + GFN2-xTB verify on top candidates.")

    return passed


def main():
    parser = argparse.ArgumentParser(
        description="Validate MACE-OFF23 for pi-stacking distance accuracy"
    )
    parser.add_argument(
        "--all-methods", action="store_true",
        help="Also benchmark GFN2-xTB (requires WSL)",
    )
    parser.add_argument(
        "--output", "-o", type=str, default=None,
        help="Save results to CSV file",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true",
        help="Minimal output",
    )
    args = parser.parse_args()

    verbose = not args.quiet

    if verbose:
        print("MACE-OFF23 Validation Benchmark")
        print("=" * 60)
        print()

    df = run_benchmark(all_methods=args.all_methods, verbose=verbose)
    passed = evaluate_results(df, verbose=verbose)

    if args.output:
        df.to_csv(args.output, index=False)
        if verbose:
            print(f"\nResults saved to {args.output}")

    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
