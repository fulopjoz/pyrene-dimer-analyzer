"""
GFN2-xTB Optimization Module
============================

Provides dispersion-corrected geometry optimization using GFN2-xTB
semi-empirical tight-binding method via the ASE interface.

Scientific Basis:
    MMFF94s force field lacks London dispersion forces, resulting in
    systematically overestimated aromatic pi-pi stacking distances
    (3.8-4.2 A vs experimental 3.3-3.5 A). GFN2-xTB includes density-
    dependent dispersion (D4-like) and reproduces pyrene dimer stacking
    distances within +/- 0.1-0.2 A of CCSD(T)/CBS benchmarks.

Performance:
    - GFN2-xTB S22 benchmark: MAE 0.31 kcal/mol
    - GFN2-xTB S66 benchmark: MAE 0.26 kcal/mol
    - Single-point ~2-5 seconds for 140-atom molecule
    - Full optimization ~1-3 minutes for 140-atom molecule

Installation:
    conda install -c conda-forge xtb-python ase

    Note: xtb-python is NOT pip-installable on Windows.
    For CREST functionality, use WSL2.

References:
    - Bannwarth et al. (2019). J. Chem. Theory Comput. 15, 1652-1671.
      DOI: 10.1021/acs.jctc.8b01176
    - PMC11476719 (2024). Pyrene dimer DLPNO-CCSD(T)/CBS benchmark.

Example:
    >>> from pyrene_analyzer.xtb_optimizer import optimize_with_gfn2xtb
    >>> mol, energy = optimize_with_gfn2xtb(mol, conf_id=0, verbose=True)
    >>> print(f"Optimized energy: {energy:.2f} kcal/mol")
"""

import logging
import time
import warnings
from typing import List, Optional, Tuple, Union

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger(__name__)

# Check for xtb-python availability
HAS_XTB = False
XTB_IMPORT_ERROR = None

try:
    from ase import Atoms
    from ase.optimize import BFGS, LBFGS
    from xtb.ase.calculator import XTB
    HAS_XTB = True
except ImportError as e:
    XTB_IMPORT_ERROR = str(e)


# Conversion factors
HARTREE_TO_KCAL = 627.5094740631  # 1 Hartree in kcal/mol
EV_TO_HARTREE = 1.0 / 27.211386245988  # 1 eV in Hartree


def has_xtb_available() -> bool:
    """
    Check if xtb-python is available for use.

    Returns:
        True if xtb-python and ASE are importable, False otherwise.

    Example:
        >>> if has_xtb_available():
        ...     mol = optimize_with_gfn2xtb(mol)
        ... else:
        ...     print("Using MMFF94s fallback")
    """
    return HAS_XTB


def get_xtb_import_error() -> Optional[str]:
    """
    Get the import error message if xtb-python is not available.

    Returns:
        Error message string, or None if xtb is available.
    """
    return XTB_IMPORT_ERROR


def _mol_to_ase_atoms(mol: Chem.Mol, conf_id: int = -1) -> "Atoms":
    """
    Convert RDKit molecule conformer to ASE Atoms object.

    Args:
        mol: RDKit molecule with conformer.
        conf_id: Conformer ID (-1 for first conformer).

    Returns:
        ASE Atoms object with positions and symbols.

    Raises:
        ValueError: If molecule has no conformers.
    """
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no conformers")

    conf = mol.GetConformer(conf_id)
    positions = conf.GetPositions()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    return Atoms(symbols=symbols, positions=positions)


def _update_mol_from_ase(mol: Chem.Mol, atoms: "Atoms", conf_id: int = -1) -> None:
    """
    Update RDKit conformer positions from ASE Atoms object.

    Args:
        mol: RDKit molecule to update.
        atoms: ASE Atoms object with optimized positions.
        conf_id: Conformer ID to update.
    """
    conf = mol.GetConformer(conf_id)
    positions = atoms.get_positions()

    for i, pos in enumerate(positions):
        conf.SetAtomPosition(i, pos)


def optimize_with_gfn2xtb(
    mol: Chem.Mol,
    conf_id: int = -1,
    method: str = "GFN2-xTB",
    max_steps: int = 500,
    fmax: float = 0.05,
    optimizer: str = "BFGS",
    verbose: bool = False,
) -> Tuple[Chem.Mol, float]:
    """
    Optimize a conformer using GFN2-xTB via ASE interface.

    GFN2-xTB provides dispersion-corrected semi-empirical energies,
    correctly reproducing aromatic pi-pi stacking distances (3.3-3.6 A)
    that MMFF94s systematically overestimates (3.8-4.2 A).

    Args:
        mol: RDKit molecule with at least one conformer.
        conf_id: Conformer ID to optimize (-1 for first conformer).
        method: xTB method to use:
            - "GFN2-xTB" (default, recommended)
            - "GFN1-xTB" (faster, less accurate)
            - "GFN0-xTB" (very fast, least accurate)
        max_steps: Maximum optimization steps (default: 500).
        fmax: Force convergence threshold in eV/A (default: 0.05).
            Typical values: 0.05 (standard), 0.01 (tight), 0.1 (loose).
        optimizer: ASE optimizer to use ("BFGS" or "LBFGS").
        verbose: Print optimization progress.

    Returns:
        Tuple of (mol, energy_kcal_mol):
            - mol: RDKit molecule with updated conformer coordinates.
            - energy_kcal_mol: Final energy in kcal/mol.

    Raises:
        ImportError: If xtb-python is not installed.
        ValueError: If molecule has no conformers or method is invalid.
        RuntimeError: If optimization fails to converge.

    Example:
        >>> mol = Chem.MolFromSmiles("c1ccccc1CCc1ccccc1")
        >>> mol = Chem.AddHs(mol)
        >>> AllChem.EmbedMolecule(mol)
        >>> mol, energy = optimize_with_gfn2xtb(mol, verbose=True)
        >>> print(f"Energy: {energy:.2f} kcal/mol")
    """
    if not HAS_XTB:
        raise ImportError(
            f"xtb-python not available: {XTB_IMPORT_ERROR}\n"
            "Install via: conda install -c conda-forge xtb-python ase"
        )

    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no conformers to optimize")

    # Map method names to xTB method strings (keys must be UPPERCASE)
    method_map = {
        "GFN2-XTB": "GFN2-xTB",
        "GFN2": "GFN2-xTB",
        "GFN1-XTB": "GFN1-xTB",
        "GFN1": "GFN1-xTB",
        "GFN0-XTB": "GFN0-xTB",
        "GFN0": "GFN0-xTB",
        "GFNFF": "GFNFF",
        "GFN-FF": "GFNFF",
    }

    xtb_method = method_map.get(method.upper().replace(" ", ""))
    if xtb_method is None:
        raise ValueError(
            f"Unknown method '{method}'. "
            f"Valid options: {list(method_map.keys())}"
        )

    if verbose:
        n_atoms = mol.GetNumAtoms()
        print(f"  Optimizing with {method} ({n_atoms} atoms)...")
        t0 = time.perf_counter()

    # Convert to ASE Atoms
    atoms = _mol_to_ase_atoms(mol, conf_id)

    # Set up xTB calculator
    calc = XTB(method=xtb_method)
    atoms.calc = calc

    # Get initial energy
    try:
        initial_energy_ev = atoms.get_potential_energy()
    except Exception as e:
        raise RuntimeError(f"xTB single-point calculation failed: {e}")

    # Choose optimizer
    if optimizer.upper() == "LBFGS":
        opt = LBFGS(atoms, logfile="-" if verbose else None)
    else:
        opt = BFGS(atoms, logfile="-" if verbose else None)

    # Run optimization
    try:
        converged = opt.run(fmax=fmax, steps=max_steps)
    except Exception as e:
        raise RuntimeError(f"xTB optimization failed: {e}")

    if not converged and verbose:
        warnings.warn(
            f"Optimization did not converge in {max_steps} steps "
            f"(fmax={fmax}). Using final geometry."
        )

    # Get final energy (in eV, convert to kcal/mol)
    final_energy_ev = atoms.get_potential_energy()
    energy_kcal = final_energy_ev * EV_TO_HARTREE * HARTREE_TO_KCAL

    # Update RDKit conformer
    _update_mol_from_ase(mol, atoms, conf_id)

    # Store energy as conformer property
    conf = mol.GetConformer(conf_id)
    conf.SetDoubleProp("energy", energy_kcal)
    conf.SetProp("optimizer", method)

    if verbose:
        elapsed = time.perf_counter() - t0
        delta_e = (final_energy_ev - initial_energy_ev) * EV_TO_HARTREE * HARTREE_TO_KCAL
        print(f"  Done in {elapsed:.1f}s: E = {energy_kcal:.2f} kcal/mol "
              f"(dE = {delta_e:.2f})")

    return mol, energy_kcal


def optimize_conformer_ensemble(
    mol: Chem.Mol,
    method: str = "GFN2-xTB",
    max_steps: int = 500,
    fmax: float = 0.05,
    verbose: bool = False,
    skip_failed: bool = True,
) -> Chem.Mol:
    """
    Optimize all conformers in a molecule with GFN2-xTB.

    Iterates through all conformers and optimizes each independently.
    Failed optimizations can be skipped or raise exceptions.

    Args:
        mol: RDKit molecule with multiple conformers.
        method: xTB method ("GFN2-xTB", "GFN1-xTB", or "GFN0-xTB").
        max_steps: Maximum optimization steps per conformer.
        fmax: Force convergence threshold (eV/A).
        verbose: Print progress for each conformer.
        skip_failed: If True, skip failed optimizations. If False, raise.

    Returns:
        Molecule with optimized conformers. Each conformer has
        "energy" (kcal/mol) and "optimizer" properties set.

    Raises:
        ImportError: If xtb-python is not installed.
        RuntimeError: If any optimization fails and skip_failed=False.

    Example:
        >>> mol = generate_conformers(mol, num_confs=50)
        >>> mol = optimize_conformer_ensemble(mol, verbose=True)
        >>> energies = [conf.GetDoubleProp("energy")
        ...             for conf in mol.GetConformers()]
    """
    if not HAS_XTB:
        raise ImportError(
            f"xtb-python not available: {XTB_IMPORT_ERROR}\n"
            "Install via: conda install -c conda-forge xtb-python ase"
        )

    n_confs = mol.GetNumConformers()
    if n_confs == 0:
        warnings.warn("Molecule has no conformers to optimize")
        return mol

    if verbose:
        print(f"Optimizing {n_confs} conformers with {method}...")
        t_total = time.perf_counter()

    n_success = 0
    n_failed = 0

    for i, conf in enumerate(mol.GetConformers()):
        cid = conf.GetId()

        if verbose:
            print(f"  Conformer {i + 1}/{n_confs} (id={cid})...", end=" ")
            t0 = time.perf_counter()

        try:
            _, energy = optimize_with_gfn2xtb(
                mol,
                conf_id=cid,
                method=method,
                max_steps=max_steps,
                fmax=fmax,
                verbose=False,
            )
            n_success += 1

            if verbose:
                elapsed = time.perf_counter() - t0
                print(f"E = {energy:.2f} kcal/mol ({elapsed:.1f}s)")

        except Exception as e:
            n_failed += 1
            if verbose:
                print(f"FAILED: {e}")

            if not skip_failed:
                raise RuntimeError(
                    f"Optimization failed for conformer {cid}: {e}"
                )

    if verbose:
        total_time = time.perf_counter() - t_total
        print(f"Optimization complete: {n_success}/{n_confs} succeeded "
              f"in {total_time:.1f}s total")

    return mol


def filter_by_energy_xtb(
    mol: Chem.Mol,
    energy_window_kcal: float = 10.0,
) -> Chem.Mol:
    """
    Filter conformers by energy window using xTB energies.

    Removes conformers whose energy exceeds the global minimum
    by more than the specified window. Uses the "energy" property
    set by optimize_with_gfn2xtb or optimize_conformer_ensemble.

    Args:
        mol: Molecule with xTB-optimized conformers.
        energy_window_kcal: Maximum energy above minimum to keep.

    Returns:
        Molecule with high-energy conformers removed.
    """
    energies = []
    for conf in mol.GetConformers():
        try:
            energy = conf.GetDoubleProp("energy")
            energies.append((conf.GetId(), energy))
        except KeyError:
            # No energy property - keep conformer
            energies.append((conf.GetId(), 0.0))

    if not energies:
        return mol

    min_energy = min(e for _, e in energies)
    keep_ids = [
        cid for cid, e in energies
        if (e - min_energy) <= energy_window_kcal
    ]

    if len(keep_ids) == len(energies):
        return mol

    # Rebuild molecule with kept conformers
    new_mol = Chem.RWMol(mol)
    new_mol.RemoveAllConformers()

    for cid in keep_ids:
        conf = mol.GetConformer(cid)
        new_mol.AddConformer(conf, assignId=True)

    return new_mol.GetMol()


def single_point_gfn2xtb(
    mol: Chem.Mol,
    conf_id: int = -1,
    method: str = "GFN2-xTB",
) -> float:
    """
    Calculate single-point energy without optimization.

    Useful for quick energy evaluation or when geometry is
    already optimized by another method.

    Args:
        mol: RDKit molecule with conformer.
        conf_id: Conformer ID (-1 for first).
        method: xTB method to use.

    Returns:
        Energy in kcal/mol.

    Raises:
        ImportError: If xtb-python is not installed.
    """
    if not HAS_XTB:
        raise ImportError(
            f"xtb-python not available: {XTB_IMPORT_ERROR}\n"
            "Install via: conda install -c conda-forge xtb-python ase"
        )

    method_map = {
        "GFN2-XTB": "GFN2-xTB", "GFN2": "GFN2-xTB",
        "GFN1-XTB": "GFN1-xTB", "GFN1": "GFN1-xTB",
        "GFN0-XTB": "GFN0-xTB", "GFN0": "GFN0-xTB",
    }
    xtb_method = method_map.get(method.upper().replace(" ", ""), "GFN2-xTB")

    atoms = _mol_to_ase_atoms(mol, conf_id)
    calc = XTB(method=xtb_method)
    atoms.calc = calc

    energy_ev = atoms.get_potential_energy()
    energy_kcal = energy_ev * EV_TO_HARTREE * HARTREE_TO_KCAL

    return energy_kcal


# Convenience function for integration with existing code
def optimize_or_fallback(
    mol: Chem.Mol,
    conf_id: int = -1,
    preferred_method: str = "GFN2-xTB",
    fallback_method: str = "MMFF94s",
    verbose: bool = False,
) -> Tuple[Chem.Mol, float, str]:
    """
    Optimize with preferred method, falling back if unavailable.

    Tries GFN2-xTB first, falls back to MMFF94s if xtb not installed.
    Useful for code that should work with or without xtb.

    Args:
        mol: RDKit molecule with conformer.
        conf_id: Conformer ID to optimize.
        preferred_method: First choice ("GFN2-xTB").
        fallback_method: Fallback if preferred unavailable ("MMFF94s").
        verbose: Print progress.

    Returns:
        Tuple of (mol, energy_kcal_mol, method_used):
            - mol: Optimized molecule.
            - energy_kcal_mol: Final energy.
            - method_used: Actual method used for optimization.

    Example:
        >>> mol, energy, method = optimize_or_fallback(mol, verbose=True)
        >>> print(f"Optimized with {method}: {energy:.2f} kcal/mol")
    """
    if HAS_XTB and preferred_method.upper().startswith("GFN"):
        try:
            mol, energy = optimize_with_gfn2xtb(
                mol, conf_id=conf_id, method=preferred_method, verbose=verbose
            )
            return mol, energy, preferred_method
        except Exception as e:
            if verbose:
                print(f"  {preferred_method} failed: {e}, trying {fallback_method}")

    # Fallback to MMFF94s
    if fallback_method.upper() == "MMFF94S":
        conf = mol.GetConformer(conf_id)

        if verbose:
            print(f"  Optimizing with MMFF94s...")

        if AllChem.MMFFHasAllMoleculeParams(mol):
            try:
                ff = AllChem.MMFFGetMoleculeForceField(
                    mol, AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s"),
                    confId=conf_id
                )
                ff.Minimize(maxIts=500)
                energy = ff.CalcEnergy()
                conf.SetDoubleProp("energy", energy)
                conf.SetProp("optimizer", "MMFF94s")
                return mol, energy, "MMFF94s"
            except Exception:
                pass

        # UFF fallback
        try:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
            ff.Minimize(maxIts=500)
            energy = ff.CalcEnergy()
            conf.SetDoubleProp("energy", energy)
            conf.SetProp("optimizer", "UFF")
            return mol, energy, "UFF"
        except Exception:
            pass

    # No optimization possible
    return mol, 0.0, "none"
