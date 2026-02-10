"""
MACE-OFF23 Optimization Module
===============================

Provides dispersion-corrected geometry optimization using MACE-OFF23
machine learning interatomic potential via the ASE interface.

Scientific Basis:
    MACE-OFF23 is a pre-trained equivariant neural network potential
    trained on ~1.1M conformations at the wB97M-D3(BJ)/def2-TZVPPD level,
    which is near gold-standard for non-covalent interactions. It correctly
    reproduces aromatic pi-pi stacking distances (3.3-3.5 A) that classical
    force fields (MMFF94s, Amber:EHT) systematically overestimate.

Performance:
    - 10-20x faster than GFN2-xTB for 130-atom molecules
    - MAE ~0.3 kcal/mol on biaryl torsion barriers
    - Supports: H, B, C, N, O, F, Si, P, S, Cl, Br, I
    - Works on CPU (pip install) or GPU (CUDA)

Installation:
    pip install mace-torch ase torch

References:
    - Batatia et al. (2024). arXiv:2401.00096 (MACE-OFF foundation model)
    - Batatia et al. (2022). NeurIPS (MACE architecture)
    - Training data: wB97M-D3(BJ)/def2-TZVPPD
      Najibi & Goerigk (2020). JCTC 16, 4479.

Example:
    >>> from pyrene_analyzer.mace_optimizer import optimize_with_mace
    >>> mol, energy = optimize_with_mace(mol, conf_id=0, verbose=True)
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

# Check for MACE availability
HAS_MACE = False
MACE_IMPORT_ERROR = None

try:
    from ase import Atoms
    from ase.optimize import BFGS, LBFGS
    from mace.calculators import mace_off
    HAS_MACE = True
except ImportError as e:
    MACE_IMPORT_ERROR = str(e)


# Conversion factors
EV_TO_KCAL = 23.060541945329334  # 1 eV in kcal/mol


def has_mace_available() -> bool:
    """
    Check if MACE-OFF23 is available for use.

    Returns:
        True if mace-torch and ASE are importable, False otherwise.

    Example:
        >>> if has_mace_available():
        ...     mol = optimize_with_mace(mol)
        ... else:
        ...     print("Using MMFF94s fallback")
    """
    return HAS_MACE


def get_mace_import_error() -> Optional[str]:
    """
    Get the import error message if MACE is not available.

    Returns:
        Error message string, or None if MACE is available.
    """
    return MACE_IMPORT_ERROR


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


# Cache for MACE calculator (loading the model is expensive)
_MACE_CALC_CACHE = {}


def _get_mace_calculator(
    model: str = "medium",
    device: str = "cpu",
    default_dtype: str = "float64",
):
    """Get or create a cached MACE calculator."""
    # Preserve the historical 2-tuple cache key for float64 to avoid
    # breaking external checks while allowing alternate dtypes.
    cache_key = (model, device) if default_dtype == "float64" else (
        model,
        device,
        default_dtype,
    )
    if cache_key not in _MACE_CALC_CACHE:
        _MACE_CALC_CACHE[cache_key] = mace_off(
            model=model, device=device, default_dtype=default_dtype
        )
    return _MACE_CALC_CACHE[cache_key]


def optimize_with_mace(
    mol: Chem.Mol,
    conf_id: int = -1,
    model: str = "medium",
    device: str = "cpu",
    default_dtype: str = "float64",
    max_steps: int = 500,
    fmax: float = 0.05,
    optimizer: str = "LBFGS",
    verbose: bool = False,
) -> Tuple[Chem.Mol, float]:
    """
    Optimize a conformer using MACE-OFF23 via ASE interface.

    MACE-OFF23 provides DFT-quality energies (trained on wB97M-D3BJ)
    at a fraction of the cost. Correctly reproduces aromatic pi-pi
    stacking distances (3.3-3.5 A).

    Args:
        mol: RDKit molecule with at least one conformer.
        conf_id: Conformer ID to optimize (-1 for first conformer).
        model: MACE-OFF23 model size:
            - "small" (fastest, ~1M params)
            - "medium" (default, balanced, ~10M params)
            - "large" (most accurate, ~30M params)
        device: Compute device ("cpu" or "cuda").
        default_dtype: MACE dtype ("float64" for highest accuracy,
            "float32" for faster GPU throughput).
        max_steps: Maximum optimization steps (default: 500).
        fmax: Force convergence threshold in eV/A (default: 0.05).
        optimizer: ASE optimizer ("LBFGS" recommended, or "BFGS").
        verbose: Print optimization progress.

    Returns:
        Tuple of (mol, energy_kcal_mol):
            - mol: RDKit molecule with updated conformer coordinates.
            - energy_kcal_mol: Final energy in kcal/mol.

    Raises:
        ImportError: If mace-torch is not installed.
        ValueError: If molecule has no conformers.
        RuntimeError: If optimization fails.

    References:
        Batatia et al. (2024). arXiv:2401.00096
        Najibi & Goerigk (2020). JCTC 16, 4479.

    Example:
        >>> mol, energy = optimize_with_mace(mol, model="medium", verbose=True)
        >>> print(f"Energy: {energy:.2f} kcal/mol")
    """
    if not HAS_MACE:
        raise ImportError(
            f"MACE-OFF23 not available: {MACE_IMPORT_ERROR}\n"
            "Install via: pip install mace-torch ase torch"
        )

    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no conformers to optimize")

    if verbose:
        n_atoms = mol.GetNumAtoms()
        print(f"  Optimizing with MACE-OFF23-{model} ({n_atoms} atoms)...")
        t0 = time.perf_counter()

    # Convert to ASE Atoms
    atoms = _mol_to_ase_atoms(mol, conf_id)

    # Set up MACE calculator (cached)
    calc = _get_mace_calculator(
        model=model,
        device=device,
        default_dtype=default_dtype,
    )
    atoms.calc = calc

    # Get initial energy
    try:
        initial_energy_ev = atoms.get_potential_energy()
    except Exception as e:
        raise RuntimeError(f"MACE single-point calculation failed: {e}")

    # Choose optimizer (LBFGS recommended for MACE)
    if optimizer.upper() == "BFGS":
        opt = BFGS(atoms, logfile="-" if verbose else None)
    else:
        opt = LBFGS(atoms, logfile="-" if verbose else None)

    # Run optimization
    try:
        converged = opt.run(fmax=fmax, steps=max_steps)
    except Exception as e:
        raise RuntimeError(f"MACE optimization failed: {e}")

    if not converged and verbose:
        warnings.warn(
            f"Optimization did not converge in {max_steps} steps "
            f"(fmax={fmax}). Using final geometry."
        )

    # Get final energy (eV → kcal/mol)
    final_energy_ev = atoms.get_potential_energy()
    energy_kcal = final_energy_ev * EV_TO_KCAL

    # Update RDKit conformer
    _update_mol_from_ase(mol, atoms, conf_id)

    # Store energy as conformer property
    conf = mol.GetConformer(conf_id)
    conf.SetDoubleProp("energy", energy_kcal)
    conf.SetProp("optimizer", f"MACE-OFF23-{model}")

    if verbose:
        elapsed = time.perf_counter() - t0
        delta_e = (final_energy_ev - initial_energy_ev) * EV_TO_KCAL
        print(f"  Done in {elapsed:.1f}s: E = {energy_kcal:.2f} kcal/mol "
              f"(dE = {delta_e:.2f})")

    return mol, energy_kcal


def optimize_conformer_ensemble_mace(
    mol: Chem.Mol,
    model: str = "medium",
    device: str = "cpu",
    default_dtype: str = "float64",
    max_steps: int = 500,
    fmax: float = 0.05,
    verbose: bool = False,
    skip_failed: bool = True,
) -> Chem.Mol:
    """
    Optimize all conformers in a molecule with MACE-OFF23.

    Iterates through all conformers and optimizes each independently.
    The MACE calculator is loaded once and reused for all conformers.

    Args:
        mol: RDKit molecule with multiple conformers.
        model: MACE-OFF23 model size ("small", "medium", "large").
        device: Compute device ("cpu" or "cuda").
        default_dtype: MACE dtype ("float64" or "float32").
        max_steps: Maximum optimization steps per conformer.
        fmax: Force convergence threshold (eV/A).
        verbose: Print progress for each conformer.
        skip_failed: If True, skip failed optimizations. If False, raise.

    Returns:
        Molecule with optimized conformers. Each conformer has
        "energy" (kcal/mol) and "optimizer" properties set.

    Raises:
        ImportError: If MACE is not installed.
        RuntimeError: If any optimization fails and skip_failed=False.

    Example:
        >>> mol = optimize_conformer_ensemble_mace(mol, verbose=True)
        >>> energies = [conf.GetDoubleProp("energy")
        ...             for conf in mol.GetConformers()]
    """
    if not HAS_MACE:
        raise ImportError(
            f"MACE-OFF23 not available: {MACE_IMPORT_ERROR}\n"
            "Install via: pip install mace-torch ase torch"
        )

    n_confs = mol.GetNumConformers()
    if n_confs == 0:
        warnings.warn("Molecule has no conformers to optimize")
        return mol

    if verbose:
        print(f"Optimizing {n_confs} conformers with MACE-OFF23-{model}...")
        t_total = time.perf_counter()

    n_success = 0
    n_failed = 0

    for i, conf in enumerate(mol.GetConformers()):
        cid = conf.GetId()

        if verbose:
            print(f"  Conformer {i + 1}/{n_confs} (id={cid})...", end=" ")
            t0 = time.perf_counter()

        try:
            _, energy = optimize_with_mace(
                mol,
                conf_id=cid,
                model=model,
                device=device,
                default_dtype=default_dtype,
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


def filter_by_energy_mace(
    mol: Chem.Mol,
    energy_window_kcal: float = 10.0,
) -> Chem.Mol:
    """
    Filter conformers by energy window using MACE energies.

    Removes conformers whose energy exceeds the global minimum
    by more than the specified window. Uses the "energy" property
    set by optimize_with_mace or optimize_conformer_ensemble_mace.

    Args:
        mol: Molecule with MACE-optimized conformers.
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

    new_mol = Chem.RWMol(mol)
    new_mol.RemoveAllConformers()

    for cid in keep_ids:
        conf = mol.GetConformer(cid)
        new_mol.AddConformer(conf, assignId=True)

    return new_mol.GetMol()


def single_point_mace(
    mol: Chem.Mol,
    conf_id: int = -1,
    model: str = "medium",
    device: str = "cpu",
    default_dtype: str = "float64",
) -> float:
    """
    Calculate single-point energy without optimization.

    Args:
        mol: RDKit molecule with conformer.
        conf_id: Conformer ID (-1 for first).
        model: MACE-OFF23 model size.
        device: Compute device.
        default_dtype: MACE dtype ("float64" or "float32").

    Returns:
        Energy in kcal/mol.

    Raises:
        ImportError: If MACE is not installed.
    """
    if not HAS_MACE:
        raise ImportError(
            f"MACE-OFF23 not available: {MACE_IMPORT_ERROR}\n"
            "Install via: pip install mace-torch ase torch"
        )

    atoms = _mol_to_ase_atoms(mol, conf_id)
    calc = _get_mace_calculator(
        model=model,
        device=device,
        default_dtype=default_dtype,
    )
    atoms.calc = calc

    energy_ev = atoms.get_potential_energy()
    energy_kcal = energy_ev * EV_TO_KCAL

    return energy_kcal


def optimize_or_fallback_mace(
    mol: Chem.Mol,
    conf_id: int = -1,
    preferred_model: str = "medium",
    device: str = "cpu",
    fallback_method: str = "MMFF94s",
    verbose: bool = False,
) -> Tuple[Chem.Mol, float, str]:
    """
    Optimize with MACE-OFF23, falling back if unavailable.

    Tries MACE-OFF23 first, falls back to MMFF94s/UFF if not installed.

    Args:
        mol: RDKit molecule with conformer.
        conf_id: Conformer ID to optimize.
        preferred_model: MACE model size ("small", "medium", "large").
        device: Compute device for MACE.
        fallback_method: Fallback if MACE unavailable ("MMFF94s").
        verbose: Print progress.

    Returns:
        Tuple of (mol, energy_kcal_mol, method_used).

    Example:
        >>> mol, energy, method = optimize_or_fallback_mace(mol, verbose=True)
        >>> print(f"Optimized with {method}: {energy:.2f} kcal/mol")
    """
    if HAS_MACE:
        try:
            mol, energy = optimize_with_mace(
                mol, conf_id=conf_id, model=preferred_model,
                device=device, verbose=verbose
            )
            return mol, energy, f"MACE-OFF23-{preferred_model}"
        except Exception as e:
            if verbose:
                print(f"  MACE-OFF23 failed: {e}, trying {fallback_method}")

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

    return mol, 0.0, "none"
