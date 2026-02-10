"""
Virtual Screening Module
========================

Substituent screening pipeline for aromatic dimer conformer analysis.

This module provides tools for virtual screening by R-group enumeration,
conformer generation (ETKDGv3 + MMFF94s), and geometric analysis using
the existing AromaticDimerAnalyzer.

Scientific Basis:
    Geometric parameters (θ, d, π-overlap) are proven proxies for excimer
    formation potential. Causality established through:
    - Ge et al. (2020): π-overlap > distance for excimer prediction
    - Basuroy et al. (2021): 42% overlap = excimer, 17% = monomer
    - Fluorene steric control (2018): geometry change → photophysics change

Limitations:
    - Static analysis only (no thermal dynamics or kinetic pathways)
    - No quantitative IE/IM ratio prediction
    - MMFF94s energies not accurate enough for Boltzmann weighting
    - R-group replacement assumes symmetric substitution on both aromatic units

Example:
    >>> from pyrene_analyzer.screening import SubstituentScreener
    >>> screener = SubstituentScreener(
    ...     template_sdf="dimer.sdf",
    ...     r_group_smarts="[CH2][CH3]",  # ethyl groups to replace
    ... )
    >>> results, summary = screener.screen(num_confs=50)
    >>> print(summary.sort_values("excimer_fraction", ascending=False))

References:
    - Wang, Witek, Landrum, Riniker (2020). J. Chem. Inf. Model. [ETKDGv3]
    - Tosco, Stiefl, Landrum (2014). J. Cheminform. [RDKit MMFF94]
    - Ge et al. (2020). J. Mater. Chem. C, 8, 10223
    - Basuroy et al. (2021). J. Chem. Phys., 155, 234304
"""

import json
import logging
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.DistanceGeometry import DoTriangleSmoothing

from pyrene_analyzer.core import AromaticDimerAnalyzer
from pyrene_analyzer.xtb_optimizer import (
    has_xtb_available,
    optimize_conformer_ensemble,
    filter_by_energy_xtb,
)
from pyrene_analyzer.mace_optimizer import (
    has_mace_available,
    optimize_conformer_ensemble_mace,
    filter_by_energy_mace,
)

logger = logging.getLogger(__name__)

# Valid optimizer choices
VALID_OPTIMIZERS = ("MMFF94s", "GFN2-xTB", "MACE-OFF23", "none")


def _estimate_pipeline_time(
    n_heavy_atoms: int,
    num_confs: int,
    use_biased: bool,
    optimizer: str = "MMFF94s",
) -> str:
    """Rough timing estimate for user feedback. Not meant to be precise."""
    base = n_heavy_atoms * 0.05
    embed_time = base * num_confs * 0.1

    opt_upper = optimizer.upper()
    if opt_upper == "MACE-OFF23":
        # MACE: ~50% of GFN2-xTB, 5-30x slower than MMFF94s
        opt_time = base * num_confs * 1.5
    elif opt_upper == "GFN2-XTB":
        # GFN2-xTB: ~200s per conformer for 130-atom dimers
        opt_time = base * num_confs * 3.0
    else:
        # MMFF94s or none
        opt_time = base * num_confs * 0.3

    if use_biased:
        embed_time *= 1.5
    total = embed_time + opt_time + 5
    low = max(10, int(total * 0.7))
    high = int(total * 1.8)
    return f"~{low}-{high}s"


# ---------------------------------------------------------------------------
# Built-in substituent library (~25 groups)
# Covers steric (H→tBu), electronic (NH2→NO2), aromatic (Ph), LC-relevant
# SMILES: first atom is the attachment point (bonds to scaffold)
# ---------------------------------------------------------------------------

COMMON_SUBSTITUENTS: Dict[str, str] = {
    # --- Alkyl (steric series) ---
    "H": "[H]",
    "Me": "C",
    "Et": "CC",
    "nPr": "CCC",
    "iPr": "C(C)C",
    "nBu": "CCCC",
    "tBu": "C(C)(C)C",
    "cHex": "C1CCCCC1",
    "cPr": "C1CC1",
    # --- Electron-donating ---
    "OMe": "OC",
    "OEt": "OCC",
    "NH2": "N",
    "NMe2": "N(C)C",
    "OH": "O",
    # --- Electron-withdrawing ---
    "F": "F",
    "Cl": "Cl",
    "Br": "Br",
    "CF3": "C(F)(F)F",
    "NO2": "[N+](=O)[O-]",
    "CN": "C#N",
    "COMe": "C(=O)C",
    # --- Aromatic ---
    "Ph": "c1ccccc1",
    "Bn": "Cc1ccccc1",
    # --- LC-relevant (mesogenic tails) ---
    "OC6H13": "OCCCCCC",
}


def validate_substituent_library(
    substituents: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """
    Validate that all substituent SMILES parse correctly.

    Args:
        substituents: Dict of name→SMILES. If None, validates COMMON_SUBSTITUENTS.

    Returns:
        Dict of valid name→SMILES pairs (invalid ones are dropped with warnings).
    """
    if substituents is None:
        substituents = COMMON_SUBSTITUENTS

    valid = {}
    for name, smiles in substituents.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            warnings.warn(f"Invalid substituent SMILES for '{name}': {smiles}")
        else:
            valid[name] = smiles
    return valid


def prepare_molecule(
    mol: Chem.Mol, verbose: bool = False
) -> Optional[Chem.Mol]:
    """
    Prepare a molecule for conformer generation (open-source MOE Wash equivalent).

    Steps:
        1. Sanitize and normalize functional groups
        2. Remove salts / keep largest fragment
        3. Add explicit hydrogens (required before embedding)
        4. Generate initial 3D coordinates if missing
        5. Energy-minimize with MMFF94s (MMFF94x equivalent)

    Protonation at specific pH is skipped because pyrene dimers with
    alkyl/aryl substituents have no ionizable groups. For ionizable
    substituents (e.g., -COOH, -NH2 at specific pH), Dimorphite-DL
    can be added as a future enhancement.

    Args:
        mol: RDKit molecule (may lack 3D coordinates).
        verbose: If True, print progress messages with timing.

    Returns:
        Prepared molecule with 3D coordinates, or None if preparation fails.
    """
    if mol is None:
        return None

    if verbose:
        print("  Sanitizing and standardizing...")
    try:
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol)
    except Exception:
        pass  # proceed with original mol if standardization fails

    mol = Chem.AddHs(mol)

    if mol.GetNumConformers() == 0:
        if verbose:
            print("  Generating initial 3D coordinates...")
            t0 = time.perf_counter()
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.useRandomCoords = True
        status = AllChem.EmbedMolecule(mol, params)
        if status == -1:
            warnings.warn("Failed to generate initial 3D coordinates")
            return None
        if verbose:
            print(f"  3D coordinates generated in {time.perf_counter() - t0:.1f}s")

    # MMFF94s minimize with UFF fallback
    if verbose:
        print("  Optimizing geometry with MMFF94s...")
        t0 = time.perf_counter()
    if AllChem.MMFFHasAllMoleculeParams(mol):
        try:
            result = AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s", maxIters=500)
            if result == -1:
                if verbose:
                    print("  MMFF94s did not converge, falling back to UFF...")
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            try:
                if verbose:
                    print("  MMFF94s failed, falling back to UFF...")
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                if verbose:
                    print("  Force field optimization failed, continuing with unoptimized geometry")
    else:
        if verbose:
            print("  MMFF94s cannot parameterize molecule, trying UFF...")
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            if verbose:
                print("  UFF also failed, continuing with unoptimized geometry")
    if verbose:
        print(f"  Optimization done in {time.perf_counter() - t0:.1f}s")

    return mol


def replace_r_groups(
    template_mol: Chem.Mol,
    r_group_smarts: str,
    new_substituent_smiles: str,
    replace_all: bool = True,
) -> Optional[Chem.Mol]:
    """
    Replace R-groups on a template molecule with a new substituent.

    Uses RDKit ReplaceSubstructs: the first atom of the replacement SMILES
    bonds to wherever the first atom of the SMARTS pattern was connected
    to the rest of the molecule.

    Args:
        template_mol: Template molecule with existing R-groups.
        r_group_smarts: SMARTS pattern matching the R-group to replace.
        new_substituent_smiles: SMILES of the new substituent. First atom
            will be the attachment point.
        replace_all: If True, replace all matches (needed for symmetric
            dimers with identical R-groups on both aromatic units).

    Returns:
        New molecule with R-groups replaced, or None if replacement fails.

    Note:
        For symmetric dimers, replaceAll=True is critical. ReplaceSubstructs
        replaces one match at a time when replaceAll=False, returning a tuple
        of molecules (one per match). With replaceAll=True, all matches are
        replaced in a single molecule.
    """
    pattern = Chem.MolFromSmarts(r_group_smarts)
    if pattern is None:
        warnings.warn(f"Invalid R-group SMARTS: {r_group_smarts}")
        return None

    replacement = Chem.MolFromSmiles(new_substituent_smiles)
    if replacement is None:
        warnings.warn(f"Invalid substituent SMILES: {new_substituent_smiles}")
        return None

    matches = list(template_mol.GetSubstructMatches(pattern))
    if not matches:
        return Chem.Mol(template_mol)

    match_sets = [set(m) for m in matches]
    if replace_all:
        for i in range(len(match_sets)):
            for j in range(i + 1, len(match_sets)):
                if match_sets[i] & match_sets[j]:
                    warnings.warn(
                        "Overlapping R-group matches detected; replacement skipped"
                    )
                    return None

    for match in matches:
        attachment_points = 0
        for idx in match:
            atom = template_mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1:
                    continue
                if neighbor.GetIdx() not in match:
                    attachment_points += 1
        if attachment_points != 1:
            warnings.warn(
                "R-group SMARTS has multiple attachment points; replacement skipped"
            )
            return None

    try:
        replaced = AllChem.ReplaceSubstructs(
            template_mol, pattern, replacement, replaceAll=replace_all
        )
    except Exception as e:
        warnings.warn(f"ReplaceSubstructs failed: {e}")
        return None

    if not replaced:
        return None

    # replaceAll=True returns a single-element tuple
    new_mol = replaced[0]

    try:
        sanitize_status = Chem.SanitizeMol(
            new_mol, catchErrors=True
        )
        if sanitize_status != Chem.SanitizeFlags.SANITIZE_NONE:
            # Expected for incompatible substituent/template combinations
            # (e.g., valence violations when F replaces a multi-bonded group).
            # Caller handles None returns gracefully.
            return None
        return new_mol
    except Exception:
        # Sanitization exception — incompatible replacement
        return None


def generate_conformers(
    mol: Chem.Mol,
    num_confs: int = 50,
    prune_rms: float = 0.5,
    random_seed: int = 42,
    optimize: bool = True,
    verbose: bool = False,
    optimizer: str = "MMFF94s",
) -> Chem.Mol:
    """
    Generate a conformer ensemble using RDKit ETKDGv3 + MMFF94s.

    Best practices from Greg Landrum, Pat Walters, and iwatobipen:
    - ETKDGv3 default since RDKit 2024.03 (Wang, Witek, Landrum, Riniker 2020)
    - Always add Hs before embedding (Landrum)
    - useRandomCoords=True for robustness (Landrum 2021)
    - pruneRmsThresh=0.5 for deduplication (Walters pattern)
    - MMFF94s optimization after embedding (Tosco et al. 2014)
    - numThreads=0 for parallel embedding

    Args:
        mol: RDKit molecule (Hs should already be added).
        num_confs: Number of conformers to generate.
        prune_rms: RMSD threshold for pruning near-duplicate conformers.
        random_seed: Random seed for reproducibility.
        optimize: If True, minimize each conformer with the chosen optimizer.
        verbose: If True, print progress messages with timing.
        optimizer: Optimization method ("MMFF94s", "GFN2-xTB", "MACE-OFF23", or "none").
            - "MMFF94s": Fast but lacks London dispersion (pi-stacking too long)
            - "GFN2-xTB": Dispersion-corrected, correct pi-stacking
            - "MACE-OFF23": DFT-quality dispersion, 10-20x faster than xTB (recommended)
            - "none": Skip optimization (fastest, but less accurate)
            Default: "MMFF94s" for backward compatibility.

    Returns:
        Molecule with embedded (and optionally optimized) conformers.
        Energies stored as conformer property "energy" (kcal/mol).
    """
    # Ensure Hs are present
    if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        mol = Chem.AddHs(mol)

    # Remove any existing conformers
    mol.RemoveAllConformers()

    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.numThreads = 0
    params.pruneRmsThresh = prune_rms
    params.useRandomCoords = True
    params.enforceChirality = True

    if verbose:
        print(f"  Embedding {num_confs} conformers (ETKDGv3)...")
        t0 = time.perf_counter()

    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

    if len(cids) == 0:
        # Retry with more attempts for difficult molecules
        if verbose:
            print("  Retrying with increased max attempts...")
        params.maxIterations = 1000
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

    if verbose:
        print(f"  Embedded {len(cids)} conformers in {time.perf_counter() - t0:.1f}s")

    if len(cids) == 0:
        warnings.warn("Could not generate any conformers")
        return mol

    # Optimization
    if optimize and len(cids) > 0 and optimizer.lower() != "none":
        opt_method = optimizer.upper()

        # Try GFN2-xTB if requested and available
        if opt_method == "GFN2-XTB":
            if has_xtb_available():
                if verbose:
                    print(f"  Optimizing {len(cids)} conformers with GFN2-xTB...")
                    t0 = time.perf_counter()
                try:
                    mol = optimize_conformer_ensemble(mol, method="GFN2-xTB", verbose=verbose)
                    if verbose:
                        print(f"  Optimization complete in {time.perf_counter() - t0:.1f}s")
                except Exception as e:
                    if verbose:
                        print(f"  GFN2-xTB failed ({e}), falling back to MMFF94s...")
                    opt_method = "MMFF94S"  # Fall back to MMFF94s
            else:
                if verbose:
                    print("  GFN2-xTB not available, falling back to MMFF94s...")
                opt_method = "MMFF94S"

        # Try MACE-OFF23 if requested and available
        if opt_method == "MACE-OFF23":
            if has_mace_available():
                if verbose:
                    print(f"  Optimizing {len(cids)} conformers with MACE-OFF23...")
                    t0 = time.perf_counter()
                try:
                    mol = optimize_conformer_ensemble_mace(mol, verbose=verbose)
                    if verbose:
                        print(f"  Optimization complete in {time.perf_counter() - t0:.1f}s")
                except Exception as e:
                    if verbose:
                        print(f"  MACE-OFF23 failed ({e}), falling back to MMFF94s...")
                    opt_method = "MMFF94S"  # Fall back to MMFF94s
            else:
                if verbose:
                    print("  MACE-OFF23 not available, falling back to MMFF94s...")
                opt_method = "MMFF94S"

        # MMFF94s optimization (default or fallback)
        if opt_method == "MMFF94S":
            if verbose:
                print(f"  Optimizing {len(cids)} conformers with MMFF94s...")
                t0 = time.perf_counter()
            has_mmff = AllChem.MMFFHasAllMoleculeParams(mol)
            if has_mmff:
                try:
                    opt_results = AllChem.MMFFOptimizeMoleculeConfs(
                        mol, mmffVariant="MMFF94s", numThreads=0
                    )
                    for i, (converged, energy) in enumerate(opt_results):
                        if i < len(cids):
                            conf = mol.GetConformer(cids[i])
                            conf.SetDoubleProp("energy", energy)
                except Exception:
                    has_mmff = False  # fall through to UFF below
            if not has_mmff:
                if verbose:
                    print("  MMFF94s unavailable, using UFF...")
                try:
                    opt_results = AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=0)
                    for i, (converged, energy) in enumerate(opt_results):
                        if i < len(cids):
                            conf = mol.GetConformer(cids[i])
                            conf.SetDoubleProp("energy", energy)
                except Exception:
                    if verbose:
                        print("  Force field optimization failed, continuing without energies")
            if verbose:
                print(f"  Optimization complete in {time.perf_counter() - t0:.1f}s")

    return mol


def filter_by_energy(mol: Chem.Mol, energy_window_kcal: float = 10.0) -> Chem.Mol:
    """
    Remove conformers outside an energy window from the global minimum.

    At 300K, kT ~ 0.6 kcal/mol. Default 10 kcal/mol window covers all
    thermally accessible states while excluding strained/unphysical geometries.

    Args:
        mol: Molecule with conformers that have "energy" property set.
        energy_window_kcal: Maximum energy above minimum to keep (default: 10.0).

    Returns:
        Molecule with high-energy conformers removed.
    """
    energies = []
    for conf in mol.GetConformers():
        try:
            energy = conf.GetDoubleProp("energy")
            energies.append((conf.GetId(), energy))
        except KeyError:
            # No energy — keep the conformer
            energies.append((conf.GetId(), 0.0))

    if not energies:
        return mol

    min_energy = min(e for _, e in energies)
    keep_ids = [
        cid for cid, e in energies if (e - min_energy) <= energy_window_kcal
    ]

    if len(keep_ids) == len(energies):
        return mol  # nothing to remove

    # Rebuild molecule with contiguous conformer IDs (0, 1, 2, ...)
    # This avoids "Bad Conformer Id" errors in downstream code that
    # iterates conformers using range(mol.GetNumConformers()).
    new_mol = Chem.RWMol(mol)
    new_mol.RemoveAllConformers()
    for cid in keep_ids:
        conf = mol.GetConformer(cid)
        new_mol.AddConformer(conf, assignId=True)

    return new_mol.GetMol()


def find_representative_atoms(
    mol: Chem.Mol,
    aromatic_atoms: List[int],
) -> List[int]:
    """
    Find bridgehead atoms (ring-junction carbons) in a fused aromatic system.

    These atoms appear in 3+ SSSR rings, making them geometrically central
    and ideal anchors for distance constraints in biased conformer generation.
    Works on 2D molecules (topological, no 3D coordinates needed).

    For pyrene (4 fused rings): returns ~4 bridgehead carbons.
    For naphthalene (2 fused rings): returns ~2 junction carbons.
    Fallback: atoms with highest ring membership count.

    Args:
        mol: RDKit molecule (2D or 3D).
        aromatic_atoms: Atom indices of one aromatic ring system.

    Returns:
        List of 1-4 representative atom indices, most central first.
    """
    aromatic_set = set(aromatic_atoms)
    ring_info = mol.GetRingInfo()

    # Count how many SSSR rings each aromatic atom belongs to
    ring_counts: Dict[int, int] = {}
    for ring in ring_info.AtomRings():
        for atom_idx in ring:
            if atom_idx in aromatic_set:
                ring_counts[atom_idx] = ring_counts.get(atom_idx, 0) + 1

    if not ring_counts:
        # Fallback: take middle atoms from the sorted list
        mid = len(aromatic_atoms) // 2
        return aromatic_atoms[max(0, mid - 1) : mid + 2]

    # Sort by ring count descending (most central first)
    sorted_atoms = sorted(ring_counts.items(), key=lambda x: -x[1])

    # Prefer atoms in 3+ rings (true bridgeheads in fused systems)
    bridgeheads = [idx for idx, count in sorted_atoms if count >= 3]
    if not bridgeheads:
        # Fallback: atoms in 2+ rings (junction atoms)
        bridgeheads = [idx for idx, count in sorted_atoms if count >= 2]
    if not bridgeheads:
        # Ultimate fallback: atom with highest ring membership
        bridgeheads = [sorted_atoms[0][0]]

    return bridgeheads[:4]


def generate_conformers_biased(
    mol: Chem.Mol,
    num_confs: int = 100,
    aromatic_system: str = "pyrene",
    constraint_levels: Optional[List[Optional[Tuple[float, float]]]] = None,
    prune_rms: float = 0.5,
    random_seed: int = 42,
    optimize: bool = True,
    verbose: bool = False,
    max_iterations: int = 1000,
    use_macrocycle_torsions: bool = True,
    bounds_mat_force_scaling: float = 1.0,
    per_level_energy_filter: bool = True,
    energy_window_kcal: float = 10.0,
    optimizer: str = "MMFF94s",
) -> Chem.Mol:
    """
    Generate conformers biased toward pi-stacking using bounds matrix constraints.

    Compensates for MMFF94s missing London dispersion forces by constraining
    inter-aromatic distances during ETKDGv3 embedding. Multiple constraint
    levels ensure diverse sampling across the stacking landscape.

    Algorithm:
        1. Identify aromatic ring systems (works on 2D molecules)
        2. Find bridgehead atoms as constraint anchors
        3. For each constraint level, modify the bounds matrix and embed
        4. Merge all conformers and optimize with MMFF94s
        5. Filter by energy per constraint level (not globally)

    Constraint levels (defaults):
        - Tight (3.0-4.5 A): targets strong excimer geometries
        - Medium (4.0-6.0 A): targets weak excimer geometries
        - Loose (5.0-8.0 A): captures intermediate geometries
        - None: unconstrained baseline for diversity

    After MMFF94s optimization (which lacks dispersion), stacking distances
    relax by ~1-2 A. This is expected and acceptable because we use geometry
    classification, not energy, as the primary metric.

    Energy filtering note:
        MMFF94s lacks London dispersion, so pi-stacked conformers (tight
        constraint levels) have artificially high energies compared to
        unconstrained ones. Global energy filtering removes exactly the
        conformers we want. Per-level filtering (default) applies the energy
        window within each constraint level independently, preserving
        conformer diversity across the stacking landscape.

    Args:
        mol: RDKit molecule with Hs. Must contain two aromatic ring systems.
        num_confs: Total conformers to generate (split across constraint levels).
        aromatic_system: Aromatic system name for ring identification.
        constraint_levels: List of (lower, upper) distance bounds in Angstroms,
            or None for unconstrained. Default: 4 levels from tight to free.
        prune_rms: RMSD threshold for deduplication within each level.
        random_seed: Random seed for reproducibility.
        optimize: If True, minimize with MMFF94s after embedding.
        verbose: If True, print progress messages with timing.
        max_iterations: Maximum distance geometry iterations per conformer.
            For large molecules (100+ atoms), 1000+ is recommended. Default: 1000.
        use_macrocycle_torsions: Use macrocycle torsion preferences in ETKDGv3.
            Recommended True for large flexible molecules. Default: True.
        bounds_mat_force_scaling: Scaling factor for distance constraint
            stiffness. Higher values (2.0-5.0) better enforce constraints
            but may reduce embedding success rate. Default: 1.0.
        per_level_energy_filter: If True, apply energy window independently
            per constraint level (prevents MMFF94s dispersion bias from
            removing pi-stacked conformers). If False, apply globally.
            Default: True.
        energy_window_kcal: Energy window in kcal/mol for per-level filtering.
            Only used when per_level_energy_filter=True. Default: 10.0.
        optimizer: Optimization method ("MMFF94s", "GFN2-xTB", "MACE-OFF23", or "none").
            - "MMFF94s": Fast but lacks London dispersion (pi-stacking too long)
            - "GFN2-xTB": Dispersion-corrected, correct pi-stacking
            - "MACE-OFF23": DFT-quality dispersion, 10-20x faster than xTB (recommended)
            - "none": Skip optimization (fastest, but less accurate)
            Default: "MMFF94s" for backward compatibility.

    Returns:
        Molecule with conformers from all constraint levels.
        Properties set: "energy" (kcal/mol), "constraint_level" (string).

    Raises:
        ValueError: If two aromatic ring systems cannot be identified.

    References:
        - Landrum (2021) "ETKDG and distance constraints" (RDKit Blog)
        - Landrum (2024) "Conformer Generation and Intramolecular H-bonds"
    """
    # Ensure Hs
    if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        mol = Chem.AddHs(mol)
    mol.RemoveAllConformers()

    # Default constraint levels
    if constraint_levels is None:
        constraint_levels = [
            (3.0, 4.5),  # tight — strong excimer target
            (4.0, 6.0),  # medium — weak excimer target
            (5.0, 8.0),  # loose — intermediate
            None,         # unconstrained baseline
        ]

    # Step 1: Identify aromatic systems (works on 2D)
    analyzer = AromaticDimerAnalyzer(
        aromatic_system=aromatic_system, verbose=False
    )
    try:
        ring1_atoms, ring2_atoms = analyzer.identify_aromatic_rings(mol)
    except Exception as e:
        raise ValueError(
            f"Cannot identify two aromatic ring systems: {e}"
        )

    # Step 2: Find representative atoms for constraints
    rep1 = find_representative_atoms(mol, ring1_atoms)
    rep2 = find_representative_atoms(mol, ring2_atoms)
    constraint_pairs = [(a, b) for a in rep1 for b in rep2]

    if not constraint_pairs:
        warnings.warn("No constraint atom pairs found, falling back to standard generation")
        return generate_conformers(mol, num_confs, prune_rms, random_seed, optimize,
                                   verbose=verbose)

    # Step 3: Generate conformers at each constraint level
    confs_per_level = max(num_confs // len(constraint_levels), 5)
    all_conformers = []  # list of (Conformer, level_str)

    if verbose:
        n_levels = len(constraint_levels)
        n_heavy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
        est = _estimate_pipeline_time(n_heavy, num_confs, True)
        print(f"  Generating {num_confs} conformers across {n_levels} constraint levels")
        print(f"  Molecule: {n_heavy} heavy atoms, {len(constraint_pairs)} constraint pairs")
        print(f"  Estimated time: {est}")

    for level_idx, level in enumerate(constraint_levels):
        seed = random_seed + level_idx * 1000
        level_str = f"{level[0]}-{level[1]}" if level else "none"
        level_name = f"{level[0]}-{level[1]} A" if level else "unconstrained"

        if verbose:
            print(f"  Level {level_idx + 1}/{len(constraint_levels)}: {level_name}...")
            t_level = time.perf_counter()

        # Create a fresh working copy for this level
        work_mol = Chem.RWMol(mol)
        work_mol.RemoveAllConformers()

        if level is None:
            # Unconstrained baseline
            params = rdDistGeom.ETKDGv3()
            params.randomSeed = seed
            params.numThreads = 0
            params.pruneRmsThresh = prune_rms
            params.useRandomCoords = True
            params.enforceChirality = True
            params.maxIterations = max_iterations
            if use_macrocycle_torsions:
                params.useMacrocycleTorsions = True
            cids = rdDistGeom.EmbedMultipleConfs(
                work_mol, numConfs=confs_per_level, params=params
            )
        else:
            lower_bound, upper_bound = level
            try:
                bounds = rdDistGeom.GetMoleculeBoundsMatrix(work_mol)
            except Exception:
                warnings.warn(f"Could not get bounds matrix for level {level_str}")
                continue

            # Set constraints for all representative atom pairs
            for a, b in constraint_pairs:
                i, j = max(a, b), min(a, b)
                # Lower triangle (i>j) = lower bound
                bounds[i, j] = max(bounds[i, j], lower_bound)
                # Upper triangle (j<i, i.e. j,i where j<i) = upper bound
                # Note: bounds[j, i] where j < i is the upper bound
                bounds[j, i] = min(bounds[j, i], upper_bound)
                # Ensure lower <= upper
                if bounds[i, j] > bounds[j, i]:
                    bounds[j, i] = bounds[i, j] + 0.1

            # Triangle smoothing — CRITICAL for constraint propagation
            DoTriangleSmoothing(bounds)

            params = rdDistGeom.ETKDGv3()
            params.randomSeed = seed
            params.numThreads = 0
            params.pruneRmsThresh = prune_rms
            params.useRandomCoords = True
            params.enforceChirality = True
            params.maxIterations = max_iterations
            if use_macrocycle_torsions:
                params.useMacrocycleTorsions = True
            if bounds_mat_force_scaling != 1.0:
                params.boundsMatForceScaling = bounds_mat_force_scaling
            params.SetBoundsMat(bounds)

            cids = rdDistGeom.EmbedMultipleConfs(
                work_mol, numConfs=confs_per_level, params=params
            )

            if len(cids) == 0:
                # Retry with relaxed bounds
                if verbose:
                    print(f"    Retrying with relaxed bounds...")
                relaxed_lower = max(lower_bound - 0.5, 2.5)
                relaxed_upper = upper_bound + 1.0
                bounds2 = rdDistGeom.GetMoleculeBoundsMatrix(work_mol)
                for a, b in constraint_pairs:
                    i, j = max(a, b), min(a, b)
                    bounds2[i, j] = max(bounds2[i, j], relaxed_lower)
                    bounds2[j, i] = min(bounds2[j, i], relaxed_upper)
                    if bounds2[i, j] > bounds2[j, i]:
                        bounds2[j, i] = bounds2[i, j] + 0.1
                DoTriangleSmoothing(bounds2)
                params2 = rdDistGeom.ETKDGv3()
                params2.randomSeed = seed + 500
                params2.numThreads = 0
                params2.useRandomCoords = True
                params2.SetBoundsMat(bounds2)
                cids = rdDistGeom.EmbedMultipleConfs(
                    work_mol, numConfs=confs_per_level, params=params2
                )

        # Collect conformers with level tags
        for cid in cids:
            conf = work_mol.GetConformer(cid)
            all_conformers.append((Chem.Conformer(conf), level_str))

        if verbose:
            elapsed = time.perf_counter() - t_level
            print(f"    -> {len(cids)} conformers in {elapsed:.1f}s")

    if not all_conformers:
        warnings.warn("Biased generation produced no conformers, trying standard")
        return generate_conformers(mol, num_confs, prune_rms, random_seed, optimize,
                                   verbose=verbose, optimizer=optimizer)

    # Step 4: Merge all conformers into the original molecule
    mol.RemoveAllConformers()
    for conf, level_str in all_conformers:
        new_id = mol.AddConformer(conf, assignId=True)
        mol.GetConformer(new_id).SetProp("constraint_level", level_str)

    # Step 5: Optimization
    if optimize and optimizer.lower() != "none":
        total_confs = mol.GetNumConformers()
        opt_method = optimizer.upper()

        # Try GFN2-xTB if requested and available
        if opt_method == "GFN2-XTB":
            if has_xtb_available():
                if verbose:
                    print(f"  Optimizing {total_confs} conformers with GFN2-xTB (slowest step)...")
                    t_opt = time.perf_counter()
                try:
                    mol = optimize_conformer_ensemble(mol, method="GFN2-xTB", verbose=verbose)
                    if verbose:
                        print(f"  Optimization complete in {time.perf_counter() - t_opt:.1f}s")
                except Exception as e:
                    if verbose:
                        print(f"  GFN2-xTB failed ({e}), falling back to MMFF94s...")
                    opt_method = "MMFF94S"  # Fall back to MMFF94s
            else:
                if verbose:
                    print("  GFN2-xTB not available, falling back to MMFF94s...")
                opt_method = "MMFF94S"

        # Try MACE-OFF23 if requested and available
        if opt_method == "MACE-OFF23":
            if has_mace_available():
                if verbose:
                    print(f"  Optimizing {total_confs} conformers with MACE-OFF23 (slowest step)...")
                    t_opt = time.perf_counter()
                try:
                    mol = optimize_conformer_ensemble_mace(mol, verbose=verbose)
                    if verbose:
                        print(f"  Optimization complete in {time.perf_counter() - t_opt:.1f}s")
                except Exception as e:
                    if verbose:
                        print(f"  MACE-OFF23 failed ({e}), falling back to MMFF94s...")
                    opt_method = "MMFF94S"  # Fall back to MMFF94s
            else:
                if verbose:
                    print("  MACE-OFF23 not available, falling back to MMFF94s...")
                opt_method = "MMFF94S"

        # MMFF94s optimization (default or fallback)
        if opt_method == "MMFF94S":
            if verbose:
                print(f"  Optimizing {total_confs} conformers with MMFF94s (slowest step)...")
                t_opt = time.perf_counter()
            has_mmff = AllChem.MMFFHasAllMoleculeParams(mol)
            if has_mmff:
                try:
                    opt_results = AllChem.MMFFOptimizeMoleculeConfs(
                        mol, mmffVariant="MMFF94s", numThreads=0
                    )
                    for i, (converged, energy) in enumerate(opt_results):
                        if i < mol.GetNumConformers():
                            conf = mol.GetConformer(i)
                            conf.SetDoubleProp("energy", energy)
                except Exception:
                    has_mmff = False  # fall through to UFF below
            if not has_mmff:
                if verbose:
                    print("  MMFF94s unavailable, using UFF...")
                try:
                    opt_results = AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=0)
                    for i, (converged, energy) in enumerate(opt_results):
                        if i < mol.GetNumConformers():
                            conf = mol.GetConformer(i)
                            conf.SetDoubleProp("energy", energy)
                except Exception:
                    if verbose:
                        print("  Force field optimization failed, continuing without energies")
            if verbose:
                print(f"  Optimization complete in {time.perf_counter() - t_opt:.1f}s")

    # Step 6: Per-level energy filtering
    # MMFF94s lacks London dispersion, so pi-stacked conformers (tight
    # constraint levels) have artificially high energies compared to
    # unconstrained ones. Global energy filtering removes exactly the
    # conformers we want to keep. Per-level filtering applies the energy
    # window within each constraint level independently.
    if per_level_energy_filter and mol.GetNumConformers() > 0:
        if verbose:
            n_before = mol.GetNumConformers()
            print(f"  Per-level energy filtering ({energy_window_kcal} kcal/mol)...")

        # Group conformer IDs by constraint level
        level_confs: Dict[str, List[Tuple[int, float]]] = {}
        for conf in mol.GetConformers():
            cid = conf.GetId()
            try:
                level = conf.GetProp("constraint_level")
            except KeyError:
                level = "unknown"
            try:
                energy = conf.GetDoubleProp("energy")
            except (KeyError, SystemError):
                energy = 0.0
            level_confs.setdefault(level, []).append((cid, energy))

        # Keep conformers within energy window per level
        keep_ids = []
        for level, confs in level_confs.items():
            if not confs:
                continue
            min_energy = min(e for _, e in confs)
            level_keep = [
                cid for cid, e in confs
                if (e - min_energy) <= energy_window_kcal
            ]
            keep_ids.extend(level_keep)
            if verbose:
                print(f"    Level {level}: {len(level_keep)}/{len(confs)} kept "
                      f"(min={min_energy:.1f}, range={max(e for _, e in confs) - min_energy:.1f})")

        # Rebuild molecule with kept conformers
        if len(keep_ids) < mol.GetNumConformers():
            new_mol = Chem.RWMol(mol)
            new_mol.RemoveAllConformers()
            for cid in keep_ids:
                conf = mol.GetConformer(cid)
                new_mol.AddConformer(conf, assignId=True)
            mol = new_mol.GetMol()

        if verbose:
            print(f"  Kept {mol.GetNumConformers()}/{n_before} conformers after per-level filter")

    return mol


def aggregate_results(results_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate per-conformer results into per-substituent summary.

    Computes comprehensive ensemble features (distributional statistics,
    threshold-based features, and Boltzmann-weighted averages) plus
    backward-compatible summary columns.

    Backward-compatible columns (always present):
        - n_conformers, n_excimer, excimer_fraction
        - mean_angle, mean_distance, mean_overlap, best_overlap
        - lowest_energy_class (if energy data available)

    New columns (from ensemble module):
        - Distributional: {descriptor}_{mean,std,min,max,median,p10,p25,p75,p90}
        - Threshold: frac_strong_excimer, frac_weak_excimer, frac_any_excimer, etc.
        - Boltzmann: {descriptor}_boltz, boltz_excimer_fraction, boltz_excimer_score

    Args:
        results_df: DataFrame with per-conformer analysis results.
            Must have columns: substituent or molecule, conformer_id,
            classification, plane_angle_deg, interplane_distance_A,
            pi_overlap_pct.

    Returns:
        Summary DataFrame indexed by group column, sorted by excimer_fraction.
    """
    if results_df.empty:
        return pd.DataFrame()

    # Auto-detect group column
    if "substituent" in results_df.columns:
        group_col = "substituent"
    else:
        group_col = "molecule"

    # Compute ensemble features via ensemble module
    from pyrene_analyzer.ensemble import compute_ensemble_features

    summary = compute_ensemble_features(
        results_df, group_col=group_col
    )

    # Add backward-compatible column aliases
    alias_map = {
        "mean_angle": "plane_angle_deg_mean",
        "mean_distance": "interplane_distance_A_mean",
        "mean_overlap": "pi_overlap_pct_mean",
        "best_overlap": "pi_overlap_pct_max",
        "excimer_fraction": "frac_any_excimer",
    }
    for alias, source in alias_map.items():
        if source in summary.columns:
            summary[alias] = summary[source]

    # Compute n_excimer from n_conformers and excimer_fraction
    if "n_conformers" in summary.columns and "excimer_fraction" in summary.columns:
        summary["n_excimer"] = (
            summary["n_conformers"] * summary["excimer_fraction"]
        ).round().astype(int)

    # Excimer score columns
    grouped = results_df.groupby(group_col)
    if "excimer_score" in results_df.columns:
        summary["mean_score"] = grouped["excimer_score"].mean()
        summary["best_score"] = grouped["excimer_score"].max()

    # Secondary metric: classification of lowest-energy conformer
    if "energy_kcal_mol" in results_df.columns:
        valid_energy = results_df.dropna(subset=["energy_kcal_mol"])
        if not valid_energy.empty:
            idx = valid_energy.groupby(group_col)["energy_kcal_mol"].idxmin()
            lowest_e = valid_energy.loc[idx]
            summary["lowest_energy_class"] = (
                lowest_e.set_index(group_col)["classification"]
            )

    if "excimer_fraction" in summary.columns:
        summary = summary.sort_values("excimer_fraction", ascending=False)
    return summary


class SubstituentScreener:
    """
    Virtual screening by R-group enumeration + conformer generation + analysis.

    Pipeline:
        1. Load template dimer molecule from SDF
        2. For each substituent, replace R-groups using ReplaceSubstructs
        3. Generate conformer ensemble with ETKDGv3 + MMFF94s
        4. Filter by energy window (default: 10 kcal/mol)
        5. Analyze each conformer with AromaticDimerAnalyzer
        6. Classify and rank substituents by excimer-forming fraction

    Scientific basis:
        R-group steric → bridge dihedrals → θ angle → π-overlap → excimer
        (Ge et al. 2020, Basuroy et al. 2021)

    Example:
        >>> screener = SubstituentScreener("dimer.sdf", "[CH2][CH3]")
        >>> results, summary = screener.screen(num_confs=50)
        >>> print(summary[["excimer_fraction", "mean_overlap", "lowest_energy_class"]])
    """

    def __init__(
        self,
        template_sdf: Union[str, Path],
        r_group_smarts: str,
        aromatic_system: str = "pyrene",
        verbose: bool = True,
    ):
        """
        Initialize the screener.

        Args:
            template_sdf: Path to SDF file with the template dimer molecule.
                The first molecule in the file is used as template.
            r_group_smarts: SMARTS pattern matching the R-group to replace.
                Example: "[CH2][CH3]" for ethyl groups.
            aromatic_system: Aromatic system for classification thresholds.
            verbose: Print progress messages.
        """
        self.template_sdf = Path(template_sdf)
        self.r_group_smarts = r_group_smarts
        self.aromatic_system = aromatic_system
        self.verbose = verbose

        # Load template molecule
        self._template = self._load_template()
        if self._template is None:
            raise ValueError(f"Could not load template from {template_sdf}")

        # Verify R-group pattern matches
        pattern = Chem.MolFromSmarts(r_group_smarts)
        if pattern is None:
            raise ValueError(f"Invalid R-group SMARTS: {r_group_smarts}")

        matches = self._template.GetSubstructMatches(pattern)
        if not matches:
            raise ValueError(
                f"R-group SMARTS '{r_group_smarts}' has no matches in template"
            )

        n_matches = len(matches)
        if self.verbose:
            print(f"Template loaded: {self._template.GetNumAtoms()} atoms")
            print(f"R-group SMARTS '{r_group_smarts}' matches {n_matches} time(s)")

    def _load_template(self) -> Optional[Chem.Mol]:
        """Load the first molecule from the template SDF."""
        try:
            supplier = Chem.SDMolSupplier(str(self.template_sdf), removeHs=False)
            for mol in supplier:
                if mol is not None:
                    return mol
        except OSError:
            return None
        return None

    def _log(self, msg: str) -> None:
        """Print message if verbose."""
        if self.verbose:
            print(msg)

    def enumerate(
        self,
        substituents: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Optional[Chem.Mol]]:
        """
        Generate new molecules by replacing R-groups with each substituent.

        Args:
            substituents: Dict of name→SMILES. If None, uses COMMON_SUBSTITUENTS.

        Returns:
            Dict of substituent_name → RDKit Mol (or None if replacement failed).
        """
        if substituents is None:
            substituents = validate_substituent_library()

        enumerated = {}
        for name, smiles in substituents.items():
            new_mol = replace_r_groups(
                self._template,
                self.r_group_smarts,
                smiles,
                replace_all=True,
            )
            if new_mol is not None:
                enumerated[name] = new_mol
                self._log(f"  {name}: {new_mol.GetNumAtoms()} atoms (OK)")
            else:
                enumerated[name] = None
                self._log(f"  {name}: FAILED")

        n_ok = sum(1 for v in enumerated.values() if v is not None)
        self._log(f"\nEnumeration: {n_ok}/{len(substituents)} succeeded")
        return enumerated

    def screen(
        self,
        substituents: Optional[Dict[str, str]] = None,
        num_confs: int = 50,
        prune_rms: float = 0.5,
        energy_window: float = 10.0,
        random_seed: int = 42,
        use_biased: bool = True,
        optimizer: str = "MMFF94s",
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Full screening pipeline: enumerate → conformers → analyze → rank.

        Args:
            substituents: Dict of name→SMILES. If None, uses COMMON_SUBSTITUENTS.
            num_confs: Number of conformers per molecule.
            prune_rms: RMSD pruning threshold (Angstroms).
            energy_window: Energy window for filtering (kcal/mol).
            random_seed: Random seed for reproducibility.
            use_biased: If True (default), use biased conformer generation with
                distance constraints to improve excimer geometry sampling.
            optimizer: Optimization method ("MMFF94s", "GFN2-xTB", "MACE-OFF23", or "none").
                - "MMFF94s": Fast but lacks London dispersion (pi-stacking too long)
                - "GFN2-xTB": Dispersion-corrected, correct pi-stacking
                - "MACE-OFF23": DFT-quality dispersion, 10-20x faster than xTB (recommended)
                - "none": Skip optimization (fastest, but less accurate)
                Default: "MMFF94s" for backward compatibility.

        Returns:
            Tuple of (all_results_df, summary_df):
            - all_results_df: Per-conformer results for all substituents
            - summary_df: Per-substituent summary ranked by excimer_fraction
        """
        self._log("=" * 60)
        self._log("VIRTUAL SCREENING PIPELINE")
        self._log("=" * 60)
        self._log(f"Template: {self.template_sdf}")
        self._log(f"R-group: {self.r_group_smarts}")
        self._log(f"System: {self.aromatic_system}")
        bias_str = "biased" if use_biased else "standard"
        opt_str = optimizer.upper()
        if opt_str == "GFN2-XTB" and not has_xtb_available():
            opt_str += " (fallback: MMFF94s)"
        elif opt_str == "MACE-OFF23" and not has_mace_available():
            opt_str += " (fallback: MMFF94s)"
        self._log(f"Conformers: {num_confs} ({bias_str}), Optimizer: {opt_str}")
        self._log(f"Energy window: {energy_window} kcal/mol")
        self._log("")

        # Step 1: Enumerate
        self._log("Step 1: Enumerating substituents...")
        enumerated = self.enumerate(substituents)

        # Step 2-5: For each substituent, generate conformers and analyze
        analyzer = AromaticDimerAnalyzer(
            aromatic_system=self.aromatic_system,
            verbose=self.verbose,
        )

        all_results = []

        for name, mol in enumerated.items():
            if mol is None:
                continue

            self._log(f"\nProcessing '{name}'...")

            # Step 2: Prepare molecule
            prepped = prepare_molecule(Chem.RWMol(mol), verbose=self.verbose)
            if prepped is None:
                self._log(f"  Preparation failed, skipping")
                continue

            # Step 3: Generate conformers
            opt_name = optimizer.upper()
            if opt_name == "GFN2-XTB" and not has_xtb_available():
                opt_name = "MMFF94s (fallback)"
            elif opt_name == "MACE-OFF23" and not has_mace_available():
                opt_name = "MMFF94s (fallback)"
            if use_biased:
                self._log(
                    f"  Generating {num_confs} conformers "
                    f"(biased ETKDGv3 + {opt_name})..."
                )
                try:
                    conf_mol = generate_conformers_biased(
                        prepped,
                        num_confs=num_confs,
                        aromatic_system=self.aromatic_system,
                        prune_rms=prune_rms,
                        random_seed=random_seed,
                        verbose=self.verbose,
                        optimizer=optimizer,
                    )
                except ValueError:
                    self._log(f"  Biased generation failed, falling back to standard")
                    conf_mol = generate_conformers(
                        prepped,
                        num_confs=num_confs,
                        prune_rms=prune_rms,
                        random_seed=random_seed,
                        verbose=self.verbose,
                        optimizer=optimizer,
                    )
            else:
                self._log(
                    f"  Generating {num_confs} conformers (ETKDGv3 + {opt_name})..."
                )
                conf_mol = generate_conformers(
                    prepped,
                    num_confs=num_confs,
                    prune_rms=prune_rms,
                    random_seed=random_seed,
                    verbose=self.verbose,
                    optimizer=optimizer,
                )

            n_before = conf_mol.GetNumConformers()

            # Step 4: Energy filter
            conf_mol = filter_by_energy(conf_mol, energy_window)
            n_after = conf_mol.GetNumConformers()
            self._log(
                f"  Conformers: {n_before} generated, "
                f"{n_after} after energy filter ({energy_window} kcal/mol)"
            )

            if n_after == 0:
                self._log(f"  No conformers survived filtering, skipping")
                continue

            # Step 5: Analyze with AromaticDimerAnalyzer
            try:
                mol_results = analyzer.analyze_molecule(
                    conf_mol, mol_name=name, show_progress=self.verbose
                )
                if not mol_results.empty:
                    mol_results = analyzer.add_classification(mol_results)
                    mol_results = analyzer.add_excimer_score(mol_results)
                    mol_results["substituent"] = name
                    all_results.append(mol_results)

                    # Quick summary
                    n_exc = sum(
                        mol_results["classification"].isin(
                            ["strong_excimer", "weak_excimer"]
                        )
                    )
                    frac = n_exc / len(mol_results) if len(mol_results) > 0 else 0
                    self._log(
                        f"  Results: {len(mol_results)} conformers, "
                        f"{n_exc} excimer ({frac:.1%})"
                    )
            except Exception as e:
                self._log(f"  Analysis failed: {e}")
                continue

        if not all_results:
            self._log("\nNo results generated.")
            return pd.DataFrame(), pd.DataFrame()

        # Combine all results
        results_df = pd.concat(all_results, ignore_index=True)

        # Step 6: Aggregate and rank
        self._log("\n" + "=" * 60)
        self._log("SCREENING RESULTS")
        self._log("=" * 60)

        summary_df = aggregate_results(results_df)

        if self.verbose and not summary_df.empty:
            show_score = "mean_score" in summary_df.columns
            header = (
                f"\n{'Substituent':<12} {'Excimer%':>8} {'MeanAngle':>9} "
                f"{'MeanDist':>8} {'MeanOvlp':>8} {'BestOvlp':>8} "
            )
            if show_score:
                header += f"{'MeanScore':>9} {'BestScore':>9} "
            header += f"{'LowEClass':<14}"
            print(header)
            print("-" * (98 if show_score else 80))

            for name, row in summary_df.iterrows():
                low_e_cls = row.get("lowest_energy_class", "N/A")
                line = (
                    f"{name:<12} {row['excimer_fraction']:>7.1%} "
                    f"{row['mean_angle']:>9.1f} {row['mean_distance']:>8.2f} "
                    f"{row['mean_overlap']:>8.1f} {row['best_overlap']:>8.1f} "
                )
                if show_score:
                    line += (
                        f"{row['mean_score']:>9.2f} {row['best_score']:>9.2f} "
                    )
                line += f"{low_e_cls:<14}"
                print(line)

        return results_df, summary_df


def analyze_from_smiles(
    smiles: str,
    aromatic_system: str = "pyrene",
    num_confs: int = 100,
    energy_window: float = 10.0,
    use_biased: bool = True,
    random_seed: int = 42,
    verbose: bool = False,
    optimizer: str = "MMFF94s",
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Analyze an aromatic dimer from SMILES: generate conformers and classify.

    Complete pipeline: SMILES -> 3D -> conformer ensemble -> geometric
    analysis -> excimer classification. No SDF template or MOE required.

    Args:
        smiles: SMILES string of the aromatic dimer molecule.
        aromatic_system: Aromatic system for classification thresholds.
        num_confs: Number of conformers to generate.
        energy_window: Energy window for filtering (kcal/mol).
        use_biased: If True, use biased conformer generation with distance
            constraints to improve excimer geometry sampling.
        random_seed: Random seed for reproducibility.
        verbose: If True, print progress messages with timing for each
            pipeline step.
        optimizer: Optimization method ("MMFF94s", "GFN2-xTB", "MACE-OFF23", or "none").
            - "MMFF94s": Fast but lacks London dispersion (pi-stacking too long)
            - "GFN2-xTB": Dispersion-corrected, correct pi-stacking
            - "MACE-OFF23": DFT-quality dispersion, 10-20x faster than xTB (recommended)
            - "none": Skip optimization (fastest, but less accurate)
            Default: "MMFF94s" for backward compatibility.

    Returns:
        Tuple of (results_df, summary_dict):
        - results_df: Per-conformer analysis with classification columns
        - summary_dict: Summary with excimer_fraction, mean_angle, etc.

    Raises:
        ValueError: If SMILES is invalid or molecule lacks two aromatic systems.
    """
    pipeline_start = time.perf_counter()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    n_heavy = mol.GetNumHeavyAtoms()

    if verbose:
        mode_str = "biased" if use_biased else "standard"
        opt_str = optimizer.upper()
        if opt_str == "GFN2-XTB" and not has_xtb_available():
            opt_str += " (fallback: MMFF94s)"
        elif opt_str == "MACE-OFF23" and not has_mace_available():
            opt_str += " (fallback: MMFF94s)"
        est = _estimate_pipeline_time(n_heavy, num_confs, use_biased, optimizer)
        print("=" * 60)
        print("SMILES ANALYSIS PIPELINE")
        print("=" * 60)
        print(f"Molecule: {n_heavy} heavy atoms")
        print(f"Mode: {num_confs} conformers ({mode_str})")
        print(f"Optimizer: {opt_str}")
        print(f"Estimated time: {est}")
        print()

    # Step 1: Prepare molecule (standardize, add Hs, initial 3D, minimize)
    if verbose:
        print("Step 1/5: Preparing molecule...")
    step_start = time.perf_counter()
    prepped = prepare_molecule(Chem.RWMol(mol), verbose=verbose)
    if prepped is None:
        raise ValueError("Failed to prepare molecule from SMILES")
    if verbose:
        print(f"  Done in {time.perf_counter() - step_start:.1f}s")
        print()

    # Step 2: Generate conformers
    if verbose:
        gen_str = "biased ETKDGv3 + MMFF94s" if use_biased else "ETKDGv3 + MMFF94s"
        print(f"Step 2/5: Generating conformers ({gen_str})...")
    step_start = time.perf_counter()
    if use_biased:
        # Biased generation includes per-level energy filtering internally
        # to prevent MMFF94s dispersion bias from removing pi-stacked conformers
        conf_mol = generate_conformers_biased(
            prepped,
            num_confs=num_confs,
            aromatic_system=aromatic_system,
            random_seed=random_seed,
            verbose=verbose,
            per_level_energy_filter=True,
            energy_window_kcal=energy_window,
            optimizer=optimizer,
        )
    else:
        conf_mol = generate_conformers(
            prepped,
            num_confs=num_confs,
            random_seed=random_seed,
            verbose=verbose,
            optimizer=optimizer,
        )
    n_generated = conf_mol.GetNumConformers()
    if verbose:
        print(f"  Generated {n_generated} conformers in "
              f"{time.perf_counter() - step_start:.1f}s")
        print()

    # Step 3: Energy filter
    # For biased generation, per-level filtering was already done inside
    # generate_conformers_biased(). Apply global filter only for standard mode.
    if use_biased:
        n_after = n_generated
        if verbose:
            print(f"Step 3/5: Energy filtering (per-level, already applied)...")
            print(f"  {n_after} conformers")
            print()
    else:
        if verbose:
            print(f"Step 3/5: Filtering by energy window ({energy_window} kcal/mol)...")
        conf_mol = filter_by_energy(conf_mol, energy_window)
        n_after = conf_mol.GetNumConformers()
        if verbose:
            print(f"  Kept {n_after} of {n_generated} conformers")
            print()

    if n_after == 0:
        if verbose:
            elapsed = time.perf_counter() - pipeline_start
            print(f"Pipeline complete in {elapsed:.1f}s total (no conformers survived)")
        return pd.DataFrame(), {"excimer_fraction": 0.0, "n_conformers": 0}

    # Step 4: Analyze with AromaticDimerAnalyzer
    if verbose:
        print("Step 4/5: Analyzing geometric parameters...")
    step_start = time.perf_counter()
    analyzer = AromaticDimerAnalyzer(
        aromatic_system=aromatic_system, verbose=verbose
    )
    results_df = analyzer.analyze_molecule(
        conf_mol, mol_name="smiles_input", show_progress=verbose
    )
    if verbose:
        print(f"  Done in {time.perf_counter() - step_start:.1f}s")
        print()

    if results_df.empty:
        if verbose:
            elapsed = time.perf_counter() - pipeline_start
            print(f"Pipeline complete in {elapsed:.1f}s total (no results)")
        return pd.DataFrame(), {"excimer_fraction": 0.0, "n_conformers": 0}

    # Step 5: Classify conformers
    if verbose:
        print("Step 5/5: Classifying conformers...")
    results_df = analyzer.add_classification(results_df)
    results_df = analyzer.add_excimer_score(results_df)

    # Build summary dict
    n_total = len(results_df)
    n_excimer = sum(
        results_df["classification"].isin(["strong_excimer", "weak_excimer"])
    )
    summary = {
        "excimer_fraction": n_excimer / n_total if n_total > 0 else 0.0,
        "n_conformers": n_total,
        "n_excimer": n_excimer,
        "mean_angle": results_df["plane_angle_deg"].mean(),
        "mean_distance": results_df["interplane_distance_A"].mean(),
        "mean_overlap": results_df["pi_overlap_pct"].mean(),
        "best_overlap": results_df["pi_overlap_pct"].max(),
    }

    if verbose:
        elapsed = time.perf_counter() - pipeline_start
        print(f"\nPipeline complete in {elapsed:.1f}s total")

    return results_df, summary


def load_substituents_from_file(filepath: Union[str, Path]) -> Dict[str, str]:
    """
    Load a custom substituent library from a JSON file.

    Expected format:
        {"Me": "C", "Et": "CC", "iPr": "C(C)C", ...}

    Args:
        filepath: Path to JSON file with name→SMILES mapping.

    Returns:
        Validated dict of substituent name→SMILES.
    """
    with open(filepath) as f:
        data = json.load(f)

    if not isinstance(data, dict):
        raise ValueError("Substituent file must contain a JSON object (dict)")

    return validate_substituent_library(data)
