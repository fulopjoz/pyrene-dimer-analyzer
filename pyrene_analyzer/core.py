"""
Core Module
===========

Main analyzer class for aromatic dimer conformational analysis.

This module provides the AromaticDimerAnalyzer class (with backward-compatible
alias PyreneDimerAnalyzer), which orchestrates the complete analysis workflow
including:
- Automatic aromatic ring detection using SMARTS patterns
- Support for multiple aromatic systems (pyrene, perylene, anthracene, etc.)
- Geometric analysis of all conformers
- Batch processing of multiple molecules
- Progress tracking and parallel processing
"""

import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from tqdm import tqdm

from pyrene_analyzer.geometry import (
    DistanceCalibration,
    calculate_centroid_distance,
    calculate_interplane_distance,
    calculate_pi_overlap,
    calculate_plane_angle,
    calculate_slip_stack_displacement,
)
from pyrene_analyzer.aromatic_systems import (
    AROMATIC_SYSTEMS,
    AromaticSystem,
    ClassificationThresholds,
    get_system,
)
from pyrene_analyzer.io import load_molecules

# SMARTS pattern for pyrene (4 fused 6-membered aromatic rings)
# Kept for backward compatibility. Canonical source is now aromatic_systems.py.
PYRENE_SMARTS = "c1cc2ccc3cccc4ccc(c1)c2c34"

# Alternative: simpler pattern for fused aromatic systems
FUSED_AROMATIC_SMARTS = "c1ccc2ccccc2c1"  # Naphthalene-like


class AromaticDimerAnalyzer:
    """
    Analyzer for aromatic dimer geometric properties.

    This class provides methods for analyzing the geometric relationships
    between two aromatic ring systems in a covalently-linked dimer, including
    plane-plane angles, inter-plane distances, and π-overlap calculations.

    Supports multiple aromatic systems (pyrene, perylene, anthracene,
    naphthalene, phenanthrene) with system-specific classification thresholds.
    Custom SMARTS patterns can also be provided for arbitrary aromatic systems.

    Attributes:
        verbose: If True, print progress messages.
        use_smarts: If True, use SMARTS patterns for ring detection.
        use_shapely: If True, use Shapely for accurate π-overlap calculation.

    Example:
        >>> from pyrene_analyzer import AromaticDimerAnalyzer
        >>> analyzer = AromaticDimerAnalyzer(aromatic_system="pyrene")
        >>> results = analyzer.analyze_file('conformers.sdf')
        >>> results.to_csv('analysis.csv')

        >>> # Use a different system
        >>> analyzer = AromaticDimerAnalyzer(aromatic_system="perylene")

        >>> # Use custom SMARTS
        >>> analyzer = AromaticDimerAnalyzer(custom_smarts="c1cc2ccc3cccc4ccc(c1)c2c34")
    """

    def __init__(
        self,
        aromatic_system: str = "pyrene",
        verbose: bool = True,
        use_smarts: bool = True,
        use_shapely: bool = True,
        custom_smarts: Optional[str] = None,
        custom_thresholds: Optional[ClassificationThresholds] = None,
        distance_calibration: Optional[DistanceCalibration] = None,
    ):
        """
        Initialize the AromaticDimerAnalyzer.

        Args:
            aromatic_system: Name of the aromatic system to analyze.
                Supported: "pyrene", "perylene", "anthracene",
                "naphthalene", "phenanthrene". Default: "pyrene".
            verbose: If True, print progress messages during analysis.
            use_smarts: If True, use SMARTS patterns for ring detection.
                        If False, use connectivity-based clustering.
            use_shapely: If True, use Shapely library for accurate π-overlap
                         calculation. If False, use centroid-based approximation.
            custom_smarts: Optional custom SMARTS pattern. If provided,
                overrides the system's default SMARTS.
            custom_thresholds: Optional custom classification thresholds.
                If provided with custom_smarts, uses these thresholds.
                Otherwise falls back to the selected system's thresholds.
            distance_calibration: Optional distance calibration to apply
                to inter-plane distances. Use :func:`make_offset_calibration`
                to create a simple offset correction.
        """
        self.verbose = verbose
        self.use_smarts = use_smarts
        self.use_shapely = use_shapely
        self._distance_calibration = distance_calibration

        # Load system definition from registry
        if custom_smarts is not None:
            self._system_name = "custom"
            self._smarts = custom_smarts
            self._pattern = Chem.MolFromSmarts(custom_smarts)
            if custom_thresholds is not None:
                self._thresholds = custom_thresholds
            else:
                self._thresholds = get_system("pyrene").thresholds
            self._min_ring_atoms = 6
        else:
            system = get_system(aromatic_system)
            self._system_name = system.name
            self._smarts = system.smarts
            self._pattern = Chem.MolFromSmarts(system.smarts)
            self._thresholds = system.thresholds
            self._min_ring_atoms = system.min_ring_atoms

        # Backward-compatible attributes
        self._pyrene_pattern = self._pattern
        self._fused_pattern = Chem.MolFromSmarts(FUSED_AROMATIC_SMARTS)

    def log(self, message: str) -> None:
        """Print message if verbose mode is enabled."""
        if self.verbose:
            print(f"[INFO] {message}")

    def identify_aromatic_rings(self, mol: Chem.Mol) -> Tuple[List[int], List[int]]:
        """
        Identify atoms belonging to each aromatic ring system.

        This method uses two strategies:
        1. SMARTS pattern matching for the configured aromatic substructure
        2. Connectivity-based clustering of aromatic rings (fallback)

        Args:
            mol: RDKit molecule object.

        Returns:
            Tuple of two lists containing atom indices for each aromatic system.

        Raises:
            ValueError: If two aromatic systems cannot be identified.

        Example:
            >>> ring1, ring2 = analyzer.identify_aromatic_rings(mol)
            >>> print(f"Ring system 1: {len(ring1)} atoms")
        """
        self.log(f"Identifying {self._system_name} ring systems...")

        # Try SMARTS-based detection first
        if self.use_smarts and self._pattern is not None:
            matches = [set(m) for m in mol.GetSubstructMatches(self._pattern)]
            if len(matches) >= 2:
                clusters = self._merge_connected_rings(matches)
                clusters = [c for c in clusters if len(c) >= self._min_ring_atoms]

                if len(clusters) >= 2:
                    clusters = sorted(clusters, key=len, reverse=True)
                    ring1 = sorted(list(clusters[0]))
                    ring2 = sorted(list(clusters[1]))
                    self.log(
                        f"Found 2 {self._system_name} systems via SMARTS matching"
                    )
                    return ring1, ring2

        # Fallback to connectivity-based detection
        return self._identify_by_connectivity(mol)

    def identify_pyrene_rings(self, mol: Chem.Mol) -> Tuple[List[int], List[int]]:
        """
        Identify atoms belonging to each pyrene ring system.

        .. deprecated::
            Use :meth:`identify_aromatic_rings` instead.

        Args:
            mol: RDKit molecule object.

        Returns:
            Tuple of two lists containing atom indices for each aromatic system.
        """
        return self.identify_aromatic_rings(mol)

    def _identify_by_connectivity(self, mol: Chem.Mol) -> Tuple[List[int], List[int]]:
        """
        Identify pyrene rings by clustering connected aromatic atoms.

        Args:
            mol: RDKit molecule object.

        Returns:
            Tuple of two lists containing atom indices for each pyrene.
        """
        # Get all aromatic 6-membered rings
        ring_info = mol.GetRingInfo()
        aromatic_rings = []

        for ring in ring_info.AtomRings():
            if len(ring) == 6:
                if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                    aromatic_rings.append(set(ring))

        if len(aromatic_rings) < 6:
            self.log(f"Warning: Found only {len(aromatic_rings)} aromatic rings")

        # Merge connected rings into clusters
        clusters = self._merge_connected_rings(aromatic_rings)

        # Filter clusters by minimum atom count for the configured system
        pyrene_clusters = [c for c in clusters if len(c) >= self._min_ring_atoms]

        # If we don't have 2 pyrene-sized clusters, try with smaller threshold
        if len(pyrene_clusters) < 2:
            pyrene_clusters = sorted(clusters, key=len, reverse=True)[:2]

        # Should have 2 clusters (2 pyrenes)
        if len(pyrene_clusters) < 2:
            raise ValueError(
                f"Could not identify 2 pyrene systems. "
                f"Found {len(clusters)} aromatic cluster(s)."
            )

        if len(pyrene_clusters) > 2:
            self.log(
                f"Warning: Found {len(pyrene_clusters)} aromatic clusters, using 2 largest"
            )
            pyrene_clusters = sorted(pyrene_clusters, key=len, reverse=True)[:2]

        clusters = pyrene_clusters

        pyrene1 = sorted(list(clusters[0]))
        pyrene2 = sorted(list(clusters[1]))

        self.log(f"Pyrene 1: {len(pyrene1)} atoms")
        self.log(f"Pyrene 2: {len(pyrene2)} atoms")

        return pyrene1, pyrene2

    def _merge_connected_rings(self, rings: List[Set[int]]) -> List[Set[int]]:
        """
        Merge rings that share atoms into connected clusters.

        Args:
            rings: List of sets, each containing atom indices of a ring.

        Returns:
            List of merged clusters.
        """
        if not rings:
            return []

        clusters = [ring.copy() for ring in rings]

        # Keep merging until no more changes
        changed = True
        while changed:
            changed = False
            new_clusters = []
            used = set()

            for i, cluster_i in enumerate(clusters):
                if i in used:
                    continue

                merged = cluster_i.copy()
                for j, cluster_j in enumerate(clusters):
                    if j <= i or j in used:
                        continue

                    if merged & cluster_j:  # If clusters share atoms
                        merged.update(cluster_j)
                        used.add(j)
                        changed = True

                new_clusters.append(merged)
                used.add(i)

            clusters = new_clusters

        return clusters

    def find_bridge_atoms(
        self, mol: Chem.Mol, pyrene1: List[int], pyrene2: List[int]
    ) -> Tuple[Optional[List[int]], Optional[List[int]]]:
        """
        Find atoms forming the bridge(s) between two pyrene rings.

        The bridge atoms are non-aromatic atoms that connect the two
        pyrene systems. For symmetric dimers, there are typically two
        bridges (left and right).

        Args:
            mol: RDKit molecule object.
            pyrene1: Atom indices of first pyrene.
            pyrene2: Atom indices of second pyrene.

        Returns:
            Tuple of two lists containing atom indices for each bridge,
            or (None, None) if bridges cannot be identified.

        Notes:
            Bridge atoms are identified by finding non-aromatic atoms
            that are connected to both pyrene systems through a path.
        """
        pyrene1_set = set(pyrene1)
        pyrene2_set = set(pyrene2)

        # Find all non-pyrene, non-hydrogen atoms
        bridge_candidates = []
        for atom_idx in range(mol.GetNumAtoms()):
            if atom_idx in pyrene1_set or atom_idx in pyrene2_set:
                continue

            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == "H":
                continue

            bridge_candidates.append(atom_idx)

        if not bridge_candidates:
            return None, None

        # Find paths connecting pyrenes through bridge atoms
        bridges = []
        visited = set()

        for start_atom in bridge_candidates:
            if start_atom in visited:
                continue

            # BFS to find connected bridge atoms
            path = []
            queue = [start_atom]
            local_visited = set()

            while queue:
                current = queue.pop(0)
                if current in local_visited:
                    continue
                local_visited.add(current)

                if current in bridge_candidates:
                    path.append(current)
                    visited.add(current)

                    atom = mol.GetAtomWithIdx(current)
                    for neighbor in atom.GetNeighbors():
                        n_idx = neighbor.GetIdx()
                        if n_idx not in local_visited:
                            queue.append(n_idx)

            if path:
                bridges.append(path)

        # Return up to 2 bridges
        if len(bridges) >= 2:
            return bridges[0][:4], bridges[1][:4]
        elif len(bridges) == 1:
            # Try to split single bridge
            if len(bridges[0]) >= 8:
                mid = len(bridges[0]) // 2
                return bridges[0][:4], bridges[0][mid : mid + 4]
            return bridges[0][:4] if len(bridges[0]) >= 4 else None, None

        return None, None

    def _get_conformer_energy(self, mol: Chem.Mol, conf_id: int) -> Optional[float]:
        """
        Extract energy value from conformer properties.

        Args:
            mol: RDKit molecule object.
            conf_id: Conformer ID.

        Returns:
            Energy in kcal/mol, or None if not available.
        """
        conf = mol.GetConformer(conf_id)

        # Try conformer properties first
        try:
            if conf.HasProp("energy"):
                return conf.GetDoubleProp("energy")
        except Exception:
            pass

        # Try molecule properties
        for prop_name in [
            "energy",
            "Energy",
            "E",
            "ENERGY",
            "mmff94_energy",
            "rel_energy",
            "dE",
            "potential_energy",
        ]:
            if mol.HasProp(prop_name):
                try:
                    return float(mol.GetProp(prop_name))
                except ValueError:
                    pass

        return None

    def analyze_conformer(
        self,
        mol: Chem.Mol,
        conf_id: int,
        pyrene1_atoms: List[int],
        pyrene2_atoms: List[int],
    ) -> Dict[str, Optional[float]]:
        """
        Analyze geometric properties of a single conformer.

        Args:
            mol: RDKit molecule object.
            conf_id: Conformer ID to analyze.
            pyrene1_atoms: Atom indices of first pyrene.
            pyrene2_atoms: Atom indices of second pyrene.

        Returns:
            Dictionary containing:
                - conformer_id: Conformer index
                - plane_angle_deg: Angle between pyrene planes (0-90°)
                - interplane_distance_A: Perpendicular distance (Å)
                - pi_overlap_pct: π-overlap percentage (0-100%)
                - centroid_distance_A: Centroid-to-centroid distance (Å)
                - slip_stack_A: Lateral displacement (Å)
                - bridge_dihedral_L_deg: Left bridge dihedral (degrees)
                - bridge_dihedral_R_deg: Right bridge dihedral (degrees)
                - energy_kcal_mol: Conformer energy (kcal/mol)
        """
        conf = mol.GetConformer(conf_id)

        # Get 3D coordinates for each pyrene
        coords1 = np.array([list(conf.GetAtomPosition(i)) for i in pyrene1_atoms])
        coords2 = np.array([list(conf.GetAtomPosition(i)) for i in pyrene2_atoms])

        # Calculate geometric properties
        theta = calculate_plane_angle(coords1, coords2)
        raw_distance = calculate_interplane_distance(
            coords1, coords2, plane_angle=theta, warn=False
        )
        if self._distance_calibration is not None:
            distance = self._distance_calibration.apply(
                raw_distance, plane_angle=theta
            )
        else:
            distance = raw_distance
        overlap = calculate_pi_overlap(coords1, coords2, use_shapely=self.use_shapely)
        centroid_dist = calculate_centroid_distance(coords1, coords2)
        slip_stack = calculate_slip_stack_displacement(coords1, coords2)

        # Build geometry warnings
        geometry_warnings = []
        if theta > self._thresholds.high_angle_warning:
            geometry_warnings.append(
                f"interplane_distance unreliable at angle {theta:.1f} deg "
                f"(>{self._thresholds.high_angle_warning:.0f} deg)"
            )
        if theta > self._thresholds.high_angle_warning and distance < 2.0:
            geometry_warnings.append(
                "CRITICAL: likely measuring edge-to-face, not pi-pi distance"
            )

        # Calculate bridge dihedrals
        phi_L = phi_R = None
        bridge1, bridge2 = self.find_bridge_atoms(mol, pyrene1_atoms, pyrene2_atoms)

        try:
            if bridge1 and len(bridge1) >= 4:
                phi_L = rdMolTransforms.GetDihedralDeg(
                    conf, bridge1[0], bridge1[1], bridge1[2], bridge1[3]
                )
            if bridge2 and len(bridge2) >= 4:
                phi_R = rdMolTransforms.GetDihedralDeg(
                    conf, bridge2[0], bridge2[1], bridge2[2], bridge2[3]
                )
        except Exception:
            pass  # Bridge dihedrals are optional

        # Get energy
        energy = self._get_conformer_energy(mol, conf_id)

        result = {
            "conformer_id": conf_id,
            "plane_angle_deg": theta,
            "interplane_distance_A": distance,
            "pi_overlap_pct": overlap,
            "centroid_distance_A": centroid_dist,
            "slip_stack_A": slip_stack,
            "bridge_dihedral_L_deg": phi_L,
            "bridge_dihedral_R_deg": phi_R,
            "energy_kcal_mol": energy,
            "geometry_warnings": (
                "; ".join(geometry_warnings) if geometry_warnings else None
            ),
        }
        if self._distance_calibration is not None:
            result["raw_interplane_distance_A"] = raw_distance
        return result

    def analyze_molecule(
        self, mol: Chem.Mol, mol_name: str = "molecule", show_progress: bool = True
    ) -> pd.DataFrame:
        """
        Analyze all conformers of a pyrene dimer molecule.

        Args:
            mol: RDKit molecule with conformers.
            mol_name: Name/ID for the molecule.
            show_progress: If True, show progress bar.

        Returns:
            DataFrame with geometric analysis for all conformers.

        Example:
            >>> df = analyzer.analyze_molecule(mol, "Et_variant")
            >>> print(df.describe())
        """
        self.log(f"\n{'='*60}")
        self.log(f"Analyzing: {mol_name}")
        self.log(f"{'='*60}")

        # Identify aromatic rings (only once per molecule)
        pyrene1, pyrene2 = self.identify_aromatic_rings(mol)

        # Analyze each conformer
        results = []
        num_confs = mol.GetNumConformers()

        self.log(f"Processing {num_confs} conformers...")

        iterator = range(num_confs)
        if show_progress:
            print(f"Analyzing {mol_name}...", flush=True)
            iterator = tqdm(
                iterator,
                desc=f"Analyzing {mol_name}",
                unit="conf",
                mininterval=0.1,
            )

        for conf_id in iterator:
            conf_results = self.analyze_conformer(mol, conf_id, pyrene1, pyrene2)
            conf_results["molecule"] = mol_name
            results.append(conf_results)

        df = pd.DataFrame(results)

        # Calculate relative energy if energy column exists
        if "energy_kcal_mol" in df.columns and df["energy_kcal_mol"].notna().any():
            min_energy = df["energy_kcal_mol"].min()
            df["rel_energy_kcal_mol"] = df["energy_kcal_mol"] - min_energy

        # Reorder columns
        col_order = [
            "molecule",
            "conformer_id",
            "plane_angle_deg",
            "interplane_distance_A",
            "pi_overlap_pct",
            "centroid_distance_A",
            "slip_stack_A",
            "bridge_dihedral_L_deg",
            "bridge_dihedral_R_deg",
            "energy_kcal_mol",
            "geometry_warnings",
        ]
        if "rel_energy_kcal_mol" in df.columns:
            col_order.append("rel_energy_kcal_mol")

        # Only include columns that exist
        col_order = [c for c in col_order if c in df.columns]
        remaining = [c for c in df.columns if c not in col_order]
        df = df[col_order + remaining]

        self.log(f"\nCompleted analysis of {mol_name}")
        self._print_summary(df)

        return df

    def _print_summary(self, df: pd.DataFrame) -> None:
        """Print summary statistics for analysis results."""
        self.log("Results summary:")
        self.log(
            f"  Plane angle range: {df['plane_angle_deg'].min():.1f} - "
            f"{df['plane_angle_deg'].max():.1f}°"
        )
        self.log(
            f"  Distance range: {df['interplane_distance_A'].min():.2f} - "
            f"{df['interplane_distance_A'].max():.2f} Å"
        )
        self.log(
            f"  pi-overlap range: {df['pi_overlap_pct'].min():.1f} - "
            f"{df['pi_overlap_pct'].max():.1f}%"
        )

        if (
            "rel_energy_kcal_mol" in df.columns
            and df["rel_energy_kcal_mol"].notna().any()
        ):
            self.log(
                f"  Energy range: {df['rel_energy_kcal_mol'].min():.2f} - "
                f"{df['rel_energy_kcal_mol'].max():.2f} kcal/mol"
            )

    def analyze_file(
        self, filename: Union[str, Path], show_progress: bool = True
    ) -> pd.DataFrame:
        """
        Analyze all molecules in a file.

        Args:
            filename: Path to molecular structure file (SDF, MOL2, PDB).
            show_progress: If True, show progress bar.

        Returns:
            DataFrame with analysis results for all molecules and conformers.

        Example:
            >>> results = analyzer.analyze_file('conformers.sdf')
            >>> results.to_csv('analysis.csv')
        """
        if show_progress:
            print("Running analysis...", flush=True)
            print(f"Loading molecules from {filename}...", flush=True)
        molecules = load_molecules(filename)
        if show_progress:
            print(f"Loaded {len(molecules)} molecule(s).", flush=True)

        all_results = []
        for mol, name in molecules:
            try:
                df = self.analyze_molecule(mol, name, show_progress)
                all_results.append(df)
            except ValueError as e:
                self.log(f"Skipping {name}: {e}")
                continue

        if not all_results:
            return pd.DataFrame()

        return pd.concat(all_results, ignore_index=True)

    def batch_analyze(
        self,
        filenames: List[Union[str, Path]],
        n_jobs: int = 1,
        show_progress: bool = True,
    ) -> pd.DataFrame:
        """
        Analyze multiple files with optional parallel processing.

        Args:
            filenames: List of file paths to analyze.
            n_jobs: Number of parallel jobs. Use -1 for all CPUs.
            show_progress: If True, show progress bar.

        Returns:
            DataFrame with combined analysis results.

        Example:
            >>> results = analyzer.batch_analyze(
            ...     ['Et.sdf', 'iPr.sdf', 'cHex.sdf', 'tBu.sdf'],
            ...     n_jobs=4
            ... )
        """
        all_results = []

        if n_jobs == 1:
            # Sequential processing
            iterator = filenames
            if show_progress:
                iterator = tqdm(
                    filenames,
                    desc="Processing files",
                    unit="file",
                    mininterval=0.1,
                )

            for filename in iterator:
                try:
                    df = self.analyze_file(filename, show_progress=False)
                    all_results.append(df)
                except Exception as e:
                    warnings.warn(f"Failed to process {filename}: {e}")
        else:
            # Parallel processing
            if n_jobs == -1:
                n_jobs = None  # Use all CPUs

            with ProcessPoolExecutor(max_workers=n_jobs) as executor:
                futures = {
                    executor.submit(self._analyze_file_worker, f): f for f in filenames
                }

                iterator = as_completed(futures)
                if show_progress:
                    iterator = tqdm(
                        iterator,
                        total=len(filenames),
                        desc="Processing files",
                        unit="file",
                        mininterval=0.1,
                    )

                for future in iterator:
                    filename = futures[future]
                    try:
                        df = future.result()
                        all_results.append(df)
                    except Exception as e:
                        warnings.warn(f"Failed to process {filename}: {e}")

        if not all_results:
            return pd.DataFrame()

        return pd.concat(all_results, ignore_index=True)

    def _analyze_file_worker(self, filename: Union[str, Path]) -> pd.DataFrame:
        """Worker function for parallel file processing."""
        # Create new analyzer instance for worker
        if self._system_name == "custom":
            analyzer = AromaticDimerAnalyzer(
                verbose=False,
                use_smarts=self.use_smarts,
                use_shapely=self.use_shapely,
                custom_smarts=self._smarts,
                custom_thresholds=self._thresholds,
                distance_calibration=self._distance_calibration,
            )
        else:
            analyzer = AromaticDimerAnalyzer(
                aromatic_system=self._system_name,
                verbose=False,
                use_smarts=self.use_smarts,
                use_shapely=self.use_shapely,
                distance_calibration=self._distance_calibration,
            )
        return analyzer.analyze_file(filename, show_progress=False)

    def classify_conformer(
        self, plane_angle: float, distance: float, overlap: float
    ) -> str:
        """
        Classify a conformer based on its geometric properties.

        Uses the classification thresholds for the configured aromatic system.

        Args:
            plane_angle: Plane-plane angle in degrees.
            distance: Inter-plane distance in Angstroms.
            overlap: π-overlap percentage.

        Returns:
            Classification string: 'strong_excimer', 'weak_excimer', or 'monomer'.

        Example:
            >>> classification = analyzer.classify_conformer(15.0, 3.5, 55.0)
            >>> print(classification)  # 'strong_excimer'

        Notes:
            Thresholds are system-specific. For pyrene (default):
            - Strong excimer: θ < 20°, d = 3.3-3.7 Å, overlap > 50%
            - Weak excimer: θ < 60°, d < 4.5 Å, overlap > 30%
            - Monomer: otherwise

            The strong overlap threshold was lowered from 70% to 50%
            based on Ge et al. 2020 and Basuroy et al. 2021.

            Eclipsed conformers (θ < 5°) are flagged with a warning
            because they may be optically dark (zero oscillator strength).

        References:
            - Ge et al. 2020 J. Mater. Chem. C 8, 10223 (overlap 40-80%)
            - Basuroy et al. 2021 J. Chem. Phys. 155, 234304 (42% = excimer)
            - Dai et al. 2024 Molecules 29, 507 (dark eclipsed excimer)
        """
        t = self._thresholds
        d_min, d_max = t.strong_distance_range

        if (
            plane_angle < t.strong_angle_max
            and d_min <= distance <= d_max
            and overlap > t.strong_overlap_min
        ):
            classification = "strong_excimer"
        elif (
            plane_angle < t.weak_angle_max
            and distance < t.weak_distance_max
            and overlap > t.weak_overlap_min
        ):
            classification = "weak_excimer"
        else:
            classification = "monomer"

        # Dark excimer warning: eclipsed geometry may have zero oscillator
        # strength (optically dark). Dai et al. 2024 Molecules 29, 507.
        if classification in ("strong_excimer", "weak_excimer") and plane_angle < 5.0:
            warnings.warn(
                f"Eclipsed excimer (θ={plane_angle:.1f}°) may be optically dark. "
                f"See Dai et al. 2024 Molecules 29, 507.",
                stacklevel=2,
            )

        return classification

    def score_conformer(
        self,
        plane_angle: float,
        distance: float,
        overlap: float,
        slip_stack: Optional[float] = None,
    ) -> float:
        """
        Compute a heuristic excimer-likeness score in [0, 1].

        The score emphasizes π-overlap, then distance and angle. A slip-stack
        penalty is applied when available. This is a ranking aid and does not
        replace classification.

        Args:
            plane_angle: Plane-plane angle in degrees.
            distance: Inter-plane distance in Angstroms.
            overlap: π-overlap percentage.
            slip_stack: Optional slip-stack displacement in Angstroms.

        Returns:
            Score between 0.0 (unlikely excimer) and 1.0 (strong excimer-like).
        """
        t = self._thresholds
        d_min, d_max = t.strong_distance_range

        # Overlap is the primary driver (Ge 2020)
        overlap_score = float(np.clip(overlap / 100.0, 0.0, 1.0))

        # Distance score: peak in strong range, taper to 0 by weak max
        if distance <= 0 or np.isnan(distance):
            distance_score = 0.0
        elif d_min <= distance <= d_max:
            distance_score = 1.0
        elif distance < d_min:
            distance_score = max(0.0, 1.0 - (d_min - distance) / d_min)
        elif distance < t.weak_distance_max:
            distance_score = max(
                0.0, 1.0 - (distance - d_max) / (t.weak_distance_max - d_max)
            )
        else:
            distance_score = 0.0

        # Angle score: strong range best, taper to 0 by weak max
        if plane_angle < 0 or np.isnan(plane_angle):
            angle_score = 0.0
        elif plane_angle <= t.strong_angle_max:
            angle_score = 1.0
        elif plane_angle <= t.weak_angle_max:
            angle_score = max(
                0.0,
                1.0
                - (plane_angle - t.strong_angle_max)
                / (t.weak_angle_max - t.strong_angle_max),
            )
        else:
            angle_score = 0.0

        # Slip-stack penalty (optional)
        if slip_stack is None or np.isnan(slip_stack):
            slip_score = 1.0
        elif slip_stack <= 2.5:
            slip_score = 1.0
        elif slip_stack <= 7.5:
            slip_score = max(0.0, 1.0 - (slip_stack - 2.5) / 5.0)
        else:
            slip_score = 0.0

        # Eclipsed penalty: near-0° conformers may be optically dark
        # (Dai et al. 2024 Molecules 29, 507). Penalize θ < 10° with a
        # Gaussian dip centered at 0°; the emissive sweet-spot is 10-50°.
        if plane_angle < 10.0 and not (np.isnan(plane_angle) or plane_angle < 0):
            eclipsed_factor = 1.0 - 0.5 * np.exp(
                -(plane_angle**2) / (2 * 10.0**2)
            )
        else:
            eclipsed_factor = 1.0

        base_score = 0.6 * overlap_score + 0.2 * distance_score + 0.2 * angle_score
        return float(np.clip(base_score * slip_score * eclipsed_factor, 0.0, 1.0))

    def add_excimer_score(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Add heuristic excimer score to analysis results.

        Args:
            df: DataFrame with analysis results.

        Returns:
            DataFrame with added 'excimer_score' column.
        """
        df = df.copy()
        df["excimer_score"] = df.apply(
            lambda row: self.score_conformer(
                row["plane_angle_deg"],
                row["interplane_distance_A"],
                row["pi_overlap_pct"],
                row.get("slip_stack_A"),
            ),
            axis=1,
        )
        return df

    def add_classification(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Add excimer/monomer classification to analysis results.

        Also adds 'dark_excimer_warning' boolean column for eclipsed
        conformers (θ < 5°) that may be optically dark.

        Args:
            df: DataFrame with analysis results.

        Returns:
            DataFrame with added 'classification' and 'dark_excimer_warning' columns.

        References:
            - Dai et al. 2024 Molecules 29, 507 (eclipsed excimer = dark)
        """
        df = df.copy()
        classifications = []
        dark_warnings = []
        for _, row in df.iterrows():
            angle = row["plane_angle_deg"]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cls = self.classify_conformer(
                    angle,
                    row["interplane_distance_A"],
                    row["pi_overlap_pct"],
                )
            classifications.append(cls)
            # Flag eclipsed excimers as potentially dark
            is_dark = cls in ("strong_excimer", "weak_excimer") and angle < 5.0
            dark_warnings.append(is_dark)
        df["classification"] = classifications
        df["dark_excimer_warning"] = dark_warnings
        return df


# Backward-compatible alias
PyreneDimerAnalyzer = AromaticDimerAnalyzer
