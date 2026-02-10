"""
Tests for the virtual screening module.

Tests cover:
    - Substituent library validation
    - R-group replacement
    - Conformer generation
    - Energy filtering
    - Molecule preparation
    - Full screening pipeline integration
    - Results aggregation and ranking
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from pyrene_analyzer.screening import (
    COMMON_SUBSTITUENTS,
    SubstituentScreener,
    aggregate_results,
    analyze_from_smiles,
    filter_by_energy,
    find_representative_atoms,
    generate_conformers,
    generate_conformers_biased,
    load_substituents_from_file,
    prepare_molecule,
    replace_r_groups,
    validate_substituent_library,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def naphthalene_dimer_smiles():
    """SMILES for a simple naphthalene dimer connected by ethylene bridge with ethyl groups."""
    # Two naphthalenes connected by -CH2-CH2- bridge, each with an ethyl substituent
    return "CCc1ccc2ccccc2c1CCc1ccc2ccccc2c1CC"


@pytest.fixture
def simple_ethyl_benzene():
    """A simple ethylbenzene molecule for R-group testing."""
    mol = Chem.MolFromSmiles("CCc1ccccc1")
    return mol


@pytest.fixture
def biphenyl_diethyl():
    """Biphenyl with two ethyl groups (one on each ring)."""
    mol = Chem.MolFromSmiles("CCc1ccc(-c2ccc(CC)cc2)cc1")
    return mol


@pytest.fixture
def test_sdf_path():
    """Path to the test SDF file."""
    path = Path(__file__).parent / "test_data" / "pyrene_dimer_set_for_MOE.sdf"
    if path.exists():
        return path
    pytest.skip("Test SDF file not found")


@pytest.fixture
def tmp_substituents_file(tmp_path):
    """Create a temporary substituents JSON file."""
    data = {"Me": "C", "Et": "CC", "iPr": "C(C)C"}
    filepath = tmp_path / "substituents.json"
    with open(filepath, "w") as f:
        json.dump(data, f)
    return filepath


# ---------------------------------------------------------------------------
# Test: Substituent Library
# ---------------------------------------------------------------------------


class TestSubstituentLibrary:
    """Tests for the built-in substituent library."""

    def test_all_smiles_parse(self):
        """All SMILES in COMMON_SUBSTITUENTS must parse correctly."""
        for name, smiles in COMMON_SUBSTITUENTS.items():
            mol = Chem.MolFromSmiles(smiles)
            assert mol is not None, f"Failed to parse '{name}': {smiles}"

    def test_library_size(self):
        """Library should have ~25 substituents."""
        assert len(COMMON_SUBSTITUENTS) >= 20
        assert len(COMMON_SUBSTITUENTS) <= 30

    def test_validate_all_valid(self):
        """validate_substituent_library returns all entries for valid library."""
        valid = validate_substituent_library()
        assert len(valid) == len(COMMON_SUBSTITUENTS)

    def test_validate_with_invalid(self):
        """validate_substituent_library drops invalid SMILES."""
        lib = {"good": "C", "bad": "XXXX_invalid_XXXX"}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            valid = validate_substituent_library(lib)
        assert "good" in valid
        assert "bad" not in valid

    def test_steric_series_present(self):
        """Library must contain key steric series: H, Me, Et, iPr, tBu."""
        for name in ["H", "Me", "Et", "iPr", "tBu"]:
            assert name in COMMON_SUBSTITUENTS

    def test_electronic_series_present(self):
        """Library must contain key electronic groups: NH2, OMe, F, Cl, NO2, CN."""
        for name in ["NH2", "OMe", "F", "Cl", "NO2", "CN"]:
            assert name in COMMON_SUBSTITUENTS


class TestLoadSubstituentsFromFile:
    """Tests for loading custom substituent files."""

    def test_load_valid_file(self, tmp_substituents_file):
        """Load a valid JSON substituent file."""
        result = load_substituents_from_file(tmp_substituents_file)
        assert "Me" in result
        assert "Et" in result
        assert "iPr" in result

    def test_load_invalid_format(self, tmp_path):
        """Reject non-dict JSON."""
        filepath = tmp_path / "bad.json"
        with open(filepath, "w") as f:
            json.dump(["not", "a", "dict"], f)
        with pytest.raises(ValueError, match="JSON object"):
            load_substituents_from_file(filepath)


# ---------------------------------------------------------------------------
# Test: R-Group Replacement
# ---------------------------------------------------------------------------


class TestRGroupReplacement:
    """Tests for R-group replacement functionality."""

    def test_replace_ethyl_with_methyl(self, simple_ethyl_benzene):
        """Replace ethyl with methyl on ethylbenzene."""
        result = replace_r_groups(simple_ethyl_benzene, "[CH2][CH3]", "C")
        assert result is not None
        # Methyl benzene (toluene) has 7 heavy atoms
        n_heavy = sum(1 for a in result.GetAtoms() if a.GetAtomicNum() != 0)
        assert n_heavy == 7

    def test_replace_ethyl_with_isopropyl(self, simple_ethyl_benzene):
        """Replace ethyl with isopropyl."""
        result = replace_r_groups(simple_ethyl_benzene, "[CH2][CH3]", "C(C)C")
        assert result is not None

    def test_replace_all_symmetric(self, biphenyl_diethyl):
        """Replace all ethyl groups on symmetric biphenyl."""
        result = replace_r_groups(
            biphenyl_diethyl, "[CH2][CH3]", "C", replace_all=True
        )
        assert result is not None

    def test_replace_with_fluorine(self, simple_ethyl_benzene):
        """Replace ethyl with fluorine."""
        result = replace_r_groups(simple_ethyl_benzene, "[CH2][CH3]", "F")
        assert result is not None

    def test_invalid_smarts_returns_none(self, simple_ethyl_benzene):
        """Invalid SMARTS should return None with warning."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = replace_r_groups(
                simple_ethyl_benzene, "INVALID_SMARTS", "C"
            )
        assert result is None

    def test_no_match_returns_molecule(self, simple_ethyl_benzene):
        """SMARTS with no match should still return a molecule."""
        # [SiH3] won't match anything in ethylbenzene
        result = replace_r_groups(simple_ethyl_benzene, "[SiH3]", "C")
        # ReplaceSubstructs returns the original mol when no match found
        assert result is not None


# ---------------------------------------------------------------------------
# Test: Conformer Generation
# ---------------------------------------------------------------------------


class TestConformerGeneration:
    """Tests for conformer generation with ETKDGv3."""

    def test_generate_basic(self):
        """Generate conformers for a simple molecule."""
        mol = Chem.MolFromSmiles("CCCCCC")  # hexane
        mol = Chem.AddHs(mol)
        result = generate_conformers(mol, num_confs=10, prune_rms=0.3)
        assert result.GetNumConformers() > 0

    def test_generate_with_optimization(self):
        """Generated conformers should have energy properties."""
        mol = Chem.MolFromSmiles("CCCCCC")
        mol = Chem.AddHs(mol)
        result = generate_conformers(mol, num_confs=10, optimize=True)
        # At least some conformers should have energy
        has_energy = False
        for conf in result.GetConformers():
            try:
                conf.GetDoubleProp("energy")
                has_energy = True
                break
            except KeyError:
                continue
        assert has_energy

    def test_generate_reproducible(self):
        """Same seed should give same results."""
        mol = Chem.MolFromSmiles("CCCCCC")
        mol1 = Chem.AddHs(Chem.Mol(mol))
        mol2 = Chem.AddHs(Chem.Mol(mol))

        r1 = generate_conformers(mol1, num_confs=5, random_seed=42)
        r2 = generate_conformers(mol2, num_confs=5, random_seed=42)

        assert r1.GetNumConformers() == r2.GetNumConformers()

    def test_generate_aromatic(self):
        """Generate conformers for a biphenyl (aromatic flexibility)."""
        mol = Chem.MolFromSmiles("c1ccc(-c2ccccc2)cc1")
        mol = Chem.AddHs(mol)
        result = generate_conformers(mol, num_confs=10)
        assert result.GetNumConformers() > 0


# ---------------------------------------------------------------------------
# Test: Energy Filtering
# ---------------------------------------------------------------------------


class TestEnergyFiltering:
    """Tests for energy-window-based conformer filtering."""

    def test_filter_removes_high_energy(self):
        """High-energy conformers should be removed."""
        mol = Chem.MolFromSmiles("CCCCCC")
        mol = Chem.AddHs(mol)
        mol = generate_conformers(mol, num_confs=20, optimize=True)

        n_before = mol.GetNumConformers()
        mol = filter_by_energy(mol, energy_window_kcal=2.0)
        n_after = mol.GetNumConformers()

        # With a tight 2 kcal/mol window, some should be removed
        # (unless all conformers have very similar energy)
        assert n_after <= n_before
        assert n_after > 0  # at least the minimum should survive

    def test_filter_keeps_all_in_wide_window(self):
        """Very wide window should keep all conformers."""
        mol = Chem.MolFromSmiles("CCCCCC")
        mol = Chem.AddHs(mol)
        mol = generate_conformers(mol, num_confs=10, optimize=True)

        n_before = mol.GetNumConformers()
        mol = filter_by_energy(mol, energy_window_kcal=1000.0)
        n_after = mol.GetNumConformers()

        assert n_after == n_before

    def test_filter_empty_molecule(self):
        """Molecule with no conformers should be returned unchanged."""
        mol = Chem.MolFromSmiles("C")
        mol = Chem.AddHs(mol)
        result = filter_by_energy(mol, energy_window_kcal=10.0)
        assert result.GetNumConformers() == 0


# ---------------------------------------------------------------------------
# Test: Molecule Preparation
# ---------------------------------------------------------------------------


class TestMoleculePreparation:
    """Tests for the MOE-replacement preparation pipeline."""

    def test_prepare_from_smiles(self):
        """Prepare a molecule from SMILES (no 3D coords)."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = prepare_molecule(mol)
        assert result is not None
        assert result.GetNumConformers() > 0

    def test_prepare_adds_hydrogens(self):
        """Prepared molecule should have explicit Hs."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = prepare_molecule(mol)
        assert result is not None
        h_count = sum(1 for a in result.GetAtoms() if a.GetAtomicNum() == 1)
        assert h_count > 0  # benzene should have 6 Hs

    def test_prepare_none_input(self):
        """None input should return None."""
        assert prepare_molecule(None) is None

    def test_prepare_with_existing_3d(self):
        """Molecule with existing 3D coords should be optimized."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        result = prepare_molecule(mol)
        assert result is not None


# ---------------------------------------------------------------------------
# Test: Results Aggregation
# ---------------------------------------------------------------------------


class TestAggregation:
    """Tests for results aggregation and ranking."""

    def test_aggregate_basic(self):
        """Basic aggregation should produce expected columns."""
        data = {
            "substituent": ["Me", "Me", "Me", "Et", "Et", "Et"],
            "conformer_id": [0, 1, 2, 0, 1, 2],
            "classification": [
                "strong_excimer", "weak_excimer", "monomer",
                "monomer", "monomer", "monomer",
            ],
            "plane_angle_deg": [10.0, 25.0, 70.0, 60.0, 65.0, 80.0],
            "interplane_distance_A": [3.5, 3.8, 5.0, 4.5, 5.0, 6.0],
            "pi_overlap_pct": [60.0, 40.0, 10.0, 15.0, 10.0, 5.0],
            "energy_kcal_mol": [0.0, 1.0, 5.0, 0.0, 2.0, 8.0],
        }
        df = pd.DataFrame(data)
        summary = aggregate_results(df)

        assert "Me" in summary.index
        assert "Et" in summary.index
        assert summary.loc["Me", "excimer_fraction"] == pytest.approx(2 / 3)
        assert summary.loc["Et", "excimer_fraction"] == pytest.approx(0.0)

    def test_aggregate_sorted_by_excimer_fraction(self):
        """Results should be sorted by excimer_fraction descending."""
        data = {
            "substituent": ["A", "A", "B", "B"],
            "conformer_id": [0, 1, 0, 1],
            "classification": ["monomer", "monomer", "strong_excimer", "strong_excimer"],
            "plane_angle_deg": [50.0, 55.0, 10.0, 12.0],
            "interplane_distance_A": [5.0, 5.5, 3.5, 3.4],
            "pi_overlap_pct": [10.0, 8.0, 60.0, 65.0],
        }
        df = pd.DataFrame(data)
        summary = aggregate_results(df)

        assert summary.index[0] == "B"  # B has 100% excimer
        assert summary.index[1] == "A"  # A has 0% excimer

    def test_aggregate_empty_df(self):
        """Empty DataFrame should return empty summary."""
        df = pd.DataFrame()
        summary = aggregate_results(df)
        assert summary.empty

    def test_aggregate_with_lowest_energy_class(self):
        """Should include lowest_energy_class when energy data is present."""
        data = {
            "substituent": ["Me", "Me", "Me"],
            "conformer_id": [0, 1, 2],
            "classification": ["strong_excimer", "weak_excimer", "monomer"],
            "plane_angle_deg": [10.0, 25.0, 70.0],
            "interplane_distance_A": [3.5, 3.8, 5.0],
            "pi_overlap_pct": [60.0, 40.0, 10.0],
            "energy_kcal_mol": [5.0, 0.5, 2.0],
        }
        df = pd.DataFrame(data)
        summary = aggregate_results(df)

        # Lowest energy is conformer 1 (0.5 kcal/mol) = weak_excimer
        assert summary.loc["Me", "lowest_energy_class"] == "weak_excimer"


class TestEnhancedAggregation:
    """Tests for enhanced aggregate_results with ensemble features."""

    def _make_test_data(self, n_per_sub=20):
        """Helper to create test data with all required columns."""
        np.random.seed(42)
        rows = []
        for sub in ["Me", "Et"]:
            for i in range(n_per_sub):
                angle = np.random.uniform(0, 90)
                dist = np.random.uniform(3.0, 6.0)
                overlap = np.random.uniform(0, 100)
                if angle < 20 and 3.3 <= dist <= 3.7 and overlap > 50:
                    cls = "strong_excimer"
                elif angle < 60 and dist < 4.5 and overlap > 30:
                    cls = "weak_excimer"
                else:
                    cls = "monomer"
                rows.append(
                    {
                        "substituent": sub,
                        "conformer_id": i,
                        "plane_angle_deg": angle,
                        "interplane_distance_A": dist,
                        "pi_overlap_pct": overlap,
                        "centroid_distance_A": dist + np.random.uniform(0, 2),
                        "slip_stack_A": np.random.uniform(0, 3),
                        "energy_kcal_mol": np.random.uniform(300, 320),
                        "classification": cls,
                        "excimer_score": np.random.uniform(0, 1),
                    }
                )
        return pd.DataFrame(rows)

    def test_backward_compatible_columns(self):
        """Enhanced aggregation should still have old column names."""
        df = self._make_test_data()
        summary = aggregate_results(df)
        for col in [
            "mean_angle", "mean_distance", "mean_overlap",
            "best_overlap", "excimer_fraction",
        ]:
            assert col in summary.columns, f"Missing backward-compatible column: {col}"

    def test_new_distributional_columns(self):
        """Should have distributional stats from ensemble module."""
        df = self._make_test_data()
        summary = aggregate_results(df)
        assert "plane_angle_deg_mean" in summary.columns
        assert "interplane_distance_A_std" in summary.columns

    def test_new_threshold_columns(self):
        """Should have threshold-based features from ensemble module."""
        df = self._make_test_data()
        summary = aggregate_results(df)
        assert "n_conformers" in summary.columns
        assert "frac_any_excimer" in summary.columns

    def test_new_boltzmann_columns(self):
        """Should have Boltzmann-weighted features from ensemble module."""
        df = self._make_test_data()
        summary = aggregate_results(df)
        assert "plane_angle_deg_boltz" in summary.columns

    def test_backward_compatible_values(self):
        """Old column values should match old behavior."""
        df = self._make_test_data()
        summary = aggregate_results(df)

        # mean_angle should equal plane_angle_deg_mean
        for sub in ["Me", "Et"]:
            if sub in summary.index:
                assert summary.loc[sub, "mean_angle"] == pytest.approx(
                    summary.loc[sub, "plane_angle_deg_mean"]
                )
                assert summary.loc[sub, "mean_distance"] == pytest.approx(
                    summary.loc[sub, "interplane_distance_A_mean"]
                )

    def test_auto_detect_molecule_col(self):
        """Should auto-detect 'molecule' column when 'substituent' not present."""
        df = self._make_test_data()
        df = df.rename(columns={"substituent": "molecule"})
        summary = aggregate_results(df)
        assert "Me" in summary.index
        assert "Et" in summary.index

    def test_excimer_score_columns(self):
        """Should include excimer score aggregations."""
        df = self._make_test_data()
        summary = aggregate_results(df)
        assert "mean_score" in summary.columns
        assert "best_score" in summary.columns

    def test_lowest_energy_class(self):
        """Should include lowest_energy_class when energy data present."""
        df = self._make_test_data()
        summary = aggregate_results(df)
        assert "lowest_energy_class" in summary.columns


# ---------------------------------------------------------------------------
# Test: SubstituentScreener (integration)
# ---------------------------------------------------------------------------


class TestSubstituentScreener:
    """Integration tests for the full screening pipeline."""

    def test_init_with_test_sdf(self, test_sdf_path):
        """Initialize screener with the test SDF file."""
        screener = SubstituentScreener(
            template_sdf=test_sdf_path,
            r_group_smarts="[CH2][CH3]",
            verbose=False,
        )
        assert screener._template is not None

    def test_init_invalid_sdf(self, tmp_path):
        """Invalid SDF should raise ValueError."""
        bad_sdf = tmp_path / "empty.sdf"
        bad_sdf.write_text("")
        with pytest.raises(ValueError, match="Could not load"):
            SubstituentScreener(bad_sdf, "[CH2][CH3]", verbose=False)

    def test_init_no_match_smarts(self, test_sdf_path):
        """SMARTS with no match should raise ValueError."""
        with pytest.raises(ValueError, match="no matches"):
            SubstituentScreener(test_sdf_path, "[SiH3][GeH3]", verbose=False)

    def test_enumerate_small_set(self, test_sdf_path):
        """Enumerate a small set of substituents."""
        screener = SubstituentScreener(
            template_sdf=test_sdf_path,
            r_group_smarts="[CH2][CH3]",
            verbose=False,
        )
        subs = {"Me": "C", "F": "F"}
        enumerated = screener.enumerate(subs)
        # At least some should succeed
        n_ok = sum(1 for v in enumerated.values() if v is not None)
        assert n_ok > 0

    @pytest.mark.slow
    def test_screen_small(self, test_sdf_path):
        """Full screening pipeline with minimal parameters."""
        screener = SubstituentScreener(
            template_sdf=test_sdf_path,
            r_group_smarts="[CH2][CH3]",
            verbose=False,
        )
        subs = {"Me": "C", "iPr": "C(C)C"}
        results_df, summary_df = screener.screen(
            substituents=subs,
            num_confs=5,  # very few for speed
            energy_window=10.0,
        )
        # May or may not produce results depending on molecule complexity
        # Just check it runs without error
        assert isinstance(results_df, pd.DataFrame)
        assert isinstance(summary_df, pd.DataFrame)

    @pytest.mark.slow
    def test_screen_with_no_bias(self, test_sdf_path):
        """Screening with use_biased=False should use standard generation."""
        screener = SubstituentScreener(
            template_sdf=test_sdf_path,
            r_group_smarts="[CH2][CH3]",
            verbose=False,
        )
        subs = {"Me": "C"}
        results_df, summary_df = screener.screen(
            substituents=subs,
            num_confs=5,
            energy_window=10.0,
            use_biased=False,
        )
        assert isinstance(results_df, pd.DataFrame)
        assert isinstance(summary_df, pd.DataFrame)


# ---------------------------------------------------------------------------
# Test: find_representative_atoms
# ---------------------------------------------------------------------------


class TestFindRepresentativeAtoms:
    """Tests for bridgehead atom identification."""

    def test_pyrene_bridgeheads(self):
        """Pyrene should have atoms in 3+ rings (bridgeheads)."""
        mol = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
        assert mol is not None
        aromatic_atoms = list(range(mol.GetNumHeavyAtoms()))
        reps = find_representative_atoms(mol, aromatic_atoms)
        assert len(reps) >= 1
        assert len(reps) <= 4
        # Bridgehead atoms should be in multiple rings
        ring_info = mol.GetRingInfo()
        for idx in reps:
            count = sum(1 for ring in ring_info.AtomRings() if idx in ring)
            assert count >= 2

    def test_naphthalene_junctions(self):
        """Naphthalene should have junction atoms in 2+ rings."""
        mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
        assert mol is not None
        aromatic_atoms = list(range(mol.GetNumHeavyAtoms()))
        reps = find_representative_atoms(mol, aromatic_atoms)
        assert len(reps) >= 1
        assert len(reps) <= 4
        # Junction atoms should be in 2 rings
        ring_info = mol.GetRingInfo()
        for idx in reps:
            count = sum(1 for ring in ring_info.AtomRings() if idx in ring)
            assert count >= 2

    def test_benzene_fallback(self):
        """Single ring should fall back gracefully."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        assert mol is not None
        aromatic_atoms = list(range(mol.GetNumHeavyAtoms()))
        reps = find_representative_atoms(mol, aromatic_atoms)
        # Should return at least 1 atom
        assert len(reps) >= 1

    def test_returns_bounded(self):
        """Should return at most 4 atoms."""
        mol = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")  # pyrene
        assert mol is not None
        aromatic_atoms = list(range(mol.GetNumHeavyAtoms()))
        reps = find_representative_atoms(mol, aromatic_atoms)
        assert len(reps) <= 4

    def test_empty_atoms_fallback(self):
        """Empty aromatic_atoms list should not crash."""
        mol = Chem.MolFromSmiles("CCCCCC")  # no aromatic atoms
        reps = find_representative_atoms(mol, [])
        # Should return empty or a small fallback list
        assert isinstance(reps, list)


# ---------------------------------------------------------------------------
# Test: generate_conformers_biased
# ---------------------------------------------------------------------------


class TestGenerateConformersBiased:
    """Tests for biased conformer generation with distance constraints."""

    def test_basic_generation(self):
        """Should produce >0 conformers for a simple dimer."""
        # Use a small naphthalene dimer
        smiles = "c1ccc2ccccc2c1CCc1ccc2ccccc2c1"
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        result = generate_conformers_biased(
            mol, num_confs=20, aromatic_system="naphthalene",
            optimize=False, random_seed=42,
        )
        assert result.GetNumConformers() > 0

    def test_constraint_levels_tagged(self):
        """Conformers should have constraint_level property."""
        smiles = "c1ccc2ccccc2c1CCc1ccc2ccccc2c1"
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        result = generate_conformers_biased(
            mol, num_confs=20, aromatic_system="naphthalene",
            optimize=False, random_seed=42,
        )
        # At least some conformers should have constraint_level property
        has_level = False
        for conf in result.GetConformers():
            try:
                conf.GetProp("constraint_level")
                has_level = True
                break
            except KeyError:
                continue
        assert has_level

    def test_custom_constraint_levels(self):
        """Should accept custom constraint bounds."""
        smiles = "c1ccc2ccccc2c1CCc1ccc2ccccc2c1"
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        custom_levels = [(3.5, 5.0), None]
        result = generate_conformers_biased(
            mol, num_confs=10, aromatic_system="naphthalene",
            constraint_levels=custom_levels,
            optimize=False, random_seed=42,
        )
        assert result.GetNumConformers() > 0

    def test_optimization_adds_energy(self):
        """Optimized conformers should have energy property."""
        smiles = "c1ccc2ccccc2c1CCc1ccc2ccccc2c1"
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        result = generate_conformers_biased(
            mol, num_confs=12, aromatic_system="naphthalene",
            optimize=True, random_seed=42,
        )
        has_energy = False
        for conf in result.GetConformers():
            try:
                conf.GetDoubleProp("energy")
                has_energy = True
                break
            except KeyError:
                continue
        assert has_energy

    def test_without_optimization(self):
        """Should work with optimize=False."""
        smiles = "c1ccc2ccccc2c1CCc1ccc2ccccc2c1"
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        result = generate_conformers_biased(
            mol, num_confs=8, aromatic_system="naphthalene",
            optimize=False, random_seed=42,
        )
        assert result.GetNumConformers() > 0

    def test_reproducible_seed(self):
        """Same seed should give same number of conformers."""
        smiles = "c1ccc2ccccc2c1CCc1ccc2ccccc2c1"

        mol1 = Chem.MolFromSmiles(smiles)
        mol1 = Chem.AddHs(mol1)
        AllChem.EmbedMolecule(mol1, AllChem.ETKDGv3())

        mol2 = Chem.MolFromSmiles(smiles)
        mol2 = Chem.AddHs(mol2)
        AllChem.EmbedMolecule(mol2, AllChem.ETKDGv3())

        r1 = generate_conformers_biased(
            mol1, num_confs=12, aromatic_system="naphthalene",
            optimize=False, random_seed=42,
        )
        r2 = generate_conformers_biased(
            mol2, num_confs=12, aromatic_system="naphthalene",
            optimize=False, random_seed=42,
        )
        assert r1.GetNumConformers() == r2.GetNumConformers()

    def test_raises_for_non_dimer(self):
        """Should raise ValueError for molecule without two aromatic systems."""
        mol = Chem.MolFromSmiles("c1ccccc1")  # single benzene
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        with pytest.raises(ValueError, match="Cannot identify"):
            generate_conformers_biased(
                mol, num_confs=5, aromatic_system="naphthalene",
                optimize=False,
            )


# ---------------------------------------------------------------------------
# Test: analyze_from_smiles
# ---------------------------------------------------------------------------


class TestAnalyzeFromSmiles:
    """Tests for SMILES-to-analysis pipeline."""

    def test_basic_pipeline(self):
        """Should return (DataFrame, dict) for a valid dimer SMILES."""
        smiles = "c1ccc2ccccc2c1CCc1ccc2ccccc2c1"
        results_df, summary = analyze_from_smiles(
            smiles, aromatic_system="naphthalene",
            num_confs=8, use_biased=False, random_seed=42,
        )
        assert isinstance(results_df, pd.DataFrame)
        assert isinstance(summary, dict)
        assert "excimer_fraction" in summary
        assert "n_conformers" in summary

    def test_invalid_smiles_raises(self):
        """Invalid SMILES should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid SMILES"):
            analyze_from_smiles("XXXX_INVALID_XXXX")

    def test_non_dimer_raises(self):
        """Molecule without two aromatic systems should raise ValueError."""
        with pytest.raises(ValueError):
            analyze_from_smiles(
                "c1ccccc1",  # single benzene
                aromatic_system="naphthalene",
                num_confs=5,
            )

    def test_biased_runs(self):
        """Biased mode should run without error."""
        smiles = "c1ccc2ccccc2c1CCc1ccc2ccccc2c1"
        results_df, summary = analyze_from_smiles(
            smiles, aromatic_system="naphthalene",
            num_confs=12, use_biased=True, random_seed=42,
        )
        assert isinstance(results_df, pd.DataFrame)
        assert summary["n_conformers"] >= 0


# ---------------------------------------------------------------------------
# Test: MACE-OFF23 integration in screening pipeline
# ---------------------------------------------------------------------------


class TestMACEScreeningPaths:
    """Tests for MACE-OFF23 code paths in the screening pipeline.

    Uses mocking to test wiring without requiring MACE installation.
    """

    def test_generate_conformers_mace_fallback_when_unavailable(self):
        """When MACE unavailable, generate_conformers should fall back to MMFF94s."""
        from unittest.mock import patch

        mol = Chem.MolFromSmiles("CCCCCC")
        mol = Chem.AddHs(mol)

        with patch("pyrene_analyzer.screening.has_mace_available", return_value=False):
            result = generate_conformers(
                mol, num_confs=5, optimizer="MACE-OFF23", verbose=False,
            )
        # Should still produce conformers (via MMFF94s fallback)
        assert result.GetNumConformers() > 0

    def test_generate_conformers_mace_fallback_verbose_message(self, capsys):
        """When MACE unavailable + verbose, should print fallback message."""
        from unittest.mock import patch

        mol = Chem.MolFromSmiles("CCCCCC")
        mol = Chem.AddHs(mol)

        with patch("pyrene_analyzer.screening.has_mace_available", return_value=False):
            result = generate_conformers(
                mol, num_confs=5, optimizer="MACE-OFF23", verbose=True,
            )
        captured = capsys.readouterr()
        assert "not available" in captured.out.lower() or "fallback" in captured.out.lower()

    def test_generate_conformers_mace_exception_fallback(self):
        """When MACE raises exception, should fall back to MMFF94s."""
        from unittest.mock import patch

        mol = Chem.MolFromSmiles("CCCCCC")
        mol = Chem.AddHs(mol)

        with patch("pyrene_analyzer.screening.has_mace_available", return_value=True), \
             patch("pyrene_analyzer.screening.optimize_conformer_ensemble_mace",
                   side_effect=RuntimeError("MACE calculation failed")):
            result = generate_conformers(
                mol, num_confs=5, optimizer="MACE-OFF23", verbose=False,
            )
        assert result.GetNumConformers() > 0

    def test_generate_conformers_biased_mace_fallback(self):
        """Biased generation with unavailable MACE should fall back."""
        from unittest.mock import patch

        smiles = "c1ccc2ccccc2c1CCc1ccc2ccccc2c1"
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        with patch("pyrene_analyzer.screening.has_mace_available", return_value=False):
            result = generate_conformers_biased(
                mol, num_confs=12, aromatic_system="naphthalene",
                optimizer="MACE-OFF23", optimize=True, random_seed=42,
            )
        assert result.GetNumConformers() > 0

    def test_analyze_from_smiles_mace_fallback_display(self, capsys):
        """analyze_from_smiles with MACE + verbose should show fallback info."""
        from unittest.mock import patch

        smiles = "c1ccc2ccccc2c1CCc1ccc2ccccc2c1"
        with patch("pyrene_analyzer.screening.has_mace_available", return_value=False):
            results_df, summary = analyze_from_smiles(
                smiles, aromatic_system="naphthalene",
                num_confs=8, use_biased=False, random_seed=42,
                verbose=True, optimizer="MACE-OFF23",
            )
        captured = capsys.readouterr()
        assert "fallback" in captured.out.lower() or "MMFF94s" in captured.out
        assert isinstance(results_df, pd.DataFrame)

    def test_valid_optimizers_tuple(self):
        """VALID_OPTIMIZERS should include all 4 options."""
        from pyrene_analyzer.screening import VALID_OPTIMIZERS
        assert "MMFF94s" in VALID_OPTIMIZERS
        assert "GFN2-xTB" in VALID_OPTIMIZERS
        assert "MACE-OFF23" in VALID_OPTIMIZERS
        assert "none" in VALID_OPTIMIZERS


class TestVerboseMode:
    """Tests that verbose=True doesn't break any function."""

    DIMER_SMILES = "c1ccc2ccccc2c1CCc1ccc2ccccc2c1"

    def test_prepare_molecule_verbose(self):
        """prepare_molecule with verbose=True should still return a valid mol."""
        mol = Chem.MolFromSmiles(self.DIMER_SMILES)
        result = prepare_molecule(Chem.RWMol(mol), verbose=True)
        assert result is not None

    def test_generate_conformers_verbose(self):
        """generate_conformers with verbose=True should produce conformers."""
        mol = Chem.MolFromSmiles(self.DIMER_SMILES)
        prepped = prepare_molecule(Chem.RWMol(mol))
        conf_mol = generate_conformers(prepped, num_confs=5, verbose=True)
        assert conf_mol.GetNumConformers() > 0

    def test_generate_conformers_biased_verbose(self):
        """generate_conformers_biased with verbose=True should produce conformers."""
        mol = Chem.MolFromSmiles(self.DIMER_SMILES)
        prepped = prepare_molecule(Chem.RWMol(mol))
        conf_mol = generate_conformers_biased(
            prepped, num_confs=8, aromatic_system="naphthalene", verbose=True,
        )
        assert conf_mol.GetNumConformers() > 0

    def test_analyze_from_smiles_verbose(self):
        """analyze_from_smiles with verbose=True should return results."""
        results_df, summary = analyze_from_smiles(
            self.DIMER_SMILES, aromatic_system="naphthalene",
            num_confs=8, use_biased=False, random_seed=42, verbose=True,
        )
        assert isinstance(results_df, pd.DataFrame)
        assert summary["n_conformers"] >= 0
