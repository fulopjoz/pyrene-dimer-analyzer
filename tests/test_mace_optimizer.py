"""
Tests for the MACE-OFF23 optimizer module.

Tests are split into:
- TestMACEAvailability: Always runs (checks import machinery)
- TestMACEScreeningIntegration: Always runs (checks pipeline wiring)
- TestMACEOptimizer: Only runs if MACE is installed
"""

import pytest
from unittest.mock import patch, MagicMock

from rdkit import Chem
from rdkit.Chem import AllChem

from pyrene_analyzer.mace_optimizer import (
    has_mace_available,
    get_mace_import_error,
    HAS_MACE,
)


# -------------------------------------------------------------------
# Fixtures
# -------------------------------------------------------------------

@pytest.fixture
def ethane_mol():
    """Simple ethane molecule for basic testing."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


@pytest.fixture
def naphthalene_dimer_mol():
    """Naphthalene dimer for pi-stacking tests."""
    mol = Chem.MolFromSmiles("c1ccc2ccccc2c1CCc1ccc2ccccc2c1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


# -------------------------------------------------------------------
# Tests that always run (no MACE required)
# -------------------------------------------------------------------

class TestMACEAvailability:
    """Tests for MACE availability checks (always run)."""

    def test_has_mace_returns_bool(self):
        result = has_mace_available()
        assert isinstance(result, bool)

    def test_import_error_type(self):
        error = get_mace_import_error()
        assert error is None or isinstance(error, str)

    def test_has_mace_consistent_with_flag(self):
        assert has_mace_available() == HAS_MACE

    def test_import_error_none_when_available(self):
        if HAS_MACE:
            assert get_mace_import_error() is None
        else:
            assert get_mace_import_error() is not None


class TestMACEScreeningIntegration:
    """Tests that MACE is wired into the screening pipeline."""

    def test_valid_optimizers_includes_mace(self):
        from pyrene_analyzer.screening import VALID_OPTIMIZERS
        assert "MACE-OFF23" in VALID_OPTIMIZERS

    def test_mace_optimizer_importable(self):
        """Module should be importable regardless of MACE installation."""
        import pyrene_analyzer.mace_optimizer as mod
        assert hasattr(mod, "optimize_with_mace")
        assert hasattr(mod, "optimize_conformer_ensemble_mace")
        assert hasattr(mod, "filter_by_energy_mace")
        assert hasattr(mod, "single_point_mace")
        assert hasattr(mod, "optimize_or_fallback_mace")

    def test_mace_raises_import_error_when_unavailable(self):
        """When MACE not installed, optimize_with_mace should raise."""
        if HAS_MACE:
            pytest.skip("MACE is installed")
        from pyrene_analyzer.mace_optimizer import optimize_with_mace
        mol = Chem.MolFromSmiles("CC")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        with pytest.raises(ImportError, match="MACE-OFF23 not available"):
            optimize_with_mace(mol)


# -------------------------------------------------------------------
# Tests that require MACE (skipped if not installed)
# -------------------------------------------------------------------

@pytest.mark.skipif(not HAS_MACE, reason="MACE-OFF23 not installed")
class TestMACEOptimizer:
    """Tests for MACE-OFF23 optimization (requires MACE installed)."""

    def test_mol_to_ase_atoms(self, ethane_mol):
        from pyrene_analyzer.mace_optimizer import _mol_to_ase_atoms
        atoms = _mol_to_ase_atoms(ethane_mol)
        assert len(atoms) == ethane_mol.GetNumAtoms()
        assert atoms.get_chemical_symbols()[0] == "C"

    def test_mol_to_ase_atoms_no_conformer(self):
        from pyrene_analyzer.mace_optimizer import _mol_to_ase_atoms
        mol = Chem.MolFromSmiles("CC")
        with pytest.raises(ValueError, match="no conformers"):
            _mol_to_ase_atoms(mol)

    def test_optimize_single_conformer(self, ethane_mol):
        from pyrene_analyzer.mace_optimizer import optimize_with_mace
        mol, energy = optimize_with_mace(ethane_mol, model="small", max_steps=50)
        assert isinstance(energy, float)
        assert energy != 0.0

    def test_optimize_sets_energy_property(self, ethane_mol):
        from pyrene_analyzer.mace_optimizer import optimize_with_mace
        mol, energy = optimize_with_mace(ethane_mol, model="small", max_steps=50)
        conf = mol.GetConformer()
        assert conf.HasProp("energy")
        assert abs(conf.GetDoubleProp("energy") - energy) < 0.01

    def test_optimize_sets_optimizer_property(self, ethane_mol):
        from pyrene_analyzer.mace_optimizer import optimize_with_mace
        mol, _ = optimize_with_mace(ethane_mol, model="small", max_steps=50)
        conf = mol.GetConformer()
        assert conf.HasProp("optimizer")
        assert "MACE-OFF23" in conf.GetProp("optimizer")

    def test_single_point(self, ethane_mol):
        from pyrene_analyzer.mace_optimizer import single_point_mace
        energy = single_point_mace(ethane_mol, model="small")
        assert isinstance(energy, float)
        assert energy != 0.0

    def test_filter_by_energy(self, ethane_mol):
        from pyrene_analyzer.mace_optimizer import filter_by_energy_mace
        # Add a second conformer
        AllChem.EmbedMultipleConfs(ethane_mol, numConfs=3, randomSeed=42)
        # Set fake energies
        for i, conf in enumerate(ethane_mol.GetConformers()):
            conf.SetDoubleProp("energy", float(i * 20))  # 0, 20, 40

        filtered = filter_by_energy_mace(ethane_mol, energy_window_kcal=10.0)
        assert filtered.GetNumConformers() == 1  # only the lowest

    def test_optimize_or_fallback(self, ethane_mol):
        from pyrene_analyzer.mace_optimizer import optimize_or_fallback_mace
        mol, energy, method = optimize_or_fallback_mace(
            ethane_mol, preferred_model="small"
        )
        assert "MACE-OFF23" in method
        assert isinstance(energy, float)

    def test_optimize_no_conformer_raises(self):
        from pyrene_analyzer.mace_optimizer import optimize_with_mace
        mol = Chem.MolFromSmiles("CC")
        mol = Chem.AddHs(mol)
        with pytest.raises(ValueError, match="no conformers"):
            optimize_with_mace(mol)

    def test_optimize_ensemble_returns_mol(self, naphthalene_dimer_mol):
        from pyrene_analyzer.mace_optimizer import optimize_conformer_ensemble_mace
        AllChem.EmbedMultipleConfs(naphthalene_dimer_mol, numConfs=3, randomSeed=42)
        result = optimize_conformer_ensemble_mace(
            naphthalene_dimer_mol, model="small", max_steps=30,
        )
        assert result.GetNumConformers() >= 1

    def test_calculator_cache(self):
        from pyrene_analyzer.mace_optimizer import _get_mace_calculator, _MACE_CALC_CACHE
        calc1 = _get_mace_calculator(model="small", device="cpu")
        calc2 = _get_mace_calculator(model="small", device="cpu")
        # Same key should return cached calculator
        assert ("small", "cpu") in _MACE_CALC_CACHE

    def test_update_mol_coordinates_roundtrip(self, ethane_mol):
        from pyrene_analyzer.mace_optimizer import _mol_to_ase_atoms, _update_mol_from_ase
        atoms = _mol_to_ase_atoms(ethane_mol)
        original_pos = atoms.get_positions().copy()
        # Shift positions slightly
        atoms.positions += 0.1
        _update_mol_from_ase(ethane_mol, atoms, conf_id=0)
        conf = ethane_mol.GetConformer(0)
        new_pos = conf.GetAtomPosition(0)
        # Should reflect the shifted coordinates
        assert abs(new_pos.x - (original_pos[0, 0] + 0.1)) < 0.01


# -------------------------------------------------------------------
# Mock-based tests (always run, no MACE required)
# -------------------------------------------------------------------

class TestMACEMockPaths:
    """Tests using mocks to verify MACE pipeline wiring without installation."""

    def test_screening_mace_path_mocked(self):
        """Mock MACE optimizer to verify screening.py calls it correctly."""
        mock_mol = Chem.MolFromSmiles("CCCCCC")
        mock_mol = Chem.AddHs(mock_mol)
        AllChem.EmbedMultipleConfs(mock_mol, numConfs=3, randomSeed=42)

        with patch("pyrene_analyzer.screening.has_mace_available", return_value=True), \
             patch("pyrene_analyzer.screening.optimize_conformer_ensemble_mace",
                   return_value=mock_mol) as mock_opt:
            from pyrene_analyzer.screening import generate_conformers
            result = generate_conformers(
                Chem.RWMol(mock_mol), num_confs=5, optimizer="MACE-OFF23",
            )
            mock_opt.assert_called_once()

    def test_screening_mace_fallback_path_mocked(self):
        """Mock MACE failure to verify MMFF94s fallback is called."""
        mol = Chem.MolFromSmiles("CCCCCC")
        mol = Chem.AddHs(mol)

        with patch("pyrene_analyzer.screening.has_mace_available", return_value=True), \
             patch("pyrene_analyzer.screening.optimize_conformer_ensemble_mace",
                   side_effect=RuntimeError("MACE exploded")):
            from pyrene_analyzer.screening import generate_conformers
            result = generate_conformers(
                mol, num_confs=5, optimizer="MACE-OFF23",
            )
        # Should still work via MMFF94s fallback
        assert result.GetNumConformers() > 0

    def test_mace_optimizer_module_always_importable(self):
        """mace_optimizer module must be importable even without MACE."""
        import pyrene_analyzer.mace_optimizer as mod
        # All public functions should exist
        assert callable(mod.has_mace_available)
        assert callable(mod.optimize_with_mace)
        assert callable(mod.optimize_conformer_ensemble_mace)
        assert callable(mod.filter_by_energy_mace)
        assert callable(mod.single_point_mace)
        assert callable(mod.optimize_or_fallback_mace)
        assert callable(mod.get_mace_import_error)

    def test_lazy_import_from_init(self):
        """MACE functions should be accessible via pyrene_analyzer.__init__."""
        import pyrene_analyzer
        assert hasattr(pyrene_analyzer, "has_mace_available")
        assert hasattr(pyrene_analyzer, "optimize_with_mace")
        assert hasattr(pyrene_analyzer, "optimize_conformer_ensemble_mace")
