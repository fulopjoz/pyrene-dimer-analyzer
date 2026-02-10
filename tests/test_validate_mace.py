"""Regression tests for validate_mace.py benchmark wiring."""

from unittest.mock import MagicMock, Mock, patch

import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

import validate_mace as vm


def _make_test_mol() -> Chem.Mol:
    """Create a small molecule with one 3D conformer."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def test_optimize_with_xtb_method_calls_optimize_with_gfn2xtb():
    """xTB wrapper must call the current optimize_with_gfn2xtb API."""
    mol = _make_test_mol()
    with patch(
        "pyrene_analyzer.xtb_optimizer.optimize_with_gfn2xtb",
        return_value=(mol, -12.34),
    ) as mock_opt:
        out_mol, energy, elapsed = vm.optimize_with_xtb_method(
            Chem.RWMol(mol), max_steps=123, fmax=0.006
        )

    mock_opt.assert_called_once()
    _, kwargs = mock_opt.call_args
    assert kwargs["method"] == "GFN2-xTB"
    assert kwargs["max_steps"] == 123
    assert kwargs["fmax"] == pytest.approx(0.006)
    assert out_mol.GetNumAtoms() == mol.GetNumAtoms()
    assert energy == pytest.approx(-12.34)
    assert elapsed >= 0.0


def test_optimize_with_mmff94s_uses_interfragment_settings():
    """MMFF benchmark must include inter-fragment nonbonded interactions."""
    mol = _make_test_mol()
    fake_props = object()
    fake_ff = MagicMock()
    fake_ff.CalcEnergy.return_value = 12.5

    with patch.object(
        vm.AllChem, "MMFFGetMoleculeProperties", return_value=fake_props
    ) as mock_props, patch.object(
        vm.AllChem, "MMFFGetMoleculeForceField", return_value=fake_ff
    ) as mock_ff:
        out_mol, energy, _ = vm.optimize_with_mmff94s(
            Chem.RWMol(mol), max_steps=321, verbose=False
        )

    mock_props.assert_called_once()
    _, props_kwargs = mock_props.call_args
    assert props_kwargs["mmffVariant"] == "MMFF94s"

    mock_ff.assert_called_once()
    _, ff_kwargs = mock_ff.call_args
    assert ff_kwargs["confId"] == -1
    assert ff_kwargs["ignoreInterfragInteractions"] is False

    fake_ff.Minimize.assert_called_once_with(maxIts=321)
    assert out_mol.GetNumAtoms() == mol.GetNumAtoms()
    assert energy == pytest.approx(12.5)


def test_run_benchmark_uses_four_angstrom_initialization():
    """Benchmark should seed dimers at 4.0 A to avoid zero-force plateaus."""
    mol = _make_test_mol()
    fake_benchmarks = [
        {
            "name": "dummy",
            "smiles": "CC",
            "reference_distance_A": 3.4,
            "reference_source": "test",
            "n_aromatic_atoms": 2,
        }
    ]

    build_mock = Mock(return_value=(mol, 1))
    mace_mock = Mock(return_value=(mol, 0.0, 0.0))

    with patch.object(vm, "BENCHMARKS", fake_benchmarks), patch(
        "pyrene_analyzer.mace_optimizer.has_mace_available", return_value=True
    ), patch.object(vm, "build_stacked_dimer", build_mock), patch.object(
        vm, "measure_interplane_distance", return_value=4.0
    ), patch.object(
        vm, "optimize_with_mmff94s", return_value=(mol, 0.0, 0.0)
    ), patch.object(
        vm, "optimize_with_mace_method", mace_mock
    ):
        df = vm.run_benchmark(all_methods=False, verbose=False)

    build_mock.assert_called_once_with("CC", initial_distance=4.0)
    assert set(df["method"]) == {"MMFF94s", "MACE-OFF23", "MACE-OFF23-medium"}
    assert mace_mock.call_count == 2
    for call in mace_mock.call_args_list:
        assert call.kwargs["max_steps"] == 500
        assert call.kwargs["fmax"] == pytest.approx(0.005)


def test_evaluate_results_threshold_behavior():
    """Pass/fail criterion should remain |error| < 0.2 A for pyrene MACE."""
    df_pass = pd.DataFrame(
        [{"method": "MACE-OFF23", "benchmark": "pyrene_dimer", "error_A": 0.19}]
    )
    df_fail = pd.DataFrame(
        [{"method": "MACE-OFF23", "benchmark": "pyrene_dimer", "error_A": 0.21}]
    )

    assert bool(vm.evaluate_results(df_pass, verbose=False))
    assert not bool(vm.evaluate_results(df_fail, verbose=False))
