"""
Pytest configuration and fixtures for pyrene-dimer-analyzer tests.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

# Test data directory
TEST_DATA_DIR = Path(__file__).parent / "test_data"


@pytest.fixture
def sample_coords_parallel():
    """Two parallel planes at z=0 and z=3.5."""
    coords1 = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 0.0],
        ]
    )
    coords2 = np.array(
        [
            [0.0, 0.0, 3.5],
            [1.0, 0.0, 3.5],
            [1.0, 1.0, 3.5],
            [0.0, 1.0, 3.5],
            [0.5, 0.5, 3.5],
        ]
    )
    return coords1, coords2


@pytest.fixture
def sample_coords_perpendicular():
    """Two perpendicular planes."""
    coords1 = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    coords2 = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [0.0, 1.0, 0.0],
        ]
    )
    return coords1, coords2


@pytest.fixture
def sample_coords_45deg():
    """Two planes at approximately 45 degrees."""
    coords1 = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    # Rotate around x-axis by 45 degrees
    angle = np.radians(45)
    coords2 = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, np.cos(angle), np.sin(angle)],
            [0.0, np.cos(angle), np.sin(angle)],
        ]
    )
    return coords1, coords2


@pytest.fixture
def sample_coords_offset():
    """Two parallel planes with lateral offset."""
    coords1 = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    coords2 = np.array(
        [
            [2.0, 0.0, 3.5],
            [3.0, 0.0, 3.5],
            [3.0, 1.0, 3.5],
            [2.0, 1.0, 3.5],
        ]
    )
    return coords1, coords2


@pytest.fixture
def benzene_mol():
    """Create a simple benzene molecule."""
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


@pytest.fixture
def naphthalene_mol():
    """Create a naphthalene molecule."""
    mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


@pytest.fixture
def pyrene_mol():
    """Create a single pyrene molecule."""
    # Pyrene SMILES
    mol = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


@pytest.fixture
def simple_dimer_mol():
    """Create a simple biphenyl-like dimer for testing."""
    # Biphenyl as a simple test case
    mol = Chem.MolFromSmiles("c1ccc(-c2ccccc2)cc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


@pytest.fixture
def sample_analysis_df():
    """Create a sample analysis DataFrame."""
    np.random.seed(42)
    n = 50
    return pd.DataFrame(
        {
            "molecule": ["test_mol"] * n,
            "conformer_id": range(n),
            "plane_angle_deg": np.random.uniform(0, 90, n),
            "interplane_distance_A": np.random.uniform(3.0, 6.0, n),
            "pi_overlap_pct": np.random.uniform(0, 100, n),
            "centroid_distance_A": np.random.uniform(3.5, 8.0, n),
            "slip_stack_A": np.random.uniform(0, 3.0, n),
            "bridge_dihedral_L_deg": np.random.uniform(-180, 180, n),
            "bridge_dihedral_R_deg": np.random.uniform(-180, 180, n),
            "energy_kcal_mol": np.random.uniform(0, 20, n),
            "rel_energy_kcal_mol": np.random.uniform(0, 15, n),
        }
    )


@pytest.fixture
def multi_variant_df():
    """Create a DataFrame with multiple variants."""
    np.random.seed(42)
    dfs = []
    for variant in ["Et", "iPr", "cHex", "tBu"]:
        n = 20
        df = pd.DataFrame(
            {
                "molecule": [variant] * n,
                "conformer_id": range(n),
                "plane_angle_deg": np.random.uniform(0, 90, n),
                "interplane_distance_A": np.random.uniform(3.0, 6.0, n),
                "pi_overlap_pct": np.random.uniform(0, 100, n),
            }
        )
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


@pytest.fixture
def test_sdf_path():
    """Path to test SDF file."""
    path = TEST_DATA_DIR / "pyrene_dimer_set_for_MOE.sdf"
    if path.exists():
        return path
    return None


@pytest.fixture
def temp_output_dir(tmp_path):
    """Temporary directory for output files."""
    return tmp_path


def create_test_sdf(output_path: Path, n_mols: int = 3) -> None:
    """Create a test SDF file with simple molecules."""
    writer = Chem.SDWriter(str(output_path))

    for i in range(n_mols):
        mol = Chem.MolFromSmiles("c1ccc(-c2ccccc2)cc1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42 + i)
        AllChem.MMFFOptimizeMolecule(mol)
        mol.SetProp("_Name", f"test_mol_{i}")
        mol.SetProp("energy", str(float(i)))
        writer.write(mol)

    writer.close()


@pytest.fixture
def simple_test_sdf(tmp_path):
    """Create a simple test SDF file."""
    sdf_path = tmp_path / "test.sdf"
    create_test_sdf(sdf_path)
    return sdf_path
