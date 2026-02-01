"""
Tests for the I/O module.

Tests cover:
- SDF file loading
- MOL2 file loading
- PDB file loading
- CSV export
- JSON export
- Excel export
"""

import json
from pathlib import Path

import pandas as pd
import pytest

from pyrene_analyzer.io import (
    export_results,
    export_to_csv,
    export_to_excel,
    export_to_json,
    load_from_mol2,
    load_from_pdb,
    load_from_sdf,
    load_molecules,
)


class TestLoadFromSDF:
    """Tests for load_from_sdf function."""

    def test_load_simple_sdf(self, simple_test_sdf):
        """Test loading a simple SDF file."""
        molecules = load_from_sdf(simple_test_sdf)

        assert len(molecules) > 0
        for mol, name in molecules:
            assert mol is not None
            assert isinstance(name, str)

    def test_load_nonexistent_file(self):
        """Should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            load_from_sdf("nonexistent.sdf")

    def test_load_with_remove_hydrogens(self, simple_test_sdf):
        """Test loading with hydrogen removal."""
        molecules_with_h = load_from_sdf(simple_test_sdf, remove_hydrogens=False)
        molecules_no_h = load_from_sdf(simple_test_sdf, remove_hydrogens=True)

        # Molecules without H should have fewer atoms
        mol_with_h, _ = molecules_with_h[0]
        mol_no_h, _ = molecules_no_h[0]

        assert mol_no_h.GetNumAtoms() <= mol_with_h.GetNumAtoms()

    def test_molecule_names(self, simple_test_sdf):
        """Test that molecule names are extracted."""
        molecules = load_from_sdf(simple_test_sdf)

        for mol, name in molecules:
            assert name is not None
            assert len(name) > 0

    def test_real_sdf_file(self, test_sdf_path):
        """Test loading real SDF file if available."""
        if test_sdf_path is None:
            pytest.skip("Test SDF file not available")

        molecules = load_from_sdf(test_sdf_path)
        assert len(molecules) > 0


class TestLoadMolecules:
    """Tests for load_molecules function."""

    def test_auto_detect_sdf(self, simple_test_sdf):
        """Test automatic format detection for SDF."""
        molecules = load_molecules(simple_test_sdf)
        assert len(molecules) > 0

    def test_unsupported_format(self, tmp_path):
        """Should raise ValueError for unsupported format."""
        fake_file = tmp_path / "test.xyz"
        fake_file.write_text("fake content")

        with pytest.raises(ValueError):
            load_molecules(fake_file)

    def test_nonexistent_file(self):
        """Should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            load_molecules("nonexistent.sdf")


class TestExportToCSV:
    """Tests for export_to_csv function."""

    def test_basic_export(self, sample_analysis_df, tmp_path):
        """Test basic CSV export."""
        output_path = tmp_path / "test.csv"
        export_to_csv(sample_analysis_df, output_path)

        assert output_path.exists()

        # Read back and verify
        df_read = pd.read_csv(output_path)
        assert len(df_read) == len(sample_analysis_df)
        assert list(df_read.columns) == list(sample_analysis_df.columns)

    def test_float_format(self, sample_analysis_df, tmp_path):
        """Test float formatting in CSV."""
        output_path = tmp_path / "test.csv"
        export_to_csv(sample_analysis_df, output_path, float_format="%.2f")

        # Read back and check formatting
        with open(output_path) as f:
            content = f.read()
            # Should not have excessive decimal places
            assert ".000000" not in content

    def test_create_parent_dirs(self, sample_analysis_df, tmp_path):
        """Test that parent directories are created."""
        output_path = tmp_path / "subdir" / "nested" / "test.csv"
        export_to_csv(sample_analysis_df, output_path)

        assert output_path.exists()


class TestExportToJSON:
    """Tests for export_to_json function."""

    def test_basic_export(self, sample_analysis_df, tmp_path):
        """Test basic JSON export."""
        output_path = tmp_path / "test.json"
        export_to_json(sample_analysis_df, output_path)

        assert output_path.exists()

        # Read back and verify
        with open(output_path) as f:
            data = json.load(f)

        assert len(data) == len(sample_analysis_df)

    def test_json_structure(self, sample_analysis_df, tmp_path):
        """Test JSON structure with records orientation."""
        output_path = tmp_path / "test.json"
        export_to_json(sample_analysis_df, output_path, orient="records")

        with open(output_path) as f:
            data = json.load(f)

        # Should be a list of dictionaries
        assert isinstance(data, list)
        assert isinstance(data[0], dict)

    def test_nan_handling(self, tmp_path):
        """Test that NaN values are handled properly."""
        import numpy as np

        df = pd.DataFrame(
            {
                "a": [1.0, np.nan, 3.0],
                "b": ["x", "y", "z"],
            }
        )

        output_path = tmp_path / "test.json"
        export_to_json(df, output_path)

        with open(output_path) as f:
            data = json.load(f)

        # NaN should be converted to None/null
        assert data[1]["a"] is None


class TestExportToExcel:
    """Tests for export_to_excel function."""

    def test_basic_export(self, sample_analysis_df, tmp_path):
        """Test basic Excel export."""
        output_path = tmp_path / "test.xlsx"
        export_to_excel(sample_analysis_df, output_path)

        assert output_path.exists()

        # Read back and verify
        df_read = pd.read_excel(output_path)
        assert len(df_read) == len(sample_analysis_df)

    def test_auto_extension(self, sample_analysis_df, tmp_path):
        """Test that .xlsx extension is added automatically."""
        output_path = tmp_path / "test"
        export_to_excel(sample_analysis_df, output_path)

        assert (tmp_path / "test.xlsx").exists()

    def test_custom_sheet_name(self, sample_analysis_df, tmp_path):
        """Test custom sheet name."""
        output_path = tmp_path / "test.xlsx"
        export_to_excel(sample_analysis_df, output_path, sheet_name="Results")

        # Read back with sheet name
        df_read = pd.read_excel(output_path, sheet_name="Results")
        assert len(df_read) == len(sample_analysis_df)


class TestExportResults:
    """Tests for export_results function."""

    def test_export_all_formats(self, sample_analysis_df, tmp_path):
        """Test exporting to all formats."""
        output_path = tmp_path / "results"
        files = export_results(sample_analysis_df, output_path)

        assert len(files) == 3
        assert (tmp_path / "results.csv").exists()
        assert (tmp_path / "results.json").exists()
        assert (tmp_path / "results.xlsx").exists()

    def test_export_specific_formats(self, sample_analysis_df, tmp_path):
        """Test exporting to specific formats."""
        output_path = tmp_path / "results"
        files = export_results(sample_analysis_df, output_path, formats=["csv", "json"])

        assert len(files) == 2
        assert (tmp_path / "results.csv").exists()
        assert (tmp_path / "results.json").exists()
        assert not (tmp_path / "results.xlsx").exists()

    def test_return_paths(self, sample_analysis_df, tmp_path):
        """Test that correct paths are returned."""
        output_path = tmp_path / "results"
        files = export_results(sample_analysis_df, output_path, formats=["csv"])

        assert len(files) == 1
        assert files[0] == tmp_path / "results.csv"
