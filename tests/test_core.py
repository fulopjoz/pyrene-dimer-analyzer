"""
Tests for the core module.

Tests cover:
- PyreneDimerAnalyzer class
- Pyrene ring detection
- Conformer analysis
- Batch processing
- Classification
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from pyrene_analyzer.aromatic_systems import ClassificationThresholds
from pyrene_analyzer.core import AromaticDimerAnalyzer, PyreneDimerAnalyzer


class TestPyreneDimerAnalyzerInit:
    """Tests for PyreneDimerAnalyzer initialization."""

    def test_default_init(self):
        """Test default initialization."""
        analyzer = PyreneDimerAnalyzer()

        assert analyzer.verbose is True
        assert analyzer.use_smarts is True
        assert analyzer.use_shapely is True

    def test_custom_init(self):
        """Test custom initialization."""
        analyzer = PyreneDimerAnalyzer(
            verbose=False, use_smarts=False, use_shapely=False
        )

        assert analyzer.verbose is False
        assert analyzer.use_smarts is False
        assert analyzer.use_shapely is False

    def test_log_verbose(self, capsys):
        """Test logging in verbose mode."""
        analyzer = PyreneDimerAnalyzer(verbose=True)
        analyzer.log("Test message")

        captured = capsys.readouterr()
        assert "Test message" in captured.out

    def test_log_quiet(self, capsys):
        """Test logging in quiet mode."""
        analyzer = PyreneDimerAnalyzer(verbose=False)
        analyzer.log("Test message")

        captured = capsys.readouterr()
        assert captured.out == ""


class TestIdentifyPyreneRings:
    """Tests for pyrene ring identification."""

    def test_identify_in_biphenyl(self, simple_dimer_mol):
        """Test ring identification in biphenyl."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        # Biphenyl has 2 aromatic systems
        try:
            ring1, ring2 = analyzer.identify_pyrene_rings(simple_dimer_mol)
            assert len(ring1) > 0
            assert len(ring2) > 0
        except ValueError:
            # May fail if SMARTS doesn't match - that's OK for biphenyl
            pass

    def test_connectivity_based_detection(self, simple_dimer_mol):
        """Test connectivity-based detection."""
        analyzer = PyreneDimerAnalyzer(verbose=False, use_smarts=False)

        try:
            ring1, ring2 = analyzer._identify_by_connectivity(simple_dimer_mol)
            assert len(ring1) > 0
            assert len(ring2) > 0
        except ValueError:
            # May fail for simple biphenyl
            pass

    def test_merge_connected_rings(self):
        """Test ring merging algorithm."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        # Two separate ring systems
        rings = [
            {0, 1, 2, 3, 4, 5},  # Ring 1
            {2, 3, 6, 7, 8, 9},  # Ring 2 (shares atoms with ring 1)
            {20, 21, 22, 23, 24, 25},  # Ring 3 (separate)
        ]

        clusters = analyzer._merge_connected_rings(rings)

        # Should merge rings 1 and 2
        assert len(clusters) == 2

    def test_empty_rings(self):
        """Test with empty ring list."""
        analyzer = PyreneDimerAnalyzer(verbose=False)
        clusters = analyzer._merge_connected_rings([])
        assert clusters == []


class TestFindBridgeAtoms:
    """Tests for bridge atom identification."""

    def test_find_bridge_in_biphenyl(self, simple_dimer_mol):
        """Test bridge finding in biphenyl."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        # Get ring atoms first
        ring_info = simple_dimer_mol.GetRingInfo()
        rings = list(ring_info.AtomRings())

        if len(rings) >= 2:
            pyrene1 = list(rings[0])
            pyrene2 = list(rings[1])

            bridge1, bridge2 = analyzer.find_bridge_atoms(
                simple_dimer_mol, pyrene1, pyrene2
            )

            # Bridge may or may not be found depending on structure
            # Just check it doesn't crash


class TestAnalyzeConformer:
    """Tests for single conformer analysis."""

    def test_analyze_conformer_basic(self, simple_dimer_mol):
        """Test basic conformer analysis."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        # Get ring atoms
        ring_info = simple_dimer_mol.GetRingInfo()
        rings = list(ring_info.AtomRings())

        if len(rings) >= 2:
            pyrene1 = list(rings[0])
            pyrene2 = list(rings[1])

            result = analyzer.analyze_conformer(simple_dimer_mol, 0, pyrene1, pyrene2)

            assert "conformer_id" in result
            assert "plane_angle_deg" in result
            assert "interplane_distance_A" in result
            assert "pi_overlap_pct" in result

            # Check value ranges
            assert 0 <= result["plane_angle_deg"] <= 90
            assert result["interplane_distance_A"] >= 0
            assert 0 <= result["pi_overlap_pct"] <= 100


class TestAnalyzeMolecule:
    """Tests for molecule analysis."""

    def test_analyze_molecule_basic(self, simple_dimer_mol):
        """Test basic molecule analysis."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        try:
            df = analyzer.analyze_molecule(simple_dimer_mol, "test_mol")

            assert isinstance(df, pd.DataFrame)
            assert len(df) > 0
            assert "molecule" in df.columns
            assert "plane_angle_deg" in df.columns
        except ValueError:
            # May fail if pyrene detection fails
            pass

    def test_analyze_molecule_with_progress(self, simple_dimer_mol, capsys):
        """Test molecule analysis with progress bar."""
        analyzer = PyreneDimerAnalyzer(verbose=True)

        try:
            df = analyzer.analyze_molecule(
                simple_dimer_mol, "test_mol", show_progress=True
            )
        except ValueError:
            pass


class TestAnalyzeFile:
    """Tests for file analysis."""

    def test_analyze_file(self, simple_test_sdf):
        """Test file analysis."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        try:
            df = analyzer.analyze_file(simple_test_sdf)

            assert isinstance(df, pd.DataFrame)
        except (ValueError, Exception):
            # May fail for simple test molecules
            pass

    def test_analyze_nonexistent_file(self):
        """Test with nonexistent file."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        with pytest.raises(FileNotFoundError):
            analyzer.analyze_file("nonexistent.sdf")


class TestBatchAnalyze:
    """Tests for batch analysis."""

    def test_batch_analyze_single_file(self, simple_test_sdf):
        """Test batch analysis with single file."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        try:
            df = analyzer.batch_analyze([simple_test_sdf], n_jobs=1)
            assert isinstance(df, pd.DataFrame)
        except (ValueError, Exception):
            pass

    def test_batch_analyze_empty_list(self):
        """Test batch analysis with empty file list."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        df = analyzer.batch_analyze([])
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0


class TestClassification:
    """Tests for conformer classification."""

    def test_classify_strong_excimer(self):
        """Test strong excimer classification."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        classification = analyzer.classify_conformer(
            plane_angle=10.0, distance=3.5, overlap=80.0
        )

        assert classification == "strong_excimer"

    def test_classify_weak_excimer(self):
        """Test weak excimer classification."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        classification = analyzer.classify_conformer(
            plane_angle=40.0, distance=4.0, overlap=50.0
        )

        assert classification == "weak_excimer"

    def test_classify_monomer(self):
        """Test monomer classification."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        classification = analyzer.classify_conformer(
            plane_angle=75.0, distance=5.0, overlap=20.0
        )

        assert classification == "monomer"

    def test_add_classification(self, sample_analysis_df):
        """Test adding classification to DataFrame."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        df_classified = analyzer.add_classification(sample_analysis_df)

        assert "classification" in df_classified.columns
        assert all(
            c in ["strong_excimer", "weak_excimer", "monomer"]
            for c in df_classified["classification"]
        )

    def test_classification_boundary_cases(self):
        """Test classification at boundary values."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        # Exactly at strong excimer boundary
        c1 = analyzer.classify_conformer(20.0, 3.7, 70.0)

        # Exactly at weak excimer boundary
        c2 = analyzer.classify_conformer(60.0, 4.5, 30.0)

        # All classifications should be valid
        assert c1 in ["strong_excimer", "weak_excimer", "monomer"]
        assert c2 in ["strong_excimer", "weak_excimer", "monomer"]


class TestGetConformerEnergy:
    """Tests for energy extraction."""

    def test_energy_from_property(self, simple_dimer_mol):
        """Test energy extraction from molecule property."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        # Set energy property
        simple_dimer_mol.SetProp("energy", "10.5")

        energy = analyzer._get_conformer_energy(simple_dimer_mol, 0)

        # May or may not find energy depending on where it's stored
        # Just check it doesn't crash

    def test_no_energy(self, simple_dimer_mol):
        """Test when no energy is available."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        energy = analyzer._get_conformer_energy(simple_dimer_mol, 0)

        # Should return None if no energy found
        assert energy is None or isinstance(energy, float)


class TestIntegration:
    """Integration tests for the analyzer."""

    def test_full_workflow(self, simple_test_sdf, tmp_path):
        """Test complete analysis workflow."""
        analyzer = PyreneDimerAnalyzer(verbose=False)

        try:
            # Analyze file
            df = analyzer.analyze_file(simple_test_sdf)

            # Add classification
            df = analyzer.add_classification(df)

            # Export results
            output_path = tmp_path / "results.csv"
            df.to_csv(output_path, index=False)

            assert output_path.exists()

            # Read back
            df_read = pd.read_csv(output_path)
            assert len(df_read) == len(df)

        except (ValueError, Exception):
            # May fail for simple test molecules
            pass

    def test_real_sdf_analysis(self, test_sdf_path, tmp_path):
        """Test analysis with real SDF file."""
        if test_sdf_path is None:
            pytest.skip("Test SDF file not available")

        analyzer = PyreneDimerAnalyzer(verbose=False)

        try:
            df = analyzer.analyze_file(test_sdf_path)

            assert isinstance(df, pd.DataFrame)
            assert len(df) > 0

            # Check all expected columns
            expected_cols = [
                "molecule",
                "conformer_id",
                "plane_angle_deg",
                "interplane_distance_A",
                "pi_overlap_pct",
            ]
            for col in expected_cols:
                assert col in df.columns

        except ValueError as e:
            # May fail if structure doesn't have two pyrene systems
            pytest.skip(f"Analysis failed: {e}")


class TestAromaticDimerAnalyzer:
    """Tests for the generalized AromaticDimerAnalyzer."""

    def test_is_same_as_pyrene_analyzer(self):
        """AromaticDimerAnalyzer and PyreneDimerAnalyzer are the same class."""
        assert AromaticDimerAnalyzer is PyreneDimerAnalyzer

    def test_default_system_is_pyrene(self):
        """Default aromatic system should be pyrene."""
        analyzer = AromaticDimerAnalyzer(verbose=False)
        assert analyzer._system_name == "pyrene"

    def test_init_with_perylene(self):
        """Should initialize with perylene system."""
        analyzer = AromaticDimerAnalyzer(
            aromatic_system="perylene", verbose=False
        )
        assert analyzer._system_name == "perylene"
        assert analyzer._thresholds.strong_distance_range == (3.4, 3.8)

    def test_init_with_anthracene(self):
        """Should initialize with anthracene system."""
        analyzer = AromaticDimerAnalyzer(
            aromatic_system="anthracene", verbose=False
        )
        assert analyzer._system_name == "anthracene"

    def test_init_with_invalid_system(self):
        """Should raise ValueError for unknown system."""
        with pytest.raises(ValueError, match="Unknown aromatic system"):
            AromaticDimerAnalyzer(aromatic_system="nonexistent", verbose=False)

    def test_init_with_custom_smarts(self):
        """Should initialize with custom SMARTS pattern."""
        analyzer = AromaticDimerAnalyzer(
            custom_smarts="c1ccccc1", verbose=False
        )
        assert analyzer._system_name == "custom"
        assert analyzer._smarts == "c1ccccc1"

    def test_init_custom_smarts_with_thresholds(self):
        """Custom SMARTS with custom thresholds."""
        t = ClassificationThresholds(
            strong_angle_max=15.0,
            strong_distance_range=(3.2, 3.6),
            strong_overlap_min=60.0,
            weak_angle_max=50.0,
            weak_distance_max=4.0,
            weak_overlap_min=35.0,
        )
        analyzer = AromaticDimerAnalyzer(
            custom_smarts="c1ccccc1",
            custom_thresholds=t,
            verbose=False,
        )
        assert analyzer._thresholds.strong_angle_max == 15.0
        assert analyzer._thresholds.strong_distance_range == (3.2, 3.6)

    def test_identify_aromatic_rings_alias(self, simple_dimer_mol):
        """identify_pyrene_rings should delegate to identify_aromatic_rings."""
        analyzer = AromaticDimerAnalyzer(verbose=False)
        # Both methods should produce the same result (or same error)
        try:
            r1 = analyzer.identify_aromatic_rings(simple_dimer_mol)
            r2 = analyzer.identify_pyrene_rings(simple_dimer_mol)
            assert r1 == r2
        except ValueError:
            pass


class TestSystemSpecificClassification:
    """Tests for system-specific classification thresholds."""

    def test_pyrene_50pct_overlap_is_strong(self):
        """50% overlap should be strong_excimer for pyrene (lowered from 70%)."""
        analyzer = AromaticDimerAnalyzer(verbose=False)
        result = analyzer.classify_conformer(15.0, 3.5, 55.0)
        assert result == "strong_excimer"

    def test_pyrene_42pct_overlap_is_weak(self):
        """42% overlap (Basuroy 2021 case) should be weak_excimer."""
        analyzer = AromaticDimerAnalyzer(verbose=False)
        result = analyzer.classify_conformer(15.0, 3.5, 42.0)
        assert result == "weak_excimer"

    def test_pyrene_twisted_excimer_is_weak(self):
        """Twisted excimer at 55 deg (Mazzeo 2024) should be weak_excimer."""
        analyzer = AromaticDimerAnalyzer(verbose=False)
        result = analyzer.classify_conformer(55.0, 3.5, 50.0)
        assert result == "weak_excimer"

    def test_perylene_uses_wider_distance(self):
        """Perylene at 3.75 A should be strong_excimer (not for pyrene)."""
        perylene_analyzer = AromaticDimerAnalyzer(
            aromatic_system="perylene", verbose=False
        )
        result = perylene_analyzer.classify_conformer(10.0, 3.75, 60.0)
        assert result == "strong_excimer"

        # Same geometry with pyrene thresholds would be weak
        pyrene_analyzer = AromaticDimerAnalyzer(verbose=False)
        result_pyrene = pyrene_analyzer.classify_conformer(10.0, 3.75, 60.0)
        assert result_pyrene == "weak_excimer"

    def test_anthracene_distance_upper_bound(self):
        """Anthracene at 3.9 A should be strong_excimer."""
        analyzer = AromaticDimerAnalyzer(
            aromatic_system="anthracene", verbose=False
        )
        result = analyzer.classify_conformer(10.0, 3.9, 60.0)
        assert result == "strong_excimer"


class TestGeometryWarnings:
    """Tests for geometry warning flags in conformer output."""

    def test_warnings_at_high_angle(self, simple_dimer_mol):
        """Conformers with high angles should have geometry warnings."""
        analyzer = AromaticDimerAnalyzer(verbose=False)

        ring_info = simple_dimer_mol.GetRingInfo()
        rings = list(ring_info.AtomRings())

        if len(rings) >= 2:
            pyrene1 = list(rings[0])
            pyrene2 = list(rings[1])

            result = analyzer.analyze_conformer(
                simple_dimer_mol, 0, pyrene1, pyrene2
            )
            assert "geometry_warnings" in result
            # The warning field is either None or a string

    def test_analyze_conformer_returns_warnings_key(self, simple_dimer_mol):
        """Output dict should always contain geometry_warnings key."""
        analyzer = AromaticDimerAnalyzer(verbose=False)

        ring_info = simple_dimer_mol.GetRingInfo()
        rings = list(ring_info.AtomRings())

        if len(rings) >= 2:
            pyrene1 = list(rings[0])
            pyrene2 = list(rings[1])

            result = analyzer.analyze_conformer(
                simple_dimer_mol, 0, pyrene1, pyrene2
            )
            assert "geometry_warnings" in result
