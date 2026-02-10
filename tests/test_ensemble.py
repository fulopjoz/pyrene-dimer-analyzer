"""
Tests for the ensemble feature engineering module.

Tests cover:
    - Distributional feature computation
    - Threshold-based feature computation
    - Boltzmann-weighted feature computation
    - Combined ensemble feature computation
    - Edge cases (empty DataFrames, missing columns, NaN energies)
"""

import warnings

import numpy as np
import pandas as pd
import pytest

from pyrene_analyzer.ensemble import (
    GEOMETRIC_DESCRIPTORS,
    KT_298K,
    compute_boltzmann_weighted_features,
    compute_distributional_features,
    compute_ensemble_features,
    compute_threshold_features,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_conformer_df():
    """Single-molecule conformer ensemble with 50 conformers."""
    np.random.seed(42)
    n = 50
    angles = np.random.uniform(0, 90, n)
    distances = np.random.uniform(3.0, 6.0, n)
    overlaps = np.random.uniform(0, 100, n)

    # Assign classification based on simple thresholds
    classifications = []
    for a, d, o in zip(angles, distances, overlaps):
        if a < 20 and 3.3 <= d <= 3.7 and o > 50:
            classifications.append("strong_excimer")
        elif a < 60 and d < 4.5 and o > 30:
            classifications.append("weak_excimer")
        else:
            classifications.append("monomer")

    return pd.DataFrame(
        {
            "molecule": ["mol_A"] * n,
            "conformer_id": range(n),
            "plane_angle_deg": angles,
            "interplane_distance_A": distances,
            "pi_overlap_pct": overlaps,
            "centroid_distance_A": distances + np.random.uniform(0, 2, n),
            "slip_stack_A": np.random.uniform(0, 3, n),
            "energy_kcal_mol": np.random.uniform(300, 320, n),
            "rel_energy_kcal_mol": np.random.uniform(0, 10, n),
            "classification": classifications,
            "excimer_score": np.random.uniform(0, 1, n),
        }
    )


@pytest.fixture
def multi_molecule_df(sample_conformer_df):
    """Two molecules, 50 conformers each."""
    df_a = sample_conformer_df.copy()
    df_b = sample_conformer_df.copy()
    df_b["molecule"] = "mol_B"
    # Shift energies for mol_B
    df_b["energy_kcal_mol"] += 5.0
    df_b["plane_angle_deg"] = np.random.uniform(0, 45, len(df_b))
    return pd.concat([df_a, df_b], ignore_index=True)


@pytest.fixture
def no_energy_df(sample_conformer_df):
    """Conformer DataFrame without energy column."""
    return sample_conformer_df.drop(columns=["energy_kcal_mol", "rel_energy_kcal_mol"])


@pytest.fixture
def no_classification_df(sample_conformer_df):
    """Conformer DataFrame without classification column."""
    return sample_conformer_df.drop(columns=["classification"])


# ---------------------------------------------------------------------------
# TestDistributionalFeatures
# ---------------------------------------------------------------------------


class TestDistributionalFeatures:
    """Tests for compute_distributional_features."""

    def test_basic_output_shape(self, sample_conformer_df):
        """Should produce one row per molecule with 45 columns (5 descriptors x 9 stats)."""
        result = compute_distributional_features(sample_conformer_df)
        assert len(result) == 1
        # 5 geometric descriptors x 9 stats = 45
        assert len(result.columns) == 45

    def test_column_naming(self, sample_conformer_df):
        """Columns should follow {descriptor}_{stat} pattern."""
        result = compute_distributional_features(sample_conformer_df)
        for desc in GEOMETRIC_DESCRIPTORS:
            for stat in ["mean", "std", "min", "max", "median", "p10", "p25", "p75", "p90"]:
                col = f"{desc}_{stat}"
                assert col in result.columns, f"Missing column: {col}"

    def test_mean_matches_pandas(self, sample_conformer_df):
        """Mean should match pandas groupby mean."""
        result = compute_distributional_features(sample_conformer_df)
        expected_mean = sample_conformer_df["plane_angle_deg"].mean()
        assert abs(result["plane_angle_deg_mean"].iloc[0] - expected_mean) < 1e-10

    def test_min_max_match(self, sample_conformer_df):
        """Min and max should match actual values."""
        result = compute_distributional_features(sample_conformer_df)
        assert result["pi_overlap_pct_min"].iloc[0] == pytest.approx(
            sample_conformer_df["pi_overlap_pct"].min()
        )
        assert result["pi_overlap_pct_max"].iloc[0] == pytest.approx(
            sample_conformer_df["pi_overlap_pct"].max()
        )

    def test_multi_molecule(self, multi_molecule_df):
        """Should produce one row per molecule."""
        result = compute_distributional_features(multi_molecule_df)
        assert len(result) == 2
        assert "mol_A" in result.index
        assert "mol_B" in result.index

    def test_empty_dataframe(self):
        """Empty DataFrame should return empty result."""
        result = compute_distributional_features(pd.DataFrame())
        assert result.empty

    def test_missing_descriptors(self):
        """DataFrame without geometric columns should warn and return index."""
        df = pd.DataFrame({"molecule": ["a", "b"], "other": [1, 2]})
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = compute_distributional_features(df)
            assert len(w) == 1
            assert "No geometric descriptor" in str(w[0].message)

    def test_custom_descriptors(self, sample_conformer_df):
        """Custom descriptors should be respected."""
        custom = ("plane_angle_deg",)
        result = compute_distributional_features(
            sample_conformer_df, descriptors=custom
        )
        assert len(result.columns) == 9  # 1 descriptor x 9 stats

    def test_percentiles_ordered(self, sample_conformer_df):
        """Percentiles should be ordered: p10 <= p25 <= median <= p75 <= p90."""
        result = compute_distributional_features(sample_conformer_df)
        for desc in GEOMETRIC_DESCRIPTORS:
            p10 = result[f"{desc}_p10"].iloc[0]
            p25 = result[f"{desc}_p25"].iloc[0]
            med = result[f"{desc}_median"].iloc[0]
            p75 = result[f"{desc}_p75"].iloc[0]
            p90 = result[f"{desc}_p90"].iloc[0]
            assert p10 <= p25 <= med <= p75 <= p90


# ---------------------------------------------------------------------------
# TestThresholdFeatures
# ---------------------------------------------------------------------------


class TestThresholdFeatures:
    """Tests for compute_threshold_features."""

    def test_basic_output(self, sample_conformer_df):
        """Should produce fractions and counts."""
        result = compute_threshold_features(sample_conformer_df)
        assert len(result) == 1
        assert "n_conformers" in result.columns
        assert "frac_any_excimer" in result.columns

    def test_n_conformers_correct(self, sample_conformer_df):
        """n_conformers should equal actual count."""
        result = compute_threshold_features(sample_conformer_df)
        assert result["n_conformers"].iloc[0] == 50

    def test_fractions_sum_to_one(self, sample_conformer_df):
        """Excimer + monomer fractions should sum to 1."""
        result = compute_threshold_features(sample_conformer_df)
        total = result["frac_any_excimer"].iloc[0] + result["frac_monomer"].iloc[0]
        assert total == pytest.approx(1.0)

    def test_strong_leq_any(self, sample_conformer_df):
        """frac_strong_excimer should be <= frac_any_excimer."""
        result = compute_threshold_features(sample_conformer_df)
        assert result["frac_strong_excimer"].iloc[0] <= result["frac_any_excimer"].iloc[0]

    def test_fractions_in_range(self, sample_conformer_df):
        """All fractions should be in [0, 1]."""
        result = compute_threshold_features(sample_conformer_df)
        frac_cols = [c for c in result.columns if c.startswith("frac_")]
        for col in frac_cols:
            val = result[col].iloc[0]
            assert 0.0 <= val <= 1.0, f"{col} = {val} out of range"

    def test_geometry_thresholds(self, sample_conformer_df):
        """Geometry-based thresholds should be present."""
        result = compute_threshold_features(sample_conformer_df)
        assert "frac_theta_lt_20" in result.columns
        assert "frac_d_in_range" in result.columns
        assert "frac_d_lt_4p5" in result.columns
        assert "frac_overlap_gt_50" in result.columns
        assert "frac_overlap_gt_30" in result.columns

    def test_no_classification_warns(self, no_classification_df):
        """Missing classification column should warn."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = compute_threshold_features(no_classification_df)
            assert any("classification" in str(x.message) for x in w)
            # Should still have geometry thresholds
            assert "frac_theta_lt_20" in result.columns

    def test_empty_dataframe(self):
        """Empty DataFrame should return empty result."""
        result = compute_threshold_features(pd.DataFrame())
        assert result.empty

    def test_custom_thresholds(self, sample_conformer_df):
        """Custom thresholds should change results."""
        # Very permissive thresholds
        loose = {
            "strong_angle_max": 90.0,
            "strong_distance_range": (0.0, 100.0),
            "strong_overlap_min": 0.0,
            "weak_distance_max": 100.0,
            "weak_overlap_min": 0.0,
        }
        result = compute_threshold_features(
            sample_conformer_df, thresholds=loose
        )
        # With weak_distance_max=100, all conformers have d < 100
        assert result["frac_d_lt_4p5"].iloc[0] == pytest.approx(1.0)
        # With strong_distance_range=(0, 100), all conformers are in range
        assert result["frac_d_in_range"].iloc[0] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# TestBoltzmannFeatures
# ---------------------------------------------------------------------------


class TestBoltzmannFeatures:
    """Tests for compute_boltzmann_weighted_features."""

    def test_basic_output(self, sample_conformer_df):
        """Should produce Boltzmann-weighted columns."""
        result = compute_boltzmann_weighted_features(sample_conformer_df)
        assert len(result) == 1
        assert "plane_angle_deg_boltz" in result.columns

    def test_boltz_columns_present(self, sample_conformer_df):
        """All descriptor _boltz columns should be present."""
        result = compute_boltzmann_weighted_features(sample_conformer_df)
        for desc in GEOMETRIC_DESCRIPTORS:
            assert f"{desc}_boltz" in result.columns

    def test_boltz_excimer_fraction(self, sample_conformer_df):
        """Should have boltz_excimer_fraction when classification present."""
        result = compute_boltzmann_weighted_features(sample_conformer_df)
        assert "boltz_excimer_fraction" in result.columns
        val = result["boltz_excimer_fraction"].iloc[0]
        assert 0.0 <= val <= 1.0

    def test_boltz_excimer_score(self, sample_conformer_df):
        """Should have boltz_excimer_score when excimer_score present."""
        result = compute_boltzmann_weighted_features(sample_conformer_df)
        assert "boltz_excimer_score" in result.columns

    def test_no_energy_warns(self, no_energy_df):
        """Missing energy column should warn and return empty-like DataFrame."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = compute_boltzmann_weighted_features(no_energy_df)
            assert any("Energy column" in str(x.message) for x in w)

    def test_empty_dataframe(self):
        """Empty DataFrame should return empty result."""
        result = compute_boltzmann_weighted_features(pd.DataFrame())
        assert result.empty

    def test_single_conformer(self):
        """Single conformer should be skipped (need >= 2)."""
        df = pd.DataFrame(
            {
                "molecule": ["A"],
                "energy_kcal_mol": [100.0],
                "plane_angle_deg": [15.0],
                "interplane_distance_A": [3.5],
                "pi_overlap_pct": [60.0],
                "centroid_distance_A": [4.0],
                "slip_stack_A": [1.0],
            }
        )
        result = compute_boltzmann_weighted_features(df)
        # Should return index but no computed values
        assert len(result) == 0 or result.isna().all().all()

    def test_nan_energies(self, sample_conformer_df):
        """NaN energies should be handled gracefully."""
        df = sample_conformer_df.copy()
        df.loc[0:5, "energy_kcal_mol"] = np.nan
        result = compute_boltzmann_weighted_features(df)
        assert len(result) == 1
        # Values should still be finite
        for col in result.columns:
            val = result[col].iloc[0]
            if not pd.isna(val):
                assert np.isfinite(val)

    def test_boltzmann_favors_low_energy(self):
        """Boltzmann weighting should favor low-energy conformers."""
        np.random.seed(99)
        # Create two conformers: low E with small angle, high E with large angle
        df = pd.DataFrame(
            {
                "molecule": ["X"] * 20,
                "energy_kcal_mol": [300.0] * 10 + [310.0] * 10,
                "plane_angle_deg": [10.0] * 10 + [80.0] * 10,
                "interplane_distance_A": [3.5] * 10 + [6.0] * 10,
                "pi_overlap_pct": [70.0] * 10 + [5.0] * 10,
                "centroid_distance_A": [4.0] * 20,
                "slip_stack_A": [1.0] * 20,
            }
        )
        result = compute_boltzmann_weighted_features(df)
        # Boltzmann-weighted angle should be closer to 10 than to 80
        boltz_angle = result["plane_angle_deg_boltz"].iloc[0]
        simple_mean = 45.0
        assert boltz_angle < simple_mean

    def test_kt_constant(self):
        """KT_298K should be approximately R*298."""
        expected = 0.001987204 * 298.0
        assert KT_298K == pytest.approx(expected, abs=0.01)


# ---------------------------------------------------------------------------
# TestComputeEnsembleFeatures
# ---------------------------------------------------------------------------


class TestComputeEnsembleFeatures:
    """Tests for compute_ensemble_features (main entry point)."""

    def test_combines_all_feature_sets(self, sample_conformer_df):
        """Should produce distributional + threshold + Boltzmann columns."""
        result = compute_ensemble_features(sample_conformer_df)
        assert len(result) == 1

        # Distributional columns
        assert "plane_angle_deg_mean" in result.columns
        assert "interplane_distance_A_std" in result.columns

        # Threshold columns
        assert "n_conformers" in result.columns
        assert "frac_any_excimer" in result.columns

        # Boltzmann columns
        assert "plane_angle_deg_boltz" in result.columns

    def test_total_column_count(self, sample_conformer_df):
        """Should have ~55+ feature columns."""
        result = compute_ensemble_features(sample_conformer_df)
        assert len(result.columns) >= 40

    def test_multi_molecule(self, multi_molecule_df):
        """Should produce one row per molecule."""
        result = compute_ensemble_features(multi_molecule_df)
        assert len(result) == 2

    def test_empty_dataframe(self):
        """Empty DataFrame should return empty result."""
        result = compute_ensemble_features(pd.DataFrame())
        assert result.empty

    def test_group_col_parameter(self, sample_conformer_df):
        """Custom group_col should work."""
        df = sample_conformer_df.copy()
        df["substituent"] = "Et"
        result = compute_ensemble_features(df, group_col="substituent")
        assert "Et" in result.index

    def test_no_energy_still_works(self, no_energy_df):
        """Should still produce distributional and threshold features without energy."""
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            result = compute_ensemble_features(no_energy_df)
            # Should have distributional and threshold columns
            assert "plane_angle_deg_mean" in result.columns
            assert "n_conformers" in result.columns
