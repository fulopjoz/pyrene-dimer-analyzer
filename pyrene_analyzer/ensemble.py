"""
Ensemble Feature Engineering Module
====================================

Statistical feature extraction from conformer ensemble analysis results.

Takes per-conformer DataFrames (from analyze_file or analyze_from_smiles)
and computes per-molecule summary features including:

- Distributional statistics (mean, std, min, max, percentiles)
- Threshold-based features (fraction excimer, fraction in range)
- Boltzmann-weighted ensemble averages at 298K

These features are suitable for downstream ML, SAR analysis, or
comparative ranking of dimer substituents.

Example::

    >>> from pyrene_analyzer.ensemble import compute_ensemble_features
    >>> results = analyzer.analyze_file("conformers.sdf")
    >>> results = analyzer.add_classification(results)
    >>> features = compute_ensemble_features(results)
"""

import warnings
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd


# Boltzmann constant * temperature at 298K in kcal/mol
KT_298K: float = 0.593

# Geometric descriptor columns produced by AromaticDimerAnalyzer
GEOMETRIC_DESCRIPTORS: Tuple[str, ...] = (
    "plane_angle_deg",
    "interplane_distance_A",
    "pi_overlap_pct",
    "centroid_distance_A",
    "slip_stack_A",
)


def compute_distributional_features(
    df: pd.DataFrame,
    group_col: str = "molecule",
    descriptors: Optional[Tuple[str, ...]] = None,
) -> pd.DataFrame:
    """Compute distributional statistics for each molecule in the ensemble.

    For each geometric descriptor, computes: mean, std, min, max, median,
    p10, p25, p75, p90.

    Args:
        df: Per-conformer analysis DataFrame.
        group_col: Column name to group by (default: ``"molecule"``).
        descriptors: Tuple of column names to compute stats for.
            If None, uses :data:`GEOMETRIC_DESCRIPTORS`.

    Returns:
        DataFrame with one row per molecule and distributional feature columns.
        Column naming: ``{descriptor}_{stat}``.
    """
    if df.empty:
        return pd.DataFrame()

    if descriptors is None:
        descriptors = GEOMETRIC_DESCRIPTORS

    # Filter to descriptors that actually exist in the DataFrame
    available = [d for d in descriptors if d in df.columns]
    if not available:
        warnings.warn(
            "No geometric descriptor columns found in DataFrame.",
            UserWarning,
            stacklevel=2,
        )
        return pd.DataFrame(index=df[group_col].unique())

    agg_funcs = {
        "mean": "mean",
        "std": "std",
        "min": "min",
        "max": "max",
        "median": "median",
        "p10": lambda x: np.nanpercentile(x.dropna(), 10) if len(x.dropna()) > 0 else np.nan,
        "p25": lambda x: np.nanpercentile(x.dropna(), 25) if len(x.dropna()) > 0 else np.nan,
        "p75": lambda x: np.nanpercentile(x.dropna(), 75) if len(x.dropna()) > 0 else np.nan,
        "p90": lambda x: np.nanpercentile(x.dropna(), 90) if len(x.dropna()) > 0 else np.nan,
    }

    result_frames = []
    grouped = df.groupby(group_col)

    for desc in available:
        for stat_name, func in agg_funcs.items():
            col_name = f"{desc}_{stat_name}"
            if callable(func):
                series = grouped[desc].agg(func)
            else:
                series = grouped[desc].agg(func)
            series.name = col_name
            result_frames.append(series)

    result = pd.concat(result_frames, axis=1)
    result.index.name = group_col
    return result


def compute_threshold_features(
    df: pd.DataFrame,
    group_col: str = "molecule",
    thresholds: Optional[Dict] = None,
) -> pd.DataFrame:
    """Compute threshold-based ensemble features.

    Calculates the fraction of conformers meeting specific geometric
    criteria relevant to excimer formation.

    Args:
        df: Per-conformer analysis DataFrame. Should contain a
            ``"classification"`` column for excimer fraction features.
        group_col: Column to group by.
        thresholds: Optional dict to override default thresholds. Keys:
            ``"strong_distance_range"`` (tuple), ``"strong_angle_max"`` (float),
            ``"weak_distance_max"`` (float). Defaults to pyrene values.

    Returns:
        DataFrame with one row per molecule and threshold feature columns.
    """
    if df.empty:
        return pd.DataFrame()

    if thresholds is None:
        thresholds = {
            "strong_angle_max": 20.0,
            "strong_distance_range": (3.3, 3.7),
            "strong_overlap_min": 50.0,
            "weak_distance_max": 4.5,
            "weak_overlap_min": 30.0,
        }

    grouped = df.groupby(group_col)
    result = pd.DataFrame(index=grouped.groups.keys())
    result.index.name = group_col

    result["n_conformers"] = grouped.size()

    # Classification-based features
    has_classification = "classification" in df.columns
    if has_classification:
        result["frac_strong_excimer"] = grouped["classification"].apply(
            lambda x: (x == "strong_excimer").mean()
        )
        result["frac_weak_excimer"] = grouped["classification"].apply(
            lambda x: (x == "weak_excimer").mean()
        )
        result["frac_any_excimer"] = grouped["classification"].apply(
            lambda x: x.isin(["strong_excimer", "weak_excimer"]).mean()
        )
        result["frac_monomer"] = grouped["classification"].apply(
            lambda x: (x == "monomer").mean()
        )
    else:
        warnings.warn(
            "No 'classification' column found. Call "
            "analyzer.add_classification() first for excimer fraction features.",
            UserWarning,
            stacklevel=2,
        )

    # Geometry-based threshold features
    if "plane_angle_deg" in df.columns:
        result["frac_theta_lt_20"] = grouped["plane_angle_deg"].apply(
            lambda x: (x < thresholds["strong_angle_max"]).mean()
        )

    if "interplane_distance_A" in df.columns:
        d_min, d_max = thresholds["strong_distance_range"]
        result["frac_d_in_range"] = grouped["interplane_distance_A"].apply(
            lambda x: ((x >= d_min) & (x <= d_max)).mean()
        )
        result["frac_d_lt_4p5"] = grouped["interplane_distance_A"].apply(
            lambda x: (x < thresholds["weak_distance_max"]).mean()
        )

    if "pi_overlap_pct" in df.columns:
        result["frac_overlap_gt_50"] = grouped["pi_overlap_pct"].apply(
            lambda x: (x > thresholds["strong_overlap_min"]).mean()
        )
        result["frac_overlap_gt_30"] = grouped["pi_overlap_pct"].apply(
            lambda x: (x > thresholds["weak_overlap_min"]).mean()
        )

    # Dark excimer fraction: eclipsed conformers (θ < 5°) that are classified
    # as excimers but may be optically dark (Dai et al. 2024 Molecules 29, 507)
    if "dark_excimer_warning" in df.columns:
        result["frac_dark_excimer"] = grouped["dark_excimer_warning"].apply(
            lambda x: x.mean()
        )

    return result


def compute_boltzmann_weighted_features(
    df: pd.DataFrame,
    group_col: str = "molecule",
    energy_col: str = "energy_kcal_mol",
    temperature_K: float = 298.0,
    descriptors: Optional[Tuple[str, ...]] = None,
) -> pd.DataFrame:
    """Compute Boltzmann-weighted ensemble averages at given temperature.

    For each geometric descriptor, computes::

        <X>_boltz = sum(X_i * w_i) / sum(w_i)

    where ``w_i = exp(-dE_i / kT)`` and ``dE_i = E_i - E_min``.

    Args:
        df: Per-conformer DataFrame with energy data.
        group_col: Column to group by.
        energy_col: Name of the energy column.
        temperature_K: Temperature in Kelvin (default: 298K).
        descriptors: Columns to weight. Defaults to :data:`GEOMETRIC_DESCRIPTORS`.

    Returns:
        DataFrame with ``{descriptor}_boltz`` columns and
        ``boltz_excimer_fraction`` if classification data is available.

    Notes:
        MMFF94s energies have ~0.5-1.0 kcal/mol accuracy, which limits
        the reliability of Boltzmann weighting. Results should be
        interpreted as qualitative trends rather than quantitative
        predictions.
    """
    if df.empty:
        return pd.DataFrame()

    if energy_col not in df.columns:
        warnings.warn(
            f"Energy column '{energy_col}' not found. "
            "Cannot compute Boltzmann-weighted features.",
            UserWarning,
            stacklevel=2,
        )
        return pd.DataFrame(index=df[group_col].unique())

    if descriptors is None:
        descriptors = GEOMETRIC_DESCRIPTORS

    available = [d for d in descriptors if d in df.columns]
    kT = 0.001987204 * temperature_K  # R in kcal/(mol*K)

    grouped = df.groupby(group_col)
    results = {}

    for name, group in grouped:
        # Filter to rows with valid energy
        energy_valid = group[energy_col].notna()
        valid_group = group[energy_valid]
        if len(valid_group) < 2:
            continue

        energies = valid_group[energy_col]
        dE = energies - energies.min()
        weights = np.exp(-dE / kT)
        total_w = weights.sum()
        if total_w == 0 or not np.isfinite(total_w):
            continue

        weights_norm = weights / total_w
        row = {}

        for desc in available:
            vals = valid_group[desc].values
            w_vals = weights_norm.values
            valid_mask = np.isfinite(vals) & np.isfinite(w_vals)
            if valid_mask.any():
                row[f"{desc}_boltz"] = np.sum(
                    vals[valid_mask] * w_vals[valid_mask]
                ) / w_vals[valid_mask].sum()
            else:
                row[f"{desc}_boltz"] = np.nan

        # Boltzmann-weighted excimer fraction
        if "classification" in valid_group.columns:
            exc_mask = valid_group["classification"].isin(
                ["strong_excimer", "weak_excimer"]
            )
            exc_binary = exc_mask.astype(float).values
            w_vals = weights_norm.values
            valid_mask = np.isfinite(w_vals)
            if valid_mask.any():
                row["boltz_excimer_fraction"] = np.sum(
                    exc_binary[valid_mask] * w_vals[valid_mask]
                ) / w_vals[valid_mask].sum()

        # Boltzmann-weighted excimer score
        if "excimer_score" in valid_group.columns:
            scores = valid_group["excimer_score"].values
            w_vals = weights_norm.values
            valid_mask = np.isfinite(scores) & np.isfinite(w_vals)
            if valid_mask.any():
                row["boltz_excimer_score"] = np.sum(
                    scores[valid_mask] * w_vals[valid_mask]
                ) / w_vals[valid_mask].sum()

        results[name] = row

    if not results:
        return pd.DataFrame(index=df[group_col].unique())

    result = pd.DataFrame.from_dict(results, orient="index")
    result.index.name = group_col
    return result


def compute_ensemble_features(
    df: pd.DataFrame,
    group_col: str = "molecule",
    energy_col: str = "energy_kcal_mol",
    temperature_K: float = 298.0,
    thresholds: Optional[Dict] = None,
    descriptors: Optional[Tuple[str, ...]] = None,
) -> pd.DataFrame:
    """Compute all ensemble features: distributional + threshold + Boltzmann.

    This is the main entry point. Combines output of the three feature
    computation functions into a single wide-format DataFrame.

    Args:
        df: Per-conformer analysis DataFrame. Should include classification
            and optionally excimer_score columns.
        group_col: Column to group by (``"molecule"`` or ``"substituent"``).
        energy_col: Energy column for Boltzmann weighting.
        temperature_K: Temperature in Kelvin.
        thresholds: Optional threshold overrides.
        descriptors: Optional descriptor column overrides.

    Returns:
        DataFrame indexed by group_col with ~55 feature columns:
        distributional + threshold + Boltzmann.

    Example::

        >>> from pyrene_analyzer.core import AromaticDimerAnalyzer
        >>> from pyrene_analyzer.ensemble import compute_ensemble_features
        >>> analyzer = AromaticDimerAnalyzer()
        >>> results = analyzer.analyze_file("Et.sdf")
        >>> results = analyzer.add_classification(results)
        >>> features = compute_ensemble_features(results)
    """
    if df.empty:
        return pd.DataFrame()

    dist_feats = compute_distributional_features(
        df, group_col=group_col, descriptors=descriptors
    )
    thresh_feats = compute_threshold_features(
        df, group_col=group_col, thresholds=thresholds
    )
    boltz_feats = compute_boltzmann_weighted_features(
        df,
        group_col=group_col,
        energy_col=energy_col,
        temperature_K=temperature_K,
        descriptors=descriptors,
    )

    # Merge all feature sets
    result = dist_feats
    if not thresh_feats.empty:
        result = result.join(thresh_feats, how="outer")
    if not boltz_feats.empty:
        result = result.join(boltz_feats, how="outer")

    return result
