"""
Substituent Physicochemical Parameters
======================================

Provides Taft steric (Es), Charton steric (v), Hammett electronic
(sigma_meta, sigma_para), and molar refractivity (MR) parameters
for common substituents used in SAR analysis.

Scientific Context:
    Literature (Wheeler-Houk 2011, Carter-Fenk & Herbert 2020) shows that
    steric effects (Taft Es) dominate over electronic effects (Hammett sigma)
    in intramolecular aromatic stacking systems like binaphthalene dimers.

    Key finding from our screening:
    - CF3 vs CN have similar sigma_para (0.54 vs 0.66)
    - But very different Es (-2.40 vs -0.51)
    - And very different excimer fractions (0% vs 12.5%)
    - => Steric control, not electronic

Sources:
    - Taft Es: Taft (1956), MacPhee & Dubois (1968)
    - Charton v: Charton (1975, 1983)
    - Hammett sigma: Hansch, Leo & Taft (1991) Chem. Rev.
    - MR: Hansch & Leo (1979) Substituent Constants for
          Correlation Analysis in Chemistry and Biology

Example:
    >>> from pyrene_analyzer.substituent_params import get_params, PARAMS_DF
    >>> params = get_params("tBu")
    >>> print(f"tBu: Es={params['Es']}, sigma_para={params['sigma_para']}")
    tBu: Es=-2.78, sigma_para=-0.20
"""

from typing import Dict, List, Optional, Union
import pandas as pd


# Parameter columns: (Es, v, sigma_meta, sigma_para, MR)
# Es: Taft steric parameter (more negative = bulkier)
# v: Charton steric parameter (more positive = bulkier)
# sigma_meta: Hammett sigma constant for meta position
# sigma_para: Hammett sigma constant for para position
# MR: Molar refractivity (cm^3/mol, related to polarizability)

SUBSTITUENT_PARAMS: Dict[str, Dict[str, float]] = {
    # === Alkyl Series (steric variation) ===
    "H": {
        "Es": 0.00, "v": 0.00,
        "sigma_meta": 0.00, "sigma_para": 0.00,
        "MR": 1.03,
        "type": "reference",
    },
    "Me": {
        "Es": -1.24, "v": 0.52,
        "sigma_meta": -0.07, "sigma_para": -0.17,
        "MR": 5.65,
        "type": "alkyl",
    },
    "Et": {
        "Es": -1.31, "v": 0.56,
        "sigma_meta": -0.07, "sigma_para": -0.15,
        "MR": 10.30,
        "type": "alkyl",
    },
    "nPr": {
        "Es": -1.60, "v": 0.68,
        "sigma_meta": -0.07, "sigma_para": -0.13,
        "MR": 14.96,
        "type": "alkyl",
    },
    "iPr": {
        "Es": -1.71, "v": 0.76,
        "sigma_meta": -0.07, "sigma_para": -0.15,
        "MR": 14.96,
        "type": "alkyl",
    },
    "nBu": {
        "Es": -1.63, "v": 0.68,
        "sigma_meta": -0.08, "sigma_para": -0.16,
        "MR": 19.61,
        "type": "alkyl",
    },
    "tBu": {
        "Es": -2.78, "v": 1.24,
        "sigma_meta": -0.10, "sigma_para": -0.20,
        "MR": 19.62,
        "type": "alkyl",
    },
    "cHex": {
        "Es": -2.03, "v": 0.87,
        "sigma_meta": -0.05, "sigma_para": -0.15,
        "MR": 25.36,
        "type": "alkyl",
    },
    "cPr": {
        "Es": -0.47, "v": 0.54,
        "sigma_meta": -0.06, "sigma_para": -0.11,
        "MR": 13.38,
        "type": "alkyl",
    },

    # === Electron-Donating Groups (EDGs) ===
    "MeO": {
        "Es": -0.55, "v": 0.36,
        "sigma_meta": 0.12, "sigma_para": -0.27,
        "MR": 7.87,
        "type": "EDG",
    },
    "OMe": {  # Alias for MeO
        "Es": -0.55, "v": 0.36,
        "sigma_meta": 0.12, "sigma_para": -0.27,
        "MR": 7.87,
        "type": "EDG",
    },
    "OEt": {
        "Es": -0.55, "v": 0.36,
        "sigma_meta": 0.10, "sigma_para": -0.24,
        "MR": 12.47,
        "type": "EDG",
    },
    "OH": {
        "Es": -0.55, "v": 0.35,
        "sigma_meta": 0.12, "sigma_para": -0.37,
        "MR": 2.85,
        "type": "EDG",
    },
    "NH2": {
        "Es": -0.61, "v": 0.35,
        "sigma_meta": -0.16, "sigma_para": -0.66,
        "MR": 5.42,
        "type": "EDG",
    },
    "NMe2": {
        "Es": -0.93, "v": 0.61,
        "sigma_meta": -0.15, "sigma_para": -0.83,
        "MR": 15.55,
        "type": "EDG",
    },

    # === Halogens ===
    "F": {
        "Es": -0.46, "v": 0.27,
        "sigma_meta": 0.34, "sigma_para": 0.06,
        "MR": 0.92,
        "type": "halogen",
    },
    "Cl": {
        "Es": -0.97, "v": 0.55,
        "sigma_meta": 0.37, "sigma_para": 0.23,
        "MR": 6.03,
        "type": "halogen",
    },
    "Br": {
        "Es": -1.16, "v": 0.65,
        "sigma_meta": 0.39, "sigma_para": 0.23,
        "MR": 8.88,
        "type": "halogen",
    },
    "I": {
        "Es": -1.40, "v": 0.78,
        "sigma_meta": 0.35, "sigma_para": 0.18,
        "MR": 13.94,
        "type": "halogen",
    },

    # === Electron-Withdrawing Groups (EWGs) ===
    "CF3": {
        "Es": -2.40, "v": 0.91,
        "sigma_meta": 0.43, "sigma_para": 0.54,
        "MR": 5.02,
        "type": "EWG",
    },
    "CN": {
        "Es": -0.51, "v": 0.40,
        "sigma_meta": 0.56, "sigma_para": 0.66,
        "MR": 6.33,
        "type": "EWG",
    },
    "NO2": {
        "Es": -1.01, "v": 0.50,
        "sigma_meta": 0.71, "sigma_para": 0.78,
        "MR": 7.36,
        "type": "EWG",
    },
    "COMe": {
        "Es": -1.10, "v": 0.50,
        "sigma_meta": 0.38, "sigma_para": 0.50,
        "MR": 11.18,
        "type": "EWG",
    },
    "CO2Me": {
        "Es": -1.05, "v": 0.45,
        "sigma_meta": 0.37, "sigma_para": 0.45,
        "MR": 12.87,
        "type": "EWG",
    },
    "CHO": {
        "Es": -0.52, "v": 0.40,
        "sigma_meta": 0.36, "sigma_para": 0.42,
        "MR": 6.88,
        "type": "EWG",
    },

    # === Aromatic ===
    "Ph": {
        "Es": -2.55, "v": 1.66,
        "sigma_meta": 0.06, "sigma_para": -0.01,
        "MR": 25.36,
        "type": "aromatic",
    },
    "Bn": {
        "Es": -1.98, "v": 0.70,
        "sigma_meta": 0.08, "sigma_para": 0.00,
        "MR": 30.01,
        "type": "aromatic",
    },
}


def get_params(substituent: str) -> Optional[Dict[str, float]]:
    """
    Get all physicochemical parameters for a substituent.

    Args:
        substituent: Substituent name (e.g., "Me", "tBu", "CF3").

    Returns:
        Dict with keys: Es, v, sigma_meta, sigma_para, MR, type.
        Returns None if substituent not found.

    Example:
        >>> params = get_params("tBu")
        >>> print(params["Es"], params["v"])
        -2.78 1.24
    """
    return SUBSTITUENT_PARAMS.get(substituent)


def get_param(substituent: str, param_name: str) -> Optional[float]:
    """
    Get a specific parameter for a substituent.

    Args:
        substituent: Substituent name.
        param_name: Parameter name (Es, v, sigma_meta, sigma_para, MR).

    Returns:
        Parameter value, or None if not found.

    Example:
        >>> get_param("CF3", "Es")
        -2.40
        >>> get_param("CF3", "sigma_para")
        0.54
    """
    params = SUBSTITUENT_PARAMS.get(substituent)
    if params is None:
        return None
    return params.get(param_name)


def list_substituents(sub_type: Optional[str] = None) -> List[str]:
    """
    List available substituents, optionally filtered by type.

    Args:
        sub_type: Filter by type (alkyl, EDG, EWG, halogen, aromatic).
            If None, returns all substituents.

    Returns:
        List of substituent names.

    Example:
        >>> list_substituents("alkyl")
        ['Me', 'Et', 'nPr', 'iPr', 'nBu', 'tBu', 'cHex', 'cPr']
    """
    if sub_type is None:
        return list(SUBSTITUENT_PARAMS.keys())

    return [
        name for name, params in SUBSTITUENT_PARAMS.items()
        if params.get("type") == sub_type
    ]


def get_params_dataframe() -> pd.DataFrame:
    """
    Get all parameters as a pandas DataFrame.

    Returns:
        DataFrame with substituent names as index and parameters as columns.

    Example:
        >>> df = get_params_dataframe()
        >>> print(df.loc[["Me", "tBu", "CF3"], ["Es", "sigma_para"]])
                  Es  sigma_para
        Me    -1.24       -0.17
        tBu   -2.78       -0.20
        CF3   -2.40        0.54
    """
    return pd.DataFrame.from_dict(SUBSTITUENT_PARAMS, orient="index")


# Pre-computed DataFrame for convenience
PARAMS_DF = get_params_dataframe()


def add_params_to_dataframe(
    df: pd.DataFrame,
    substituent_col: str = "substituent",
    params: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Add substituent parameters as new columns to a DataFrame.

    Args:
        df: DataFrame with a substituent column.
        substituent_col: Name of column containing substituent names.
        params: List of parameters to add. If None, adds all:
            ["Es", "v", "sigma_meta", "sigma_para", "MR"].

    Returns:
        DataFrame with parameter columns added.

    Example:
        >>> results = pd.DataFrame({
        ...     "substituent": ["Me", "Et", "tBu"],
        ...     "excimer_fraction": [0.40, 0.35, 0.05],
        ... })
        >>> results = add_params_to_dataframe(results)
        >>> print(results[["substituent", "Es", "excimer_fraction"]])
    """
    if params is None:
        params = ["Es", "v", "sigma_meta", "sigma_para", "MR"]

    df = df.copy()

    for param in params:
        df[param] = df[substituent_col].apply(
            lambda x: get_param(x, param)
        )

    return df


def calculate_steric_correlation(
    df: pd.DataFrame,
    substituent_col: str = "substituent",
    target_col: str = "excimer_fraction",
) -> Dict[str, float]:
    """
    Calculate correlation between steric parameters and a target variable.

    Uses Spearman rank correlation for robustness to outliers.

    Args:
        df: DataFrame with substituent and target columns.
        substituent_col: Column with substituent names.
        target_col: Column with target variable (e.g., excimer_fraction).

    Returns:
        Dict with correlation coefficients and p-values for Es and v.

    Example:
        >>> results = calculate_steric_correlation(screening_df)
        >>> print(f"Es correlation: r={results['Es_r']:.3f}, p={results['Es_p']:.4f}")
    """
    from scipy.stats import spearmanr, pearsonr

    # Add parameters if not present
    if "Es" not in df.columns:
        df = add_params_to_dataframe(df, substituent_col, ["Es", "v"])

    # Drop rows with missing values
    df_clean = df.dropna(subset=["Es", "v", target_col])

    results = {}

    # Es correlation
    r_es, p_es = spearmanr(df_clean["Es"], df_clean[target_col])
    results["Es_r"] = r_es
    results["Es_p"] = p_es

    # Charton v correlation
    r_v, p_v = spearmanr(df_clean["v"], df_clean[target_col])
    results["v_r"] = r_v
    results["v_p"] = p_v

    # Pearson for comparison
    r_es_pearson, _ = pearsonr(df_clean["Es"], df_clean[target_col])
    results["Es_r_pearson"] = r_es_pearson

    return results


# Convenience: mapping of screening group names to canonical names
SCREENING_GROUP_ALIASES = {
    "H": "H",
    "Me": "Me",
    "Et": "Et",
    "nPr": "nPr",
    "iPr": "iPr",
    "nBu": "nBu",
    "tBu": "tBu",
    "cHex": "cHex",
    "MeO": "MeO",
    "OMe": "MeO",
    "OEt": "OEt",
    "F": "F",
    "Cl": "Cl",
    "CF3": "CF3",
    "CN": "CN",
    "NMe2": "NMe2",
    "Ph": "Ph",
}


def normalize_substituent_name(name: str) -> str:
    """
    Normalize substituent name to canonical form.

    Args:
        name: Substituent name (possibly with variations).

    Returns:
        Canonical name used in SUBSTITUENT_PARAMS.

    Example:
        >>> normalize_substituent_name("OMe")
        'MeO'
    """
    return SCREENING_GROUP_ALIASES.get(name, name)
