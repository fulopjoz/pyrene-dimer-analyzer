#!/usr/bin/env python
"""Select a TD-DFT candidate subset from MACE re-optimization outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


EXCIMER_CLASSES = {"strong_excimer", "weak_excimer"}


def _require_columns(df: pd.DataFrame, required: list[str]) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")


def _clip01(series: pd.Series) -> pd.Series:
    return series.clip(lower=0.0, upper=1.0)


def compute_priority_score(df: pd.DataFrame) -> pd.Series:
    """Compute a ranking score based on geometric excimer propensity."""
    score = pd.Series(0.0, index=df.index, dtype=float)

    if "classification" in df.columns:
        cls = df["classification"].fillna("")
        score += (cls == "strong_excimer").astype(float) * 2.0
        score += (cls == "weak_excimer").astype(float) * 1.0

    if "pi_overlap_pct" in df.columns:
        score += _clip01(df["pi_overlap_pct"].fillna(0.0) / 100.0) * 2.0

    if "interplane_distance_A" in df.columns:
        score -= (df["interplane_distance_A"].fillna(99.0) - 3.45).abs() * 0.8

    if "plane_angle_deg" in df.columns:
        score -= _clip01(df["plane_angle_deg"].fillna(90.0) / 90.0) * 0.5

    # Boost thermally significant conformers when present.
    for col in ("boltz_weight", "boltz_weight_298K", "weight_298K"):
        if col in df.columns:
            score += df[col].fillna(0.0) * 1.5
            break

    # Prefer lower energies when available.
    for col in ("mace_energy_kcal", "energy_kcal_mol"):
        if col in df.columns:
            centered = df[col].fillna(df[col].median())
            score -= (centered - centered.min()) * 0.02
            break

    return score


def select_candidates(
    df: pd.DataFrame,
    top_per_molecule: int,
    max_total: int,
    max_theta: float,
    max_distance: float,
    min_overlap: float,
    require_excimer_class: bool,
) -> pd.DataFrame:
    work = df.copy()

    mask = pd.Series(True, index=work.index)
    if "plane_angle_deg" in work.columns:
        mask &= work["plane_angle_deg"] <= max_theta
    if "interplane_distance_A" in work.columns:
        mask &= work["interplane_distance_A"] <= max_distance
    if "pi_overlap_pct" in work.columns:
        mask &= work["pi_overlap_pct"] >= min_overlap
    if require_excimer_class and "classification" in work.columns:
        mask &= work["classification"].isin(EXCIMER_CLASSES)

    filtered = work.loc[mask].copy()
    if filtered.empty:
        return filtered

    filtered["orca_priority_score"] = compute_priority_score(filtered)
    filtered = filtered.sort_values("orca_priority_score", ascending=False)

    if "name" in filtered.columns:
        filtered = (
            filtered.groupby("name", group_keys=False).head(top_per_molecule).copy()
        )
        filtered = filtered.sort_values("orca_priority_score", ascending=False)

    filtered = filtered.head(max_total).copy()
    filtered["orca_candidate_rank"] = range(1, len(filtered) + 1)
    return filtered


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--input-csv",
        default="runs/moe_mace_reopt_all_conformers.csv",
        help="Merged conformer CSV from reoptimize workflow",
    )
    p.add_argument(
        "--output-csv",
        default="runs/orca_candidates.csv",
        help="Output candidate CSV for ORCA preparation",
    )
    p.add_argument(
        "--top-per-molecule",
        type=int,
        default=2,
        help="Max candidates retained per molecule name (default: 2)",
    )
    p.add_argument(
        "--max-total",
        type=int,
        default=96,
        help="Global cap on selected conformers (default: 96)",
    )
    p.add_argument(
        "--max-theta",
        type=float,
        default=50.0,
        help="Maximum plane angle in degrees (default: 50.0)",
    )
    p.add_argument(
        "--max-distance",
        type=float,
        default=4.6,
        help="Maximum interplane distance in Angstrom (default: 4.6)",
    )
    p.add_argument(
        "--min-overlap",
        type=float,
        default=30.0,
        help="Minimum pi-overlap percentage (default: 30.0)",
    )
    p.add_argument(
        "--require-excimer-class",
        action="store_true",
        help="Keep only weak/strong excimer classifications",
    )
    return p


def main() -> int:
    args = build_parser().parse_args()

    in_path = Path(args.input_csv)
    if not in_path.exists():
        raise SystemExit(f"Input CSV not found: {in_path}")

    df = pd.read_csv(in_path)
    _require_columns(df, ["conformer_id"])

    selected = select_candidates(
        df=df,
        top_per_molecule=max(args.top_per_molecule, 1),
        max_total=max(args.max_total, 1),
        max_theta=args.max_theta,
        max_distance=args.max_distance,
        min_overlap=args.min_overlap,
        require_excimer_class=bool(args.require_excimer_class),
    )
    if selected.empty:
        raise SystemExit("No candidates passed the selected filters.")

    out_path = Path(args.output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(out_path, index=False, encoding="utf-8-sig")

    n_mols = selected["name"].nunique() if "name" in selected.columns else "N/A"
    print(f"Wrote {out_path} with {len(selected)} candidates (molecules={n_mols})")
    print(
        "Score range:",
        f"{selected['orca_priority_score'].min():.3f}",
        "to",
        f"{selected['orca_priority_score'].max():.3f}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
