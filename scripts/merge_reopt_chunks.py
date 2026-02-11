#!/usr/bin/env python
"""Merge chunked re-optimization outputs into final tables."""

from __future__ import annotations

import argparse
import glob
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def _safe_col_mean(df: pd.DataFrame, col: str) -> float:
    if col not in df.columns:
        return float("nan")
    return float(df[col].mean())


def _safe_col_max(df: pd.DataFrame, col: str) -> float:
    if col not in df.columns:
        return float("nan")
    return float(df[col].max())


def _build_summary(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()

    records = []
    for name, group in df.groupby("name"):
        n_conformers = int(len(group))
        if "classification" in group.columns:
            n_excimer = int(
                group["classification"].isin(["strong_excimer", "weak_excimer"]).sum()
            )
        else:
            n_excimer = 0
        excimer_fraction = n_excimer / n_conformers if n_conformers else 0.0

        records.append(
            {
                "name": name,
                "r_group": group["r_group"].iloc[0] if "r_group" in group.columns else "",
                "screen_group": (
                    group["screen_group"].iloc[0] if "screen_group" in group.columns else ""
                ),
                "n_conformers": n_conformers,
                "n_excimer": n_excimer,
                "excimer_fraction": excimer_fraction,
                "mean_distance_moe": _safe_col_mean(group, "moe_interplane_distance_A"),
                "mean_distance_mace": _safe_col_mean(group, "interplane_distance_A"),
                "mean_distance_delta": _safe_col_mean(group, "distance_delta_A"),
                "mean_angle": _safe_col_mean(group, "plane_angle_deg"),
                "mean_overlap": _safe_col_mean(group, "pi_overlap_pct"),
                "best_overlap": _safe_col_max(group, "pi_overlap_pct"),
                "mean_energy": _safe_col_mean(group, "mace_energy_kcal"),
            }
        )

    summary_df = pd.DataFrame(records)
    if "excimer_fraction" in summary_df.columns:
        summary_df = summary_df.sort_values("excimer_fraction", ascending=False)
    return summary_df


def _write_boltzmann_tables(df: pd.DataFrame, output_prefix: Path) -> None:
    if df.empty or "mace_energy_kcal" not in df.columns:
        return
    try:
        from pyrene_analyzer.ensemble import compute_boltzmann_weighted_features
    except ImportError:
        print("Skipping Boltzmann tables: pyrene_analyzer.ensemble not importable")
        return

    for temperature in (200, 298, 350, 400):
        try:
            boltz = compute_boltzmann_weighted_features(
                df,
                group_col="name",
                energy_col="mace_energy_kcal",
                temperature_K=temperature,
            )
            out = output_prefix.parent / f"{output_prefix.name}_boltzmann_{temperature}K.csv"
            boltz.to_csv(out, index=False, encoding="utf-8-sig")
            print(f"Wrote {out}")
        except Exception as exc:
            print(f"Boltzmann {temperature}K failed: {exc}")


def _read_parts(paths: Iterable[Path]) -> pd.DataFrame:
    dfs = []
    for path in paths:
        try:
            df = pd.read_csv(path)
            df["source_file"] = str(path)
            dfs.append(df)
            print(f"Loaded {path} ({len(df)} rows)")
        except Exception as exc:
            print(f"Skipping {path}: {exc}")

    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-glob",
        required=True,
        help="Glob for chunk outputs, e.g. runs/moe_mace_reopt_gpu*_all_conformers.csv",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for merged outputs (default: current directory)",
    )
    parser.add_argument(
        "--output-prefix",
        default="moe_mace_reopt",
        help="Prefix for merged files (default: moe_mace_reopt)",
    )
    parser.add_argument(
        "--no-boltzmann",
        action="store_true",
        help="Skip Boltzmann table generation",
    )
    args = parser.parse_args()

    part_paths = [Path(p) for p in sorted(glob.glob(args.input_glob))]
    if not part_paths:
        print(f"No files matched glob: {args.input_glob}")
        return 1

    merged = _read_parts(part_paths)
    if merged.empty:
        print("No valid rows found in chunk inputs")
        return 1

    dedupe_cols = [c for c in ("name", "conformer_id") if c in merged.columns]
    if dedupe_cols:
        merged = merged.drop_duplicates(subset=dedupe_cols, keep="last")

    if (
        "name" in merged.columns
        and "conformer_id" in merged.columns
        and np.issubdtype(merged["conformer_id"].dtype, np.number)
    ):
        merged = merged.sort_values(["name", "conformer_id"])

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / args.output_prefix

    all_csv = out_dir / f"{args.output_prefix}_all_conformers.csv"
    merged.to_csv(all_csv, index=False, encoding="utf-8-sig")
    print(f"Wrote {all_csv} ({len(merged)} rows)")

    summary = _build_summary(merged)
    summary_csv = out_dir / f"{args.output_prefix}_summary.csv"
    summary.to_csv(summary_csv, index=False, encoding="utf-8-sig")
    print(f"Wrote {summary_csv} ({len(summary)} rows)")

    if not args.no_boltzmann:
        _write_boltzmann_tables(merged, out_prefix)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
