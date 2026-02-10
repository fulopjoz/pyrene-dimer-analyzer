"""
Three-way comparison of conformer generation methods.

Compares geometric analysis results across three pipelines:
1. RDKit ETKDGv3 + MMFF94s (fast, no dispersion)
2. MOE LowModeMD + Amber:EHT (commercial, partial dispersion)
3. MOE + MACE-OFF23 re-optimization (DFT-quality dispersion)

Generates comparison plots and summary statistics to evaluate method
agreement and identify force field artifacts.

Scientific basis:
    Classical force fields (MMFF94s, Amber) lack or underestimate London
    dispersion, systematically overestimating pi-stacking distances by
    ~1-2 A. MACE-OFF23 (trained on wB97M-D3BJ/def2-TZVPPD) corrects
    this, giving distances consistent with CCSD(T)/CBS benchmarks.

    Ge et al. (2020) J. Mater. Chem. C showed that pi-overlap is more
    predictive than distance for excimer formation. Method comparison
    reveals which pipeline provides the most reliable geometric analysis.

Usage:
    python compare_methods.py
    python compare_methods.py --rdkit binaph_screening_all_conformers.csv
    python compare_methods.py --output-dir plots/comparison

References:
    - Batatia et al. (2024) arXiv:2401.00096 (MACE-OFF23)
    - Ge et al. (2020) J. Mater. Chem. C 8, 10223
    - PMC11476719 (2024) DLPNO-CCSD(T)/CBS benchmark
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def load_all_results(
    rdkit_csv: str,
    moe_csv: str,
    mace_csv: str = None,
) -> dict:
    """Load results from all available methods.

    Returns dict of method_name -> DataFrame.
    """
    data = {}

    if Path(rdkit_csv).exists():
        df = pd.read_csv(rdkit_csv)
        df["method"] = "RDKit/MMFF94s"
        data["RDKit/MMFF94s"] = df
        print(f"  RDKit/MMFF94s: {len(df)} conformers")

    if Path(moe_csv).exists():
        df = pd.read_csv(moe_csv)
        df["method"] = "MOE/Amber:EHT"
        data["MOE/Amber:EHT"] = df
        print(f"  MOE/Amber:EHT: {len(df)} conformers")

    if mace_csv and Path(mace_csv).exists():
        df = pd.read_csv(mace_csv)
        df["method"] = "MOE+MACE-OFF23"
        data["MOE+MACE-OFF23"] = df
        print(f"  MOE+MACE-OFF23: {len(df)} conformers")

    return data


def compute_comparison_summary(data: dict) -> pd.DataFrame:
    """Compute per-method summary statistics.

    Returns DataFrame with method-level comparison.
    """
    records = []
    for method, df in data.items():
        has_class = "classification" in df.columns
        n_exc = 0
        if has_class:
            n_exc = sum(df["classification"].isin(["strong_excimer", "weak_excimer"]))

        dist_col = "interplane_distance_A"
        angle_col = "plane_angle_deg"
        overlap_col = "pi_overlap_pct"

        records.append({
            "method": method,
            "n_conformers": len(df),
            "n_molecules": df["name"].nunique() if "name" in df.columns else df.get("molecule", df.get("substituent", pd.Series())).nunique(),
            "n_excimer": n_exc,
            "excimer_fraction": n_exc / len(df) if len(df) > 0 else 0,
            "mean_distance": df[dist_col].mean() if dist_col in df.columns else np.nan,
            "std_distance": df[dist_col].std() if dist_col in df.columns else np.nan,
            "median_distance": df[dist_col].median() if dist_col in df.columns else np.nan,
            "mean_angle": df[angle_col].mean() if angle_col in df.columns else np.nan,
            "mean_overlap": df[overlap_col].mean() if overlap_col in df.columns else np.nan,
            "frac_distance_lt_4": (
                (df[dist_col] < 4.0).mean() if dist_col in df.columns else np.nan
            ),
            "frac_overlap_gt_40": (
                (df[overlap_col] > 40).mean() if overlap_col in df.columns else np.nan
            ),
        })

    return pd.DataFrame(records)


def plot_distance_kde(data: dict, output: Path) -> None:
    """KDE plot of interplane distances by method."""
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {"RDKit/MMFF94s": "#e74c3c", "MOE/Amber:EHT": "#3498db", "MOE+MACE-OFF23": "#2ecc71"}

    for method, df in data.items():
        col = "interplane_distance_A"
        if col not in df.columns:
            continue
        vals = df[col].dropna()
        if len(vals) < 2:
            continue
        color = colors.get(method, "#95a5a6")
        ax.hist(vals, bins=50, alpha=0.3, color=color, density=True, label=f"{method} (n={len(vals)})")
        try:
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(vals)
            x = np.linspace(vals.min(), vals.max(), 200)
            ax.plot(x, kde(x), color=color, linewidth=2)
        except Exception:
            pass

    # Reference lines
    ax.axvline(3.43, color="black", linestyle="--", linewidth=1, alpha=0.7,
               label="CCSD(T)/CBS: 3.43 A")
    ax.axvline(4.5, color="gray", linestyle=":", linewidth=1, alpha=0.7,
               label="Weak excimer: 4.5 A")

    ax.set_xlabel("Interplane Distance (A)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title("Interplane Distance Distribution by Method", fontsize=14)
    ax.legend(fontsize=9)
    ax.set_xlim(0, 15)
    fig.tight_layout()
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output}")


def plot_excimer_heatmap(data: dict, output: Path) -> None:
    """Heatmap of excimer fraction by R-group x screening group per method."""
    methods = list(data.keys())
    n_methods = len(methods)
    fig, axes = plt.subplots(1, n_methods, figsize=(7 * n_methods, 6), squeeze=False)

    for i, method in enumerate(methods):
        ax = axes[0, i]
        df = data[method]

        # Need r_group and screen_group columns
        if "r_group" not in df.columns or "screen_group" not in df.columns:
            ax.set_title(f"{method}\n(no R-group data)")
            continue

        has_class = "classification" in df.columns
        if not has_class:
            ax.set_title(f"{method}\n(no classification)")
            continue

        # Compute excimer fraction per (r_group, screen_group)
        df["is_excimer"] = df["classification"].isin(["strong_excimer", "weak_excimer"])
        pivot = df.groupby(["r_group", "screen_group"])["is_excimer"].mean().unstack(fill_value=0)

        sns.heatmap(
            pivot, annot=True, fmt=".0%", cmap="YlOrRd", vmin=0, vmax=0.5,
            ax=ax, cbar_kws={"label": "Excimer Fraction"}
        )
        ax.set_title(f"{method}", fontsize=12)
        ax.set_xlabel("Screening Group")
        ax.set_ylabel("R-Group")

    fig.suptitle("Excimer Fraction by R-Group x Screening Group", fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output}")


def plot_distance_correction(moe_df: pd.DataFrame, mace_df: pd.DataFrame, output: Path) -> None:
    """Scatter plot of MOE vs MACE distances (shows systematic correction)."""
    if "moe_interplane_distance_A" not in mace_df.columns:
        print("  Skipping distance correction plot (no MOE distance in MACE data)")
        return

    fig, ax = plt.subplots(figsize=(8, 8))

    x = mace_df["moe_interplane_distance_A"]
    y = mace_df["interplane_distance_A"]

    # Color by classification if available
    if "classification" in mace_df.columns:
        colors = {
            "strong_excimer": "#2ecc71",
            "weak_excimer": "#f39c12",
            "monomer": "#95a5a6",
        }
        for cls, color in colors.items():
            mask = mace_df["classification"] == cls
            ax.scatter(x[mask], y[mask], c=color, alpha=0.4, s=10, label=cls)
    else:
        ax.scatter(x, y, alpha=0.3, s=10, c="#3498db")

    # Diagonal (no change)
    lim = max(x.max(), y.max()) * 1.05
    ax.plot([0, lim], [0, lim], "k--", alpha=0.3, label="No change")

    # Reference lines
    ax.axhline(3.43, color="green", linestyle=":", alpha=0.5, label="CCSD(T)/CBS: 3.43 A")
    ax.axhline(4.5, color="orange", linestyle=":", alpha=0.5, label="Weak threshold: 4.5 A")

    ax.set_xlabel("MOE/Amber:EHT Distance (A)", fontsize=12)
    ax.set_ylabel("MACE-OFF23 Distance (A)", fontsize=12)
    ax.set_title("Distance Correction: MOE vs MACE-OFF23", fontsize=14)
    ax.legend(fontsize=9)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output}")


def plot_r_group_comparison(data: dict, output: Path) -> None:
    """Bar plot of excimer fraction by R-group, one group of bars per method."""
    records = []
    for method, df in data.items():
        if "r_group" not in df.columns or "classification" not in df.columns:
            continue
        df["is_excimer"] = df["classification"].isin(["strong_excimer", "weak_excimer"])
        for rg, group in df.groupby("r_group"):
            records.append({
                "method": method,
                "r_group": rg,
                "excimer_fraction": group["is_excimer"].mean(),
                "n_conformers": len(group),
            })

    if not records:
        print("  Skipping R-group comparison (no classification data)")
        return

    plot_df = pd.DataFrame(records)

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(
        data=plot_df, x="r_group", y="excimer_fraction", hue="method",
        ax=ax, palette="Set2",
    )
    ax.set_xlabel("Luminescent R-Group", fontsize=12)
    ax.set_ylabel("Excimer Fraction", fontsize=12)
    ax.set_title("Excimer Fraction by R-Group Across Methods", fontsize=14)
    ax.legend(title="Method", fontsize=9)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.0%}"))
    fig.tight_layout()
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output}")


def main():
    parser = argparse.ArgumentParser(
        description="Three-way comparison of conformer methods"
    )
    parser.add_argument(
        "--rdkit",
        default="binaph_screening_all_conformers.csv",
        help="RDKit/MMFF94s results CSV",
    )
    parser.add_argument(
        "--moe",
        default="moe_screening_all_conformers.csv",
        help="MOE/Amber:EHT results CSV",
    )
    parser.add_argument(
        "--mace",
        default="moe_mace_reopt_all_conformers.csv",
        help="MOE+MACE-OFF23 results CSV",
    )
    parser.add_argument(
        "--output-dir",
        default="plots/comparison",
        help="Output directory for plots",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading results...")
    data = load_all_results(args.rdkit, args.moe, args.mace)

    if not data:
        print("ERROR: No result files found!")
        sys.exit(1)

    # Compute summary
    print("\nComputing comparison summary...")
    summary = compute_comparison_summary(data)
    summary_csv = output_dir / "comparison_summary.csv"
    summary.to_csv(summary_csv, index=False, encoding="utf-8-sig")
    print(f"  Saved: {summary_csv}")

    # Print summary table
    print(f"\n{'Method':<20} {'n_conf':>7} {'n_mol':>5} {'n_exc':>5} "
          f"{'frac':>6} {'d_mean':>6} {'d_med':>5} {'olap':>5}")
    print("-" * 72)
    for _, row in summary.iterrows():
        print(
            f"{row['method']:<20} {row['n_conformers']:>7d} "
            f"{row['n_molecules']:>5.0f} {row['n_excimer']:>5d} "
            f"{row['excimer_fraction']:>5.1%} {row['mean_distance']:>6.2f} "
            f"{row['median_distance']:>5.2f} {row['mean_overlap']:>5.1f}"
        )

    # Generate plots
    print("\nGenerating plots...")

    plot_distance_kde(data, output_dir / "comparison_distance_kde.png")
    plot_excimer_heatmap(data, output_dir / "comparison_excimer_heatmap.png")
    plot_r_group_comparison(data, output_dir / "comparison_r_group_bars.png")

    # Distance correction plot (needs MACE data with MOE distances)
    if "MOE+MACE-OFF23" in data and "MOE/Amber:EHT" in data:
        plot_distance_correction(
            data["MOE/Amber:EHT"],
            data["MOE+MACE-OFF23"],
            output_dir / "comparison_distance_correction.png",
        )

    print(f"\nAll outputs saved to {output_dir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
