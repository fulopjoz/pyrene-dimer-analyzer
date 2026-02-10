"""Structure-Activity Relationship (SAR) Analysis of Screening Results.

Correlates substituent descriptors (Hammett sigma, Taft Es, Charton v, MR)
with excimer formation metrics from the binaphthalene dimer screening.

Usage:
    python sar_analysis.py
    python sar_analysis.py --output-dir plots/sar
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


# ============================================================
# Substituent Descriptors (Literature Values)
# ============================================================
# Sources:
#   Hammett sigma_para, sigma_meta: Hansch, Leo & Taft (1991)
#   Taft Es (steric): Taft (1956), updated by MacPhee (1978)
#   Charton v (steric): Charton (1975)
#   Molar Refractivity (MR): Hansch & Leo (1979)
#   Swain-Lupton F (field), R (resonance): Swain & Lupton (1968)

SUBSTITUENT_DESCRIPTORS = {
    "H":    {"sigma_para": 0.00, "sigma_meta": 0.00, "Es": 0.00,  "v": 0.00, "MR": 0.10,  "F": 0.00,  "R": 0.00,  "type": "reference"},
    "Me":   {"sigma_para": -0.17, "sigma_meta": -0.07, "Es": -1.24, "v": 0.52, "MR": 5.65,  "F": 0.01,  "R": -0.18, "type": "alkyl_small"},
    "Et":   {"sigma_para": -0.15, "sigma_meta": -0.07, "Es": -1.31, "v": 0.56, "MR": 10.30, "F": -0.01, "R": -0.15, "type": "alkyl_small"},
    "nPr":  {"sigma_para": -0.13, "sigma_meta": -0.06, "Es": -1.60, "v": 0.68, "MR": 14.93, "F": -0.01, "R": -0.12, "type": "alkyl_medium"},
    "iPr":  {"sigma_para": -0.15, "sigma_meta": -0.07, "Es": -1.71, "v": 0.76, "MR": 14.93, "F": -0.01, "R": -0.15, "type": "alkyl_medium"},
    "nBu":  {"sigma_para": -0.16, "sigma_meta": -0.07, "Es": -1.63, "v": 0.68, "MR": 19.58, "F": -0.01, "R": -0.16, "type": "alkyl_medium"},
    "tBu":  {"sigma_para": -0.20, "sigma_meta": -0.10, "Es": -2.78, "v": 1.24, "MR": 19.62, "F": -0.01, "R": -0.20, "type": "alkyl_bulky"},
    "cHex": {"sigma_para": -0.15, "sigma_meta": -0.05, "Es": -2.03, "v": 0.87, "MR": 25.36, "F": -0.02, "R": -0.13, "type": "alkyl_bulky"},
    "MeO":  {"sigma_para": -0.27, "sigma_meta": 0.12,  "Es": -0.55, "v": 0.36, "MR": 7.87,  "F": 0.29,  "R": -0.56, "type": "donor"},
    "OEt":  {"sigma_para": -0.24, "sigma_meta": 0.10,  "Es": -0.55, "v": 0.36, "MR": 12.47, "F": 0.22,  "R": -0.46, "type": "donor"},
    "F":    {"sigma_para": 0.06,  "sigma_meta": 0.34,  "Es": -0.46, "v": 0.27, "MR": 0.92,  "F": 0.45,  "R": -0.39, "type": "halogen"},
    "Cl":   {"sigma_para": 0.23,  "sigma_meta": 0.37,  "Es": -0.97, "v": 0.55, "MR": 6.03,  "F": 0.42,  "R": -0.19, "type": "halogen"},
    "CF3":  {"sigma_para": 0.54,  "sigma_meta": 0.43,  "Es": -2.40, "v": 0.91, "MR": 5.02,  "F": 0.38,  "R": 0.16,  "type": "ewg"},
    "CN":   {"sigma_para": 0.66,  "sigma_meta": 0.56,  "Es": -0.51, "v": 0.40, "MR": 6.33,  "F": 0.51,  "R": 0.15,  "type": "ewg"},
    "NMe2": {"sigma_para": -0.83, "sigma_meta": -0.16, "Es": -0.93, "v": 0.61, "MR": 15.55, "F": 0.15,  "R": -0.98, "type": "donor"},
    "Ph":   {"sigma_para": -0.01, "sigma_meta": 0.06,  "Es": -2.55, "v": 1.66, "MR": 25.36, "F": 0.12,  "R": -0.13, "type": "aromatic"},
}


def load_screening_data():
    """Load and merge screening results with substituent descriptors."""
    # Load summary
    summary = pd.read_csv("binaph_screening_summary.csv")

    # Parse substituent name (format: RGroup_ScreenGroup)
    summary["r_group"] = summary["substituent"].str.split("_").str[0:2].str.join("_")
    summary["screen_group"] = summary["substituent"].str.split("_").str[-1]

    # Also extract the luminescent R-group prefix
    lum_map = {}
    for name in summary["substituent"]:
        parts = name.split("_")
        if len(parts) >= 2:
            # Find the screen group (last part)
            for sg in SUBSTITUENT_DESCRIPTORS:
                if name.endswith(f"_{sg}"):
                    lum = name[: -(len(sg) + 1)]
                    lum_map[name] = (lum, sg)
                    break

    summary["lum_r_group"] = summary["substituent"].map(
        lambda x: lum_map.get(x, ("?", "?"))[0]
    )
    summary["screen_group"] = summary["substituent"].map(
        lambda x: lum_map.get(x, ("?", "?"))[1]
    )

    # Merge substituent descriptors
    desc_df = pd.DataFrame(SUBSTITUENT_DESCRIPTORS).T
    desc_df.index.name = "screen_group"
    desc_df = desc_df.reset_index()

    merged = summary.merge(desc_df, on="screen_group", how="left")

    # Ensure numeric columns are float (merge can produce object dtype)
    numeric_cols = ["sigma_para", "sigma_meta", "Es", "v", "MR", "F", "R",
                    "excimer_fraction", "mean_angle", "mean_distance",
                    "mean_overlap", "best_overlap"]
    for col in numeric_cols:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")

    return merged


def compute_correlations(df):
    """Compute Spearman correlations between descriptors and excimer metrics."""
    desc_cols = ["sigma_para", "sigma_meta", "Es", "v", "MR", "F", "R"]
    metric_cols = ["excimer_fraction", "mean_angle", "mean_distance",
                   "mean_overlap", "best_overlap"]

    results = []
    for desc in desc_cols:
        for metric in metric_cols:
            valid = df[[desc, metric]].dropna()
            if len(valid) < 5:
                continue
            rho, p = stats.spearmanr(valid[desc], valid[metric])
            r_pearson, p_pearson = stats.pearsonr(valid[desc], valid[metric])
            results.append({
                "descriptor": desc,
                "metric": metric,
                "spearman_rho": round(rho, 3),
                "spearman_p": round(p, 4),
                "pearson_r": round(r_pearson, 3),
                "pearson_p": round(p_pearson, 4),
                "n": len(valid),
            })

    return pd.DataFrame(results)


def plot_heatmap(df, output_dir):
    """4 x 16 heatmap of excimer fraction: R-groups vs screening groups."""
    pivot = df.pivot_table(
        values="excimer_fraction",
        index="lum_r_group",
        columns="screen_group",
        aggfunc="first",
    )

    # Order columns by substituent type
    col_order = ["H", "Me", "Et", "nPr", "iPr", "nBu", "tBu", "cHex",
                 "MeO", "OEt", "F", "Cl", "CF3", "CN", "NMe2", "Ph"]
    col_order = [c for c in col_order if c in pivot.columns]
    pivot = pivot[col_order]

    # Order rows
    row_order = ["Pyr", "EtynPyr", "DCV_Th", "CNPh_Th"]
    row_order = [r for r in row_order if r in pivot.index]
    pivot = pivot.loc[row_order]

    fig, ax = plt.subplots(figsize=(14, 4))
    sns.heatmap(
        pivot,
        annot=True,
        fmt=".0%",
        cmap="RdYlGn",
        vmin=0, vmax=0.5,
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": "Excimer Fraction"},
    )
    ax.set_title("Excimer Fraction: Luminescent R-Group vs Screening Substituent")
    ax.set_ylabel("Luminescent R-Group")
    ax.set_xlabel("Screening Substituent")
    plt.tight_layout()
    fig.savefig(output_dir / "heatmap_excimer_fraction.png", dpi=200)
    plt.close(fig)
    print(f"Saved heatmap to {output_dir / 'heatmap_excimer_fraction.png'}")


def plot_r_group_comparison(df, output_dir):
    """Bar chart comparing excimer fraction across R-groups."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: Mean excimer fraction by R-group
    r_stats = df.groupby("lum_r_group")["excimer_fraction"].agg(["mean", "std", "count"])
    r_order = ["Pyr", "EtynPyr", "DCV_Th", "CNPh_Th"]
    r_stats = r_stats.reindex([r for r in r_order if r in r_stats.index])

    colors = ["#2196F3", "#FF9800", "#4CAF50", "#9C27B0"]
    axes[0].bar(range(len(r_stats)), r_stats["mean"], yerr=r_stats["std"],
                color=colors[:len(r_stats)], capsize=5, edgecolor="black", linewidth=0.5)
    axes[0].set_xticks(range(len(r_stats)))
    axes[0].set_xticklabels(r_stats.index, fontsize=10)
    axes[0].set_ylabel("Mean Excimer Fraction")
    axes[0].set_title("Excimer Fraction by R-Group")
    axes[0].set_ylim(0, max(0.3, r_stats["mean"].max() * 1.5))

    # Add n variants with excimer
    r_exc = df[df["excimer_fraction"] > 0].groupby("lum_r_group").size()
    r_total = df.groupby("lum_r_group").size()
    for i, rg in enumerate(r_stats.index):
        n_exc = r_exc.get(rg, 0)
        n_tot = r_total.get(rg, 0)
        axes[0].text(i, r_stats.loc[rg, "mean"] + r_stats.loc[rg, "std"] + 0.01,
                     f"{n_exc}/{n_tot}", ha="center", fontsize=9)

    # Panel 2: Hit rate (fraction of variants showing any excimer)
    hit_rates = []
    for rg in r_stats.index:
        rg_df = df[df["lum_r_group"] == rg]
        n_hit = sum(rg_df["excimer_fraction"] > 0)
        hit_rates.append(n_hit / len(rg_df) if len(rg_df) > 0 else 0)

    axes[1].bar(range(len(r_stats)), hit_rates, color=colors[:len(r_stats)],
                edgecolor="black", linewidth=0.5)
    axes[1].set_xticks(range(len(r_stats)))
    axes[1].set_xticklabels(r_stats.index, fontsize=10)
    axes[1].set_ylabel("Fraction of Variants with Excimer")
    axes[1].set_title("Excimer Hit Rate by R-Group")
    axes[1].set_ylim(0, 0.6)

    plt.tight_layout()
    fig.savefig(output_dir / "r_group_comparison.png", dpi=200)
    plt.close(fig)
    print(f"Saved R-group comparison to {output_dir / 'r_group_comparison.png'}")


def plot_sar_correlations(df, output_dir):
    """Scatter plots of excimer fraction vs substituent descriptors."""
    desc_cols = [
        ("Es", "Taft Steric Parameter (Es)"),
        ("sigma_para", "Hammett $\\sigma_{para}$"),
        ("v", "Charton Steric Parameter (v)"),
        ("MR", "Molar Refractivity (MR)"),
    ]

    r_groups = df["lum_r_group"].unique()
    colors = {"Pyr": "#2196F3", "EtynPyr": "#FF9800",
              "DCV_Th": "#4CAF50", "CNPh_Th": "#9C27B0"}

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for ax, (desc, label) in zip(axes.flat, desc_cols):
        for rg in r_groups:
            rg_df = df[df["lum_r_group"] == rg]
            ax.scatter(rg_df[desc], rg_df["excimer_fraction"],
                       color=colors.get(rg, "gray"), label=rg,
                       s=40, alpha=0.7, edgecolors="black", linewidths=0.3)

        # Overall trend line
        valid = df[[desc, "excimer_fraction"]].dropna()
        if len(valid) > 5:
            rho, p = stats.spearmanr(valid[desc], valid["excimer_fraction"])
            slope, intercept, r, p_lr, se = stats.linregress(
                valid[desc], valid["excimer_fraction"]
            )
            x_range = np.linspace(valid[desc].min(), valid[desc].max(), 100)
            ax.plot(x_range, slope * x_range + intercept, "k--", alpha=0.5, linewidth=1)
            ax.set_title(f"{label}\n$\\rho$={rho:.2f}, p={p:.3f}")
        else:
            ax.set_title(label)

        ax.set_xlabel(label)
        ax.set_ylabel("Excimer Fraction")
        ax.set_ylim(-0.03, max(0.5, df["excimer_fraction"].max() * 1.2))

    axes[0, 0].legend(fontsize=8, loc="upper right")
    plt.tight_layout()
    fig.savefig(output_dir / "sar_correlations.png", dpi=200)
    plt.close(fig)
    print(f"Saved SAR correlations to {output_dir / 'sar_correlations.png'}")


def plot_geometric_distributions(df, output_dir):
    """Box plots of geometric parameters by R-group and screening group type."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    metrics = [
        ("mean_angle", "Mean Plane Angle (deg)"),
        ("mean_distance", "Mean Inter-plane Distance (A)"),
        ("mean_overlap", "Mean Pi Overlap (%)"),
    ]

    for ax, (col, label) in zip(axes, metrics):
        r_order = ["Pyr", "EtynPyr", "DCV_Th", "CNPh_Th"]
        data_frames = []
        for rg in r_order:
            rg_df = df[df["lum_r_group"] == rg]
            if not rg_df.empty:
                data_frames.append(rg_df)

        if data_frames:
            plot_df = pd.concat(data_frames)
            sns.boxplot(data=plot_df, x="lum_r_group", y=col,
                        order=[r for r in r_order if r in plot_df["lum_r_group"].values],
                        ax=ax, palette=["#2196F3", "#FF9800", "#4CAF50", "#9C27B0"])
        ax.set_xlabel("R-Group")
        ax.set_ylabel(label)
        ax.set_title(label)

    plt.tight_layout()
    fig.savefig(output_dir / "geometric_distributions.png", dpi=200)
    plt.close(fig)
    print(f"Saved geometric distributions to {output_dir / 'geometric_distributions.png'}")


def plot_substituent_type_analysis(df, output_dir):
    """Excimer fraction by substituent electronic/steric type."""
    type_order = ["reference", "alkyl_small", "alkyl_medium", "alkyl_bulky",
                  "donor", "halogen", "ewg", "aromatic"]
    type_labels = {"reference": "H", "alkyl_small": "Small Alkyl",
                   "alkyl_medium": "Med Alkyl", "alkyl_bulky": "Bulky Alkyl",
                   "donor": "e-Donor", "halogen": "Halogen",
                   "ewg": "e-Withdrawing", "aromatic": "Aromatic"}

    fig, ax = plt.subplots(figsize=(10, 5))

    present_types = [t for t in type_order if t in df["type"].values]
    type_stats = df.groupby("type")["excimer_fraction"].agg(["mean", "std", "count"])
    type_stats = type_stats.reindex(present_types)

    colors = {
        "reference": "#9E9E9E", "alkyl_small": "#81C784",
        "alkyl_medium": "#4CAF50", "alkyl_bulky": "#2E7D32",
        "donor": "#42A5F5", "halogen": "#EF5350",
        "ewg": "#FFA726", "aromatic": "#AB47BC",
    }

    bars = ax.bar(
        range(len(type_stats)),
        type_stats["mean"],
        yerr=type_stats["std"],
        color=[colors.get(t, "gray") for t in type_stats.index],
        capsize=4, edgecolor="black", linewidth=0.5,
    )
    ax.set_xticks(range(len(type_stats)))
    ax.set_xticklabels([type_labels.get(t, t) for t in type_stats.index],
                       rotation=30, ha="right")
    ax.set_ylabel("Mean Excimer Fraction")
    ax.set_title("Excimer Fraction by Substituent Electronic/Steric Type")

    # Add count annotations
    for i, (t, row) in enumerate(type_stats.iterrows()):
        ax.text(i, row["mean"] + row["std"] + 0.005,
                f"n={int(row['count'])}", ha="center", fontsize=8)

    plt.tight_layout()
    fig.savefig(output_dir / "substituent_type_analysis.png", dpi=200)
    plt.close(fig)
    print(f"Saved substituent type analysis to {output_dir / 'substituent_type_analysis.png'}")


def main():
    parser = argparse.ArgumentParser(description="SAR Analysis of Screening Results")
    parser.add_argument("--output-dir", type=str, default="plots/sar",
                        help="Output directory for plots")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("SAR ANALYSIS OF BINAPHTHALENE DIMER SCREENING")
    print("=" * 70)

    # Load data
    df = load_screening_data()
    print(f"Loaded {len(df)} variants")
    print(f"R-groups: {sorted(df['lum_r_group'].unique())}")
    print(f"Screen groups: {sorted(df['screen_group'].unique())}")
    print()

    # Correlation analysis
    print("--- CORRELATION ANALYSIS ---")
    corr_df = compute_correlations(df)
    exc_corrs = corr_df[corr_df["metric"] == "excimer_fraction"].sort_values(
        "spearman_rho", key=abs, ascending=False
    )
    print("\nCorrelations with excimer_fraction (sorted by |rho|):")
    print(exc_corrs.to_string(index=False))

    # Per-R-group correlations
    print("\n--- PER-R-GROUP CORRELATIONS ---")
    for rg in sorted(df["lum_r_group"].unique()):
        rg_df = df[df["lum_r_group"] == rg]
        rg_corrs = compute_correlations(rg_df)
        rg_exc = rg_corrs[rg_corrs["metric"] == "excimer_fraction"]
        sig = rg_exc[rg_exc["spearman_p"] < 0.1]
        if not sig.empty:
            print(f"\n{rg} (significant correlations with excimer_fraction):")
            print(sig.to_string(index=False))

    # Summary statistics by R-group
    print("\n--- R-GROUP SUMMARY ---")
    for rg in ["Pyr", "EtynPyr", "DCV_Th", "CNPh_Th"]:
        rg_df = df[df["lum_r_group"] == rg]
        n_exc = sum(rg_df["excimer_fraction"] > 0)
        mean_exc = rg_df["excimer_fraction"].mean()
        max_exc = rg_df["excimer_fraction"].max()
        best = rg_df.loc[rg_df["excimer_fraction"].idxmax(), "substituent"] if max_exc > 0 else "none"
        print(f"  {rg:12s}: {n_exc:>2d}/16 hit, mean={mean_exc:.1%}, "
              f"max={max_exc:.1%} ({best})")

    # Save correlation table
    corr_df.to_csv(output_dir / "correlation_analysis.csv", index=False)
    print(f"\nSaved correlation analysis to {output_dir / 'correlation_analysis.csv'}")

    # Save merged data with descriptors
    df.to_csv(output_dir / "screening_with_descriptors.csv", index=False)
    print(f"Saved merged data to {output_dir / 'screening_with_descriptors.csv'}")

    # Generate plots
    print("\n--- GENERATING PLOTS ---")
    plot_heatmap(df, output_dir)
    plot_r_group_comparison(df, output_dir)
    plot_sar_correlations(df, output_dir)
    plot_geometric_distributions(df, output_dir)
    plot_substituent_type_analysis(df, output_dir)

    print("\nDone!")


if __name__ == "__main__":
    main()
