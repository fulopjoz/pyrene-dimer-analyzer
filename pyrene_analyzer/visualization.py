"""
Visualization Module
====================

Functions for creating publication-quality plots of pyrene dimer analysis results.

This module provides various plotting functions for visualizing:
- Angle vs energy relationships
- Distance vs overlap correlations
- Conformer distributions
- Variant comparisons

All plots are designed to be publication-ready with proper labels,
legends, and formatting.
"""

from pathlib import Path
from typing import List, Optional, Tuple, Union

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Set default style
plt.style.use("seaborn-v0_8-whitegrid")
sns.set_palette("husl")


def _setup_figure(
    figsize: Tuple[float, float] = (8, 6), dpi: int = 150
) -> Tuple[plt.Figure, plt.Axes]:
    """Create a figure with standard settings."""
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    return fig, ax


def _save_figure(
    fig: plt.Figure, output: Optional[Union[str, Path]], bbox_inches: str = "tight"
) -> None:
    """Save figure if output path is provided."""
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, bbox_inches=bbox_inches, facecolor="white")


def plot_angle_vs_energy(
    df: pd.DataFrame,
    output: Optional[Union[str, Path]] = None,
    title: Optional[str] = None,
    color_by: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 7),
    show_excimer_threshold: bool = True,
    alpha: float = 0.7,
    marker_size: int = 60,
) -> plt.Figure:
    """
    Create a scatter plot of plane angle (θ) vs relative energy.

    This visualization helps identify the relationship between
    conformational geometry and energetic stability.

    Args:
        df: DataFrame with analysis results.
        output: Optional path to save the figure.
        title: Plot title. If None, uses default.
        color_by: Column name to use for coloring points (e.g., 'molecule').
        figsize: Figure size in inches.
        show_excimer_threshold: If True, show vertical lines at excimer thresholds.
        alpha: Point transparency.
        marker_size: Size of scatter points.

    Returns:
        Matplotlib Figure object.

    Example:
        >>> fig = plot_angle_vs_energy(results, output='angle_energy.png')

    Notes:
        Excimer formation thresholds:
        - θ < 20°: Strong excimer (green region)
        - θ = 20-60°: Weak excimer (yellow region)
        - θ > 60°: Monomer (red region)
    """
    fig, ax = _setup_figure(figsize=figsize)

    # Determine energy column
    energy_col = "rel_energy_kcal_mol"
    if energy_col not in df.columns:
        energy_col = "energy_kcal_mol"

    if energy_col not in df.columns or df[energy_col].isna().all():
        # Create dummy energy if not available
        df = df.copy()
        df["_energy"] = range(len(df))
        energy_col = "_energy"

    # Create scatter plot
    if color_by and color_by in df.columns:
        unique_values = df[color_by].unique()
        colors = sns.color_palette("husl", len(unique_values))
        color_map = dict(zip(unique_values, colors))

        for value in unique_values:
            mask = df[color_by] == value
            ax.scatter(
                df.loc[mask, "plane_angle_deg"],
                df.loc[mask, energy_col],
                c=[color_map[value]],
                label=value,
                alpha=alpha,
                s=marker_size,
                edgecolors="white",
                linewidths=0.5,
            )
        ax.legend(title=color_by.replace("_", " ").title(), loc="best")
    else:
        scatter = ax.scatter(
            df["plane_angle_deg"],
            df[energy_col],
            c=df["pi_overlap_pct"],
            cmap="viridis",
            alpha=alpha,
            s=marker_size,
            edgecolors="white",
            linewidths=0.5,
        )
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("π-overlap (%)", fontsize=11)

    # Add excimer threshold lines
    if show_excimer_threshold:
        ax.axvline(
            x=20,
            color="green",
            linestyle="--",
            alpha=0.5,
            linewidth=2,
            label="Strong excimer threshold",
        )
        ax.axvline(
            x=60,
            color="orange",
            linestyle="--",
            alpha=0.5,
            linewidth=2,
            label="Weak excimer threshold",
        )

        # Add shaded regions
        ylim = ax.get_ylim()
        ax.axvspan(0, 20, alpha=0.1, color="green")
        ax.axvspan(20, 60, alpha=0.1, color="yellow")
        ax.axvspan(60, 90, alpha=0.1, color="red")
        ax.set_ylim(ylim)

    # Labels and title
    ax.set_xlabel("Plane-Plane Angle θ (degrees)", fontsize=12)
    ylabel = "Relative Energy (kcal/mol)" if "rel" in energy_col else "Conformer Index"
    ax.set_ylabel(ylabel, fontsize=12)

    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(
            "Plane Angle vs Energy for Pyrene Dimer Conformers",
            fontsize=14,
            fontweight="bold",
        )

    ax.set_xlim(0, 90)

    plt.tight_layout()
    _save_figure(fig, output)

    return fig


def plot_distance_vs_overlap(
    df: pd.DataFrame,
    output: Optional[Union[str, Path]] = None,
    title: Optional[str] = None,
    color_by: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 7),
    show_optimal_region: bool = True,
    alpha: float = 0.7,
    marker_size: int = 60,
) -> plt.Figure:
    """
    Create a scatter plot of inter-plane distance vs π-overlap.

    This visualization shows the relationship between stacking distance
    and orbital overlap, which are both critical for excimer formation.

    Args:
        df: DataFrame with analysis results.
        output: Optional path to save the figure.
        title: Plot title.
        color_by: Column name to use for coloring points.
        figsize: Figure size in inches.
        show_optimal_region: If True, highlight optimal excimer region.
        alpha: Point transparency.
        marker_size: Size of scatter points.

    Returns:
        Matplotlib Figure object.

    Example:
        >>> fig = plot_distance_vs_overlap(results, color_by='molecule')

    Notes:
        Optimal excimer formation region:
        - Distance: 3.3-3.7 Å
        - Overlap: > 70%
    """
    fig, ax = _setup_figure(figsize=figsize)

    # Create scatter plot
    if color_by and color_by in df.columns:
        unique_values = df[color_by].unique()
        colors = sns.color_palette("husl", len(unique_values))
        color_map = dict(zip(unique_values, colors))

        for value in unique_values:
            mask = df[color_by] == value
            ax.scatter(
                df.loc[mask, "interplane_distance_A"],
                df.loc[mask, "pi_overlap_pct"],
                c=[color_map[value]],
                label=value,
                alpha=alpha,
                s=marker_size,
                edgecolors="white",
                linewidths=0.5,
            )
        ax.legend(title=color_by.replace("_", " ").title(), loc="best")
    else:
        scatter = ax.scatter(
            df["interplane_distance_A"],
            df["pi_overlap_pct"],
            c=df["plane_angle_deg"],
            cmap="coolwarm",
            alpha=alpha,
            s=marker_size,
            edgecolors="white",
            linewidths=0.5,
        )
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Plane Angle θ (degrees)", fontsize=11)

    # Highlight optimal excimer region
    if show_optimal_region:
        rect = mpatches.Rectangle(
            (3.3, 70),
            0.4,
            30,
            linewidth=2,
            edgecolor="green",
            facecolor="green",
            alpha=0.2,
            linestyle="--",
            label="Optimal excimer region",
        )
        ax.add_patch(rect)
        ax.legend(loc="best")

    # Labels and title
    ax.set_xlabel("Inter-plane Distance d (Å)", fontsize=12)
    ax.set_ylabel("π-Overlap (%)", fontsize=12)

    if title:
        ax.set_title(title, fontsize=14, fontweight="bold")
    else:
        ax.set_title(
            "Inter-plane Distance vs π-Overlap", fontsize=14, fontweight="bold"
        )

    ax.set_ylim(0, 100)

    plt.tight_layout()
    _save_figure(fig, output)

    return fig


def plot_conformer_distribution(
    df: pd.DataFrame,
    output: Optional[Union[str, Path]] = None,
    parameters: Optional[List[str]] = None,
    figsize: Tuple[float, float] = (14, 10),
    bins: int = 20,
) -> plt.Figure:
    """
    Create histograms showing the distribution of geometric parameters.

    Args:
        df: DataFrame with analysis results.
        output: Optional path to save the figure.
        parameters: List of parameters to plot. If None, uses defaults.
        figsize: Figure size in inches.
        bins: Number of histogram bins.

    Returns:
        Matplotlib Figure object.

    Example:
        >>> fig = plot_conformer_distribution(results)
    """
    if parameters is None:
        parameters = [
            "plane_angle_deg",
            "interplane_distance_A",
            "pi_overlap_pct",
            "centroid_distance_A",
        ]

    # Filter to available parameters
    parameters = [p for p in parameters if p in df.columns]
    n_params = len(parameters)

    if n_params == 0:
        raise ValueError("No valid parameters found in DataFrame")

    # Create subplot grid
    n_cols = min(2, n_params)
    n_rows = (n_params + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_params == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # Parameter labels
    labels = {
        "plane_angle_deg": "Plane Angle θ (degrees)",
        "interplane_distance_A": "Inter-plane Distance d (Å)",
        "pi_overlap_pct": "π-Overlap (%)",
        "centroid_distance_A": "Centroid Distance (Å)",
        "slip_stack_A": "Slip-Stack Displacement (Å)",
        "bridge_dihedral_L_deg": "Bridge Dihedral φL (degrees)",
        "bridge_dihedral_R_deg": "Bridge Dihedral φR (degrees)",
        "rel_energy_kcal_mol": "Relative Energy (kcal/mol)",
    }

    for idx, param in enumerate(parameters):
        ax = axes[idx]
        data = df[param].dropna()

        ax.hist(
            data,
            bins=bins,
            edgecolor="white",
            alpha=0.8,
            color=sns.color_palette("husl")[idx % 10],
        )

        ax.set_xlabel(labels.get(param, param), fontsize=11)
        ax.set_ylabel("Count", fontsize=11)
        ax.set_title(
            f"Distribution of {labels.get(param, param)}",
            fontsize=12,
            fontweight="bold",
        )

        # Add statistics
        mean_val = data.mean()
        ax.axvline(
            mean_val,
            color="red",
            linestyle="--",
            linewidth=2,
            label=f"Mean: {mean_val:.2f}",
        )
        ax.legend(loc="best")

    # Hide unused axes
    for idx in range(n_params, len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle("Conformer Parameter Distributions", fontsize=14, fontweight="bold")
    plt.tight_layout()
    _save_figure(fig, output)

    return fig


def plot_variant_comparison(
    df: pd.DataFrame,
    variants: Optional[List[str]] = None,
    output: Optional[Union[str, Path]] = None,
    parameter: str = "plane_angle_deg",
    figsize: Tuple[float, float] = (12, 6),
) -> plt.Figure:
    """
    Create box plots comparing geometric parameters across R-group variants.

    Args:
        df: DataFrame with analysis results.
        variants: List of variant names to compare. If None, uses all unique molecules.
        output: Optional path to save the figure.
        parameter: Parameter to compare.
        figsize: Figure size in inches.

    Returns:
        Matplotlib Figure object.

    Example:
        >>> fig = plot_variant_comparison(
        ...     results,
        ...     variants=['Et', 'iPr', 'cHex', 'tBu'],
        ...     parameter='plane_angle_deg'
        ... )
    """
    fig, ax = _setup_figure(figsize=figsize)

    if variants is None:
        variants = df["molecule"].unique().tolist()

    # Filter data for specified variants
    plot_data = df[df["molecule"].isin(variants)].copy()

    if plot_data.empty:
        raise ValueError("No data found for specified variants")

    # Create box plot
    sns.boxplot(
        data=plot_data,
        x="molecule",
        y=parameter,
        hue="molecule",
        order=variants,
        palette="husl",
        ax=ax,
        legend=False,
    )

    # Add individual points
    sns.stripplot(
        data=plot_data,
        x="molecule",
        y=parameter,
        order=variants,
        color="black",
        alpha=0.3,
        size=4,
        ax=ax,
    )

    # Parameter labels
    labels = {
        "plane_angle_deg": "Plane Angle θ (degrees)",
        "interplane_distance_A": "Inter-plane Distance d (Å)",
        "pi_overlap_pct": "π-Overlap (%)",
        "centroid_distance_A": "Centroid Distance (Å)",
        "rel_energy_kcal_mol": "Relative Energy (kcal/mol)",
    }

    ax.set_xlabel("R-Group Variant", fontsize=12)
    ax.set_ylabel(labels.get(parameter, parameter), fontsize=12)
    ax.set_title(
        f"Comparison of {labels.get(parameter, parameter)} Across Variants",
        fontsize=14,
        fontweight="bold",
    )

    # Add excimer threshold for angle
    if parameter == "plane_angle_deg":
        ax.axhline(
            y=20,
            color="green",
            linestyle="--",
            alpha=0.5,
            label="Strong excimer threshold",
        )
        ax.axhline(
            y=60,
            color="orange",
            linestyle="--",
            alpha=0.5,
            label="Weak excimer threshold",
        )
        ax.legend(loc="best")

    plt.tight_layout()
    _save_figure(fig, output)

    return fig


def plot_energy_landscape(
    df: pd.DataFrame,
    output: Optional[Union[str, Path]] = None,
    figsize: Tuple[float, float] = (12, 8),
) -> plt.Figure:
    """
    Create a 2D energy landscape plot showing angle and distance vs energy.

    Args:
        df: DataFrame with analysis results.
        output: Optional path to save the figure.
        figsize: Figure size in inches.

    Returns:
        Matplotlib Figure object.
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    energy_col = "rel_energy_kcal_mol"
    if energy_col not in df.columns or df[energy_col].isna().all():
        energy_col = "energy_kcal_mol"

    if energy_col not in df.columns or df[energy_col].isna().all():
        raise ValueError("No energy data available")

    # Angle vs Energy
    ax1 = axes[0]
    scatter1 = ax1.scatter(
        df["plane_angle_deg"],
        df[energy_col],
        c=df["pi_overlap_pct"],
        cmap="viridis",
        alpha=0.7,
        s=50,
        edgecolors="white",
        linewidths=0.5,
    )
    ax1.set_xlabel("Plane Angle θ (degrees)", fontsize=11)
    ax1.set_ylabel("Relative Energy (kcal/mol)", fontsize=11)
    ax1.set_title("Angle vs Energy", fontsize=12, fontweight="bold")
    cbar1 = plt.colorbar(scatter1, ax=ax1)
    cbar1.set_label("π-overlap (%)")

    # Distance vs Energy
    ax2 = axes[1]
    scatter2 = ax2.scatter(
        df["interplane_distance_A"],
        df[energy_col],
        c=df["plane_angle_deg"],
        cmap="coolwarm",
        alpha=0.7,
        s=50,
        edgecolors="white",
        linewidths=0.5,
    )
    ax2.set_xlabel("Inter-plane Distance d (Å)", fontsize=11)
    ax2.set_ylabel("Relative Energy (kcal/mol)", fontsize=11)
    ax2.set_title("Distance vs Energy", fontsize=12, fontweight="bold")
    cbar2 = plt.colorbar(scatter2, ax=ax2)
    cbar2.set_label("Plane Angle θ (degrees)")

    plt.suptitle("Energy Landscape Analysis", fontsize=14, fontweight="bold")
    plt.tight_layout()
    _save_figure(fig, output)

    return fig


def plot_correlation_matrix(
    df: pd.DataFrame,
    output: Optional[Union[str, Path]] = None,
    figsize: Tuple[float, float] = (10, 8),
) -> plt.Figure:
    """
    Create a correlation matrix heatmap of geometric parameters.

    Args:
        df: DataFrame with analysis results.
        output: Optional path to save the figure.
        figsize: Figure size in inches.

    Returns:
        Matplotlib Figure object.
    """
    fig, ax = _setup_figure(figsize=figsize)

    # Select numeric columns
    numeric_cols = [
        "plane_angle_deg",
        "interplane_distance_A",
        "pi_overlap_pct",
        "centroid_distance_A",
        "slip_stack_A",
        "rel_energy_kcal_mol",
    ]
    numeric_cols = [c for c in numeric_cols if c in df.columns]

    if len(numeric_cols) < 2:
        raise ValueError("Not enough numeric columns for correlation matrix")

    # Calculate correlation matrix
    corr_matrix = df[numeric_cols].corr()

    # Create heatmap
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)
    sns.heatmap(
        corr_matrix,
        mask=mask,
        annot=True,
        fmt=".2f",
        cmap="coolwarm",
        center=0,
        square=True,
        linewidths=0.5,
        ax=ax,
        cbar_kws={"shrink": 0.8},
    )

    ax.set_title(
        "Correlation Matrix of Geometric Parameters", fontsize=14, fontweight="bold"
    )

    plt.tight_layout()
    _save_figure(fig, output)

    return fig


def create_summary_figure(
    df: pd.DataFrame,
    output: Optional[Union[str, Path]] = None,
    figsize: Tuple[float, float] = (16, 12),
) -> plt.Figure:
    """
    Create a comprehensive summary figure with multiple panels.

    Args:
        df: DataFrame with analysis results.
        output: Optional path to save the figure.
        figsize: Figure size in inches.

    Returns:
        Matplotlib Figure object.
    """
    fig = plt.figure(figsize=figsize)

    # Create grid
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

    # Panel 1: Angle distribution
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(
        df["plane_angle_deg"].dropna(),
        bins=20,
        edgecolor="white",
        alpha=0.8,
        color=sns.color_palette("husl")[0],
    )
    ax1.set_xlabel("Plane Angle θ (degrees)")
    ax1.set_ylabel("Count")
    ax1.set_title("Angle Distribution")
    ax1.axvline(
        df["plane_angle_deg"].mean(),
        color="red",
        linestyle="--",
        label=f'Mean: {df["plane_angle_deg"].mean():.1f}°',
    )
    ax1.legend()

    # Panel 2: Distance distribution
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.hist(
        df["interplane_distance_A"].dropna(),
        bins=20,
        edgecolor="white",
        alpha=0.8,
        color=sns.color_palette("husl")[1],
    )
    ax2.set_xlabel("Inter-plane Distance d (Å)")
    ax2.set_ylabel("Count")
    ax2.set_title("Distance Distribution")
    ax2.axvline(
        df["interplane_distance_A"].mean(),
        color="red",
        linestyle="--",
        label=f'Mean: {df["interplane_distance_A"].mean():.2f} Å',
    )
    ax2.legend()

    # Panel 3: Overlap distribution
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.hist(
        df["pi_overlap_pct"].dropna(),
        bins=20,
        edgecolor="white",
        alpha=0.8,
        color=sns.color_palette("husl")[2],
    )
    ax3.set_xlabel("π-Overlap (%)")
    ax3.set_ylabel("Count")
    ax3.set_title("Overlap Distribution")
    ax3.axvline(
        df["pi_overlap_pct"].mean(),
        color="red",
        linestyle="--",
        label=f'Mean: {df["pi_overlap_pct"].mean():.1f}%',
    )
    ax3.legend()

    # Panel 4: Angle vs Distance scatter
    ax4 = fig.add_subplot(gs[1, 0])
    scatter = ax4.scatter(
        df["plane_angle_deg"],
        df["interplane_distance_A"],
        c=df["pi_overlap_pct"],
        cmap="viridis",
        alpha=0.7,
        s=40,
    )
    ax4.set_xlabel("Plane Angle θ (degrees)")
    ax4.set_ylabel("Inter-plane Distance d (Å)")
    ax4.set_title("Angle vs Distance")
    plt.colorbar(scatter, ax=ax4, label="π-overlap (%)")

    # Panel 5: Distance vs Overlap scatter
    ax5 = fig.add_subplot(gs[1, 1])
    scatter = ax5.scatter(
        df["interplane_distance_A"],
        df["pi_overlap_pct"],
        c=df["plane_angle_deg"],
        cmap="coolwarm",
        alpha=0.7,
        s=40,
    )
    ax5.set_xlabel("Inter-plane Distance d (Å)")
    ax5.set_ylabel("π-Overlap (%)")
    ax5.set_title("Distance vs Overlap")
    plt.colorbar(scatter, ax=ax5, label="Angle (°)")

    # Panel 6: Classification pie chart (if classification exists)
    ax6 = fig.add_subplot(gs[1, 2])
    if "classification" in df.columns:
        class_counts = df["classification"].value_counts()
        colors = {"strong_excimer": "green", "weak_excimer": "yellow", "monomer": "red"}
        ax6.pie(
            class_counts.values,
            labels=class_counts.index,
            autopct="%1.1f%%",
            colors=[colors.get(c, "gray") for c in class_counts.index],
        )
        ax6.set_title("Conformer Classification")
    else:
        # Show statistics table instead
        stats = df[
            ["plane_angle_deg", "interplane_distance_A", "pi_overlap_pct"]
        ].describe()
        ax6.axis("off")
        table = ax6.table(
            cellText=stats.round(2).values,
            rowLabels=stats.index,
            colLabels=["θ (°)", "d (Å)", "Overlap (%)"],
            loc="center",
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        ax6.set_title("Summary Statistics")

    plt.suptitle(
        "Pyrene Dimer Conformational Analysis Summary", fontsize=16, fontweight="bold"
    )

    _save_figure(fig, output)

    return fig
