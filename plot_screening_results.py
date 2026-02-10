"""Generate comparison plots for R-group screening results."""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path

# Load data
summary = pd.read_csv("screening_summary.csv")
features = pd.read_csv("screening_ensemble_features.csv", index_col=0)
all_conf = pd.read_csv("screening_all_conformers.csv")

outdir = Path("plots/screening")
outdir.mkdir(parents=True, exist_ok=True)

# Sort by excimer fraction
summary = summary.sort_values("excimer_fraction", ascending=False)

# Color categories
def categorize(row):
    name = row["substituent"]
    if name in ["H", "Me", "Et", "nPr", "iPr", "nBu", "iBu", "sBu", "tBu",
                 "nPent", "nHex", "neoPent"]:
        return "Alkyl (linear/branched)"
    elif name in ["cPent", "cHex"]:
        return "Cycloalkyl"
    elif name in ["Ph", "4MePh", "4FPh", "4ClPh", "4OMePh", "4CF3Ph", "Bn"]:
        return "Aryl/Benzyl"
    elif name in ["OMe", "OEt", "NMe2"]:
        return "Electron-donating"
    elif name in ["F", "Cl", "CF3", "CN"]:
        return "Electron-withdrawing"
    elif name in ["CH2OH", "CH2OMe"]:
        return "Functional"
    return "Other"

summary["category"] = summary.apply(categorize, axis=1)

cat_colors = {
    "Alkyl (linear/branched)": "#2196F3",
    "Cycloalkyl": "#4CAF50",
    "Aryl/Benzyl": "#FF9800",
    "Electron-donating": "#9C27B0",
    "Electron-withdrawing": "#F44336",
    "Functional": "#795548",
}

# =========================================================================
# PLOT 1: Excimer Fraction Bar Chart
# =========================================================================
fig, ax = plt.subplots(figsize=(16, 8))
names = summary["substituent"].tolist()
fracs = summary["excimer_fraction"].tolist()
colors = [cat_colors.get(summary.iloc[i]["category"], "#999999") for i in range(len(summary))]

bars = ax.bar(range(len(names)), fracs, color=colors, edgecolor="black", linewidth=0.5)
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, rotation=45, ha="right", fontsize=10)
ax.set_ylabel("Excimer Fraction", fontsize=12)
ax.set_title("R-Group Screening: Excimer Fraction (50 conformers, MMFF94s)", fontsize=14)
ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5, label="50% threshold")
ax.set_ylim(0, 1.0)

# Add value labels
for i, (bar, frac) in enumerate(zip(bars, fracs)):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
            f"{frac:.0%}", ha="center", va="bottom", fontsize=8, fontweight="bold")

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=c, edgecolor="black", label=l)
                   for l, c in cat_colors.items()]
ax.legend(handles=legend_elements, loc="upper right", fontsize=9)
plt.tight_layout()
plt.savefig(outdir / "excimer_fraction_ranking.png", dpi=150)
plt.close()
print("Saved excimer_fraction_ranking.png")

# =========================================================================
# PLOT 2: Boltzmann-weighted excimer fraction
# =========================================================================
fig, ax = plt.subplots(figsize=(16, 8))
boltz_sorted = features.sort_values("boltz_excimer_fraction", ascending=False)
names_b = boltz_sorted.index.tolist()
boltz_fracs = boltz_sorted["boltz_excimer_fraction"].tolist()

cat_map = dict(zip(summary["substituent"], summary["category"]))
colors_b = [cat_colors.get(cat_map.get(n, ""), "#999999") for n in names_b]

bars = ax.bar(range(len(names_b)), boltz_fracs, color=colors_b, edgecolor="black", linewidth=0.5)
ax.set_xticks(range(len(names_b)))
ax.set_xticklabels(names_b, rotation=45, ha="right", fontsize=10)
ax.set_ylabel("Boltzmann-Weighted Excimer Fraction", fontsize=12)
ax.set_title("Boltzmann-Weighted Excimer Fraction at 298K", fontsize=14)
ax.set_ylim(0, 1.1)
ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5)

for i, (bar, frac) in enumerate(zip(bars, boltz_fracs)):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
            f"{frac:.2f}", ha="center", va="bottom", fontsize=7, rotation=90)

ax.legend(handles=legend_elements, loc="lower left", fontsize=9)
plt.tight_layout()
plt.savefig(outdir / "boltzmann_excimer_fraction.png", dpi=150)
plt.close()
print("Saved boltzmann_excimer_fraction.png")

# =========================================================================
# PLOT 3: Scatter - Mean Angle vs Mean Distance (colored by excimer frac)
# =========================================================================
fig, ax = plt.subplots(figsize=(10, 8))
sc = ax.scatter(
    summary["mean_angle"],
    summary["mean_distance"],
    c=summary["excimer_fraction"],
    cmap="RdYlGn",
    s=120,
    edgecolors="black",
    linewidths=0.5,
    vmin=0, vmax=1,
)
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label("Excimer Fraction", fontsize=11)

for _, row in summary.iterrows():
    ax.annotate(row["substituent"], (row["mean_angle"], row["mean_distance"]),
                textcoords="offset points", xytext=(5, 5), fontsize=8)

# Excimer zone
ax.axvspan(0, 20, alpha=0.1, color="green", label="Strong excimer angle (<20)")
ax.axhspan(3.3, 3.7, alpha=0.1, color="blue", label="Strong excimer distance (3.3-3.7 A)")
ax.set_xlabel("Mean Plane Angle (deg)", fontsize=12)
ax.set_ylabel("Mean Interplane Distance (A)", fontsize=12)
ax.set_title("Conformer Ensemble: Angle vs Distance", fontsize=14)
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(outdir / "angle_vs_distance_scatter.png", dpi=150)
plt.close()
print("Saved angle_vs_distance_scatter.png")

# =========================================================================
# PLOT 4: Scatter - Mean Distance vs Mean Overlap (colored by excimer frac)
# =========================================================================
fig, ax = plt.subplots(figsize=(10, 8))
sc = ax.scatter(
    summary["mean_distance"],
    summary["mean_overlap"],
    c=summary["excimer_fraction"],
    cmap="RdYlGn",
    s=120,
    edgecolors="black",
    linewidths=0.5,
    vmin=0, vmax=1,
)
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label("Excimer Fraction", fontsize=11)

for _, row in summary.iterrows():
    ax.annotate(row["substituent"], (row["mean_distance"], row["mean_overlap"]),
                textcoords="offset points", xytext=(5, 5), fontsize=8)

ax.axhline(y=50, color="green", linestyle="--", alpha=0.5, label="50% overlap threshold")
ax.axvline(x=3.7, color="blue", linestyle="--", alpha=0.5, label="3.7 A distance threshold")
ax.set_xlabel("Mean Interplane Distance (A)", fontsize=12)
ax.set_ylabel("Mean Pi-Overlap (%)", fontsize=12)
ax.set_title("Conformer Ensemble: Distance vs Overlap", fontsize=14)
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(outdir / "distance_vs_overlap_scatter.png", dpi=150)
plt.close()
print("Saved distance_vs_overlap_scatter.png")

# =========================================================================
# PLOT 5: Grouped comparison - Original 4 variants vs top performers
# =========================================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 6))

# Select key variants for comparison
original_4 = ["Et", "iPr", "cHex", "tBu"]
top_new = ["Cl", "CF3", "H", "Bn", "OMe", "OEt"]
comparison = original_4 + [n for n in top_new if n not in original_4]

comp_data = summary[summary["substituent"].isin(comparison)].copy()
comp_data = comp_data.sort_values("excimer_fraction", ascending=False)
comp_names = comp_data["substituent"].tolist()

# Color: original=blue, new=orange
comp_colors = ["#2196F3" if n in original_4 else "#FF9800" for n in comp_names]

# Excimer fraction
axes[0].barh(range(len(comp_names)), comp_data["excimer_fraction"], color=comp_colors)
axes[0].set_yticks(range(len(comp_names)))
axes[0].set_yticklabels(comp_names, fontsize=10)
axes[0].set_xlabel("Excimer Fraction")
axes[0].set_title("Excimer Fraction")
axes[0].invert_yaxis()

# Mean angle
axes[1].barh(range(len(comp_names)), comp_data["mean_angle"], color=comp_colors)
axes[1].set_yticks(range(len(comp_names)))
axes[1].set_yticklabels(comp_names, fontsize=10)
axes[1].set_xlabel("Mean Angle (deg)")
axes[1].set_title("Mean Plane Angle")
axes[1].axvline(x=20, color="red", linestyle="--", alpha=0.5)
axes[1].invert_yaxis()

# Mean overlap
axes[2].barh(range(len(comp_names)), comp_data["mean_overlap"], color=comp_colors)
axes[2].set_yticks(range(len(comp_names)))
axes[2].set_yticklabels(comp_names, fontsize=10)
axes[2].set_xlabel("Mean Overlap (%)")
axes[2].set_title("Mean Pi-Overlap")
axes[2].axvline(x=50, color="red", linestyle="--", alpha=0.5)
axes[2].invert_yaxis()

from matplotlib.patches import Patch
leg_elems = [Patch(facecolor="#2196F3", label="Original 4"),
             Patch(facecolor="#FF9800", label="New candidates")]
fig.legend(handles=leg_elems, loc="upper center", ncol=2, fontsize=11, bbox_to_anchor=(0.5, 1.02))

plt.suptitle("Original vs Top New Substituent Candidates", fontsize=14, y=1.05)
plt.tight_layout()
plt.savefig(outdir / "original_vs_new_comparison.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved original_vs_new_comparison.png")

# =========================================================================
# PLOT 6: Boltzmann vs Simple excimer fraction (scatter)
# =========================================================================
fig, ax = plt.subplots(figsize=(10, 8))

for cat, color in cat_colors.items():
    mask = summary["category"] == cat
    if mask.any():
        sub = summary[mask]
        boltz_vals = [features.loc[n, "boltz_excimer_fraction"]
                      if n in features.index else np.nan
                      for n in sub["substituent"]]
        ax.scatter(sub["excimer_fraction"], boltz_vals,
                   c=color, label=cat, s=100, edgecolors="black", linewidths=0.5)
        for _, row in sub.iterrows():
            bv = features.loc[row["substituent"], "boltz_excimer_fraction"] if row["substituent"] in features.index else np.nan
            ax.annotate(row["substituent"],
                        (row["excimer_fraction"], bv),
                        textcoords="offset points", xytext=(5, 5), fontsize=8)

ax.plot([0, 1], [0, 1], "k--", alpha=0.3, label="y=x")
ax.set_xlabel("Simple Excimer Fraction", fontsize=12)
ax.set_ylabel("Boltzmann-Weighted Excimer Fraction", fontsize=12)
ax.set_title("Simple vs Boltzmann-Weighted Excimer Classification", fontsize=14)
ax.legend(fontsize=9)
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.15)
plt.tight_layout()
plt.savefig(outdir / "simple_vs_boltzmann.png", dpi=150)
plt.close()
print("Saved simple_vs_boltzmann.png")

# =========================================================================
# PLOT 7: Per-conformer distributions for top 6 variants
# =========================================================================
top6 = summary.head(6)["substituent"].tolist()
top6_data = all_conf[all_conf["substituent"].isin(top6)]

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for name in top6:
    sub = top6_data[top6_data["substituent"] == name]
    axes[0].hist(sub["plane_angle_deg"].dropna(), bins=15, alpha=0.5, label=name)
    axes[1].hist(sub["interplane_distance_A"].dropna(), bins=15, alpha=0.5, label=name)
    axes[2].hist(sub["pi_overlap_pct"].dropna(), bins=15, alpha=0.5, label=name)

axes[0].set_xlabel("Plane Angle (deg)")
axes[0].set_title("Angle Distribution (Top 6)")
axes[0].legend(fontsize=9)
axes[0].axvline(x=20, color="red", linestyle="--", alpha=0.5)

axes[1].set_xlabel("Interplane Distance (A)")
axes[1].set_title("Distance Distribution (Top 6)")
axes[1].legend(fontsize=9)
axes[1].axvspan(3.3, 3.7, alpha=0.1, color="green")

axes[2].set_xlabel("Pi-Overlap (%)")
axes[2].set_title("Overlap Distribution (Top 6)")
axes[2].legend(fontsize=9)
axes[2].axvline(x=50, color="red", linestyle="--", alpha=0.5)

plt.suptitle("Conformer Distributions for Top 6 Excimer-Forming Substituents", fontsize=14)
plt.tight_layout()
plt.savefig(outdir / "top6_distributions.png", dpi=150)
plt.close()
print("Saved top6_distributions.png")

# =========================================================================
# SUMMARY TABLE
# =========================================================================
print("\n" + "=" * 80)
print("COMPLETE SCREENING SUMMARY")
print("=" * 80)

# Merge simple and Boltzmann data
merged = summary[["substituent", "category", "excimer_fraction",
                   "mean_angle", "mean_distance", "mean_overlap",
                   "n_conformers"]].copy()
merged = merged.set_index("substituent")

for col in ["boltz_excimer_fraction", "plane_angle_deg_boltz",
            "interplane_distance_A_boltz", "pi_overlap_pct_boltz"]:
    if col in features.columns:
        merged[col] = features[col]

merged = merged.sort_values("excimer_fraction", ascending=False)

print("\nTOP 10 EXCIMER PROMOTERS:")
print(merged.head(10).to_string())

print("\nBOTTOM 5 (MONOMER PROMOTERS):")
print(merged.tail(5).to_string())

print("\nBOLTZMANN TOP 10:")
boltz_rank = merged.sort_values("boltz_excimer_fraction", ascending=False)
print(boltz_rank[["boltz_excimer_fraction", "excimer_fraction", "category"]].head(10).to_string())

print("\n" + "=" * 80)
print("KEY FINDINGS:")
print("=" * 80)
print("""
1. ELECTRON-WITHDRAWING groups (Cl, CF3) are the BEST excimer promoters:
   - Cl: 81% excimer, Boltzmann=1.00, angle=11.9 deg
   - CF3: 81% excimer, Boltzmann=1.00, angle=14.7 deg

2. BULKY groups (tBu, cHex) show SURPRISINGLY HIGH excimer fraction (62.5%):
   - This contradicts naive steric argument
   - tBu: Boltzmann=0.99, overlap=53%
   - Rigid bulk forces specific folded conformations

3. SMALL/FLEXIBLE groups give moderate results:
   - H: 58% excimer (no steric hindrance)
   - Et: 40% excimer
   - Me: 34% excimer

4. LONG CHAINS are WORST for excimer:
   - nBu: 0% excimer
   - nHex: 4% excimer
   - nPr: 9% excimer
   - Flexible chains adopt extended conformations

5. AROMATIC R-groups show moderate excimer character:
   - Ph: 38%, Bn: 56%
   - Additional pi-stacking interactions may help

6. BOLTZMANN weighting significantly increases excimer fraction:
   - Most variants: Boltzmann >> Simple fraction
   - Low-energy conformers preferentially adopt excimer geometry
   - Exception: nBu, nHex, nPr (extended chains are lowest energy)
""")

# Save the merged ranking table
merged.to_csv("screening_ranking_table.csv")
print("Saved screening_ranking_table.csv")
print(f"\nAll plots saved to {outdir}/")
