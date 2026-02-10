"""Generate comprehensive MOE vs RDKit comparison plots."""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Load all data
et_moe = pd.read_csv('Et_moe_analyzed.csv')
ipr_moe = pd.read_csv('iPr_moe_analyzed.csv')
et_rdkit = pd.read_csv('Et_rdkit_conformers.csv')
ipr_rdkit = pd.read_csv('iPr_rdkit_conformers.csv')
chex_rdkit = pd.read_csv('cHex_rdkit_conformers.csv')
tbu_rdkit = pd.read_csv('tBu_rdkit_conformers.csv')

def classify(row):
    theta = row['plane_angle_deg']
    d = row['interplane_distance_A']
    overlap = row['pi_overlap_pct']
    if theta < 20 and 3.0 <= d <= 4.0 and overlap > 50:
        return 'strong_excimer'
    elif theta < 60 and d < 4.5 and overlap > 30:
        return 'weak_excimer'
    else:
        return 'monomer'

# Add classification where missing
for df in [et_moe, ipr_moe, chex_rdkit, tbu_rdkit]:
    if 'classification' not in df.columns:
        df['classification'] = df.apply(classify, axis=1)

colors_method = {'MOE': '#2196F3', 'RDKit': '#FF5722'}

# ===== FIGURE 1: Side-by-side geometric distributions (8-panel) =====
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
fig.suptitle('MOE (Amber10:EHT/LowModeMD) vs RDKit (MMFF94s/ETKDGv3) - Geometric Comparison',
             fontsize=14, fontweight='bold')

bins_theta = np.linspace(0, 90, 19)
bins_d = np.linspace(0, 8, 17)
bins_ov = np.linspace(0, 80, 17)

# Row 1: Et
ax = axes[0, 0]
ax.hist(et_moe['plane_angle_deg'], bins=bins_theta, alpha=0.6, color='#2196F3',
        label=f'MOE (n={len(et_moe)})', density=True, edgecolor='white')
ax.hist(et_rdkit['plane_angle_deg'], bins=bins_theta, alpha=0.6, color='#FF5722',
        label=f'RDKit (n={len(et_rdkit)})', density=True, edgecolor='white')
ax.axvline(20, color='green', ls='--', lw=1, alpha=0.7, label='Strong cutoff')
ax.axvline(60, color='red', ls='--', lw=1, alpha=0.7, label='Monomer cutoff')
ax.set_xlabel('Plane Angle (deg)')
ax.set_ylabel('Density')
ax.set_title('Et: Plane Angle Distribution')
ax.legend(fontsize=7)

ax = axes[0, 1]
ax.hist(et_moe['interplane_distance_A'], bins=bins_d, alpha=0.6, color='#2196F3',
        label='MOE', density=True, edgecolor='white')
ax.hist(et_rdkit['interplane_distance_A'], bins=bins_d, alpha=0.6, color='#FF5722',
        label='RDKit', density=True, edgecolor='white')
ax.axvline(3.3, color='green', ls='--', lw=1, alpha=0.7)
ax.axvline(4.5, color='red', ls='--', lw=1, alpha=0.7)
ax.axvspan(3.3, 3.7, alpha=0.1, color='green')
ax.set_xlabel('Inter-plane Distance (A)')
ax.set_ylabel('Density')
ax.set_title('Et: Distance Distribution')
ax.legend(fontsize=8)

ax = axes[0, 2]
ax.hist(et_moe['pi_overlap_pct'], bins=bins_ov, alpha=0.6, color='#2196F3',
        label='MOE', density=True, edgecolor='white')
ax.hist(et_rdkit['pi_overlap_pct'], bins=bins_ov, alpha=0.6, color='#FF5722',
        label='RDKit', density=True, edgecolor='white')
ax.axvline(30, color='orange', ls='--', lw=1, alpha=0.7, label='Weak threshold')
ax.axvline(50, color='green', ls='--', lw=1, alpha=0.7, label='Strong threshold')
ax.set_xlabel('Pi-Overlap (%)')
ax.set_ylabel('Density')
ax.set_title('Et: Pi-Overlap Distribution')
ax.legend(fontsize=8)

ax = axes[0, 3]
sc1 = ax.scatter(et_moe['plane_angle_deg'], et_moe['interplane_distance_A'],
                 c=et_moe['pi_overlap_pct'], cmap='YlOrRd', s=20, alpha=0.5, marker='o',
                 edgecolors='blue', linewidths=0.3, vmin=0, vmax=80)
ax.scatter(et_rdkit['plane_angle_deg'], et_rdkit['interplane_distance_A'],
           c=et_rdkit['pi_overlap_pct'], cmap='YlOrRd', s=30, alpha=0.8, marker='^',
           edgecolors='red', linewidths=0.5, vmin=0, vmax=80)
ax.axhline(4.5, color='red', ls=':', alpha=0.5)
ax.axvline(20, color='green', ls=':', alpha=0.5)
ax.set_xlabel('Plane Angle (deg)')
ax.set_ylabel('Distance (A)')
ax.set_title('Et: Angle vs Distance')
cb = plt.colorbar(sc1, ax=ax, shrink=0.8)
cb.set_label('Overlap %', fontsize=8)
ax.legend([plt.Line2D([0], [0], marker='o', color='blue', ls='', ms=5),
           plt.Line2D([0], [0], marker='^', color='red', ls='', ms=5)],
          ['MOE', 'RDKit'], fontsize=8)

# Row 2: iPr
ax = axes[1, 0]
ax.hist(ipr_moe['plane_angle_deg'], bins=bins_theta, alpha=0.6, color='#2196F3',
        label=f'MOE (n={len(ipr_moe)})', density=True, edgecolor='white')
ax.hist(ipr_rdkit['plane_angle_deg'], bins=bins_theta, alpha=0.6, color='#FF5722',
        label=f'RDKit (n={len(ipr_rdkit)})', density=True, edgecolor='white')
ax.axvline(20, color='green', ls='--', lw=1, alpha=0.7)
ax.axvline(60, color='red', ls='--', lw=1, alpha=0.7)
ax.set_xlabel('Plane Angle (deg)')
ax.set_ylabel('Density')
ax.set_title('iPr: Plane Angle Distribution')
ax.legend(fontsize=8)

ax = axes[1, 1]
ax.hist(ipr_moe['interplane_distance_A'], bins=bins_d, alpha=0.6, color='#2196F3',
        label='MOE', density=True, edgecolor='white')
ax.hist(ipr_rdkit['interplane_distance_A'], bins=bins_d, alpha=0.6, color='#FF5722',
        label='RDKit', density=True, edgecolor='white')
ax.axvline(3.3, color='green', ls='--', lw=1, alpha=0.7)
ax.axvline(4.5, color='red', ls='--', lw=1, alpha=0.7)
ax.axvspan(3.3, 3.7, alpha=0.1, color='green')
ax.set_xlabel('Inter-plane Distance (A)')
ax.set_ylabel('Density')
ax.set_title('iPr: Distance Distribution')
ax.legend(fontsize=8)

ax = axes[1, 2]
ax.hist(ipr_moe['pi_overlap_pct'], bins=bins_ov, alpha=0.6, color='#2196F3',
        label='MOE', density=True, edgecolor='white')
ax.hist(ipr_rdkit['pi_overlap_pct'], bins=bins_ov, alpha=0.6, color='#FF5722',
        label='RDKit', density=True, edgecolor='white')
ax.axvline(30, color='orange', ls='--', lw=1, alpha=0.7)
ax.axvline(50, color='green', ls='--', lw=1, alpha=0.7)
ax.set_xlabel('Pi-Overlap (%)')
ax.set_ylabel('Density')
ax.set_title('iPr: Pi-Overlap Distribution')
ax.legend(fontsize=8)

ax = axes[1, 3]
sc2 = ax.scatter(ipr_moe['plane_angle_deg'], ipr_moe['interplane_distance_A'],
                 c=ipr_moe['pi_overlap_pct'], cmap='YlOrRd', s=20, alpha=0.5, marker='o',
                 edgecolors='blue', linewidths=0.3, vmin=0, vmax=80)
ax.scatter(ipr_rdkit['plane_angle_deg'], ipr_rdkit['interplane_distance_A'],
           c=ipr_rdkit['pi_overlap_pct'], cmap='YlOrRd', s=30, alpha=0.8, marker='^',
           edgecolors='red', linewidths=0.5, vmin=0, vmax=80)
ax.axhline(4.5, color='red', ls=':', alpha=0.5)
ax.axvline(20, color='green', ls=':', alpha=0.5)
ax.set_xlabel('Plane Angle (deg)')
ax.set_ylabel('Distance (A)')
ax.set_title('iPr: Angle vs Distance')
cb = plt.colorbar(sc2, ax=ax, shrink=0.8)
cb.set_label('Overlap %', fontsize=8)
ax.legend([plt.Line2D([0], [0], marker='o', color='blue', ls='', ms=5),
           plt.Line2D([0], [0], marker='^', color='red', ls='', ms=5)],
          ['MOE', 'RDKit'], fontsize=8)

plt.tight_layout()
plt.savefig('plots/moe_vs_rdkit_geometry.png', dpi=200, bbox_inches='tight')
print("Saved: plots/moe_vs_rdkit_geometry.png")
plt.close()

# ===== FIGURE 2: Classification comparison + energy-geometry =====
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
fig.suptitle('MOE vs RDKit: Classification & Energy-Geometry Relationships',
             fontsize=14, fontweight='bold')

# Panel 1: Classification bar chart
ax = axes[0, 0]
datasets = {
    'Et\nMOE': et_moe, 'Et\nRDKit': et_rdkit,
    'iPr\nMOE': ipr_moe, 'iPr\nRDKit': ipr_rdkit,
}
x_pos = np.arange(len(datasets))
width = 0.25
strong_vals, weak_vals, mono_vals = [], [], []
for label, df in datasets.items():
    n = len(df)
    strong_vals.append((df['classification'] == 'strong_excimer').sum() / n * 100)
    weak_vals.append((df['classification'] == 'weak_excimer').sum() / n * 100)
    mono_vals.append((df['classification'] == 'monomer').sum() / n * 100)

ax.bar(x_pos - width, strong_vals, width, label='Strong Excimer', color='#4CAF50')
ax.bar(x_pos, weak_vals, width, label='Weak Excimer', color='#FFC107')
ax.bar(x_pos + width, mono_vals, width, label='Monomer', color='#9E9E9E')
ax.set_ylabel('Fraction (%)')
ax.set_title('Classification Distribution')
ax.set_xticks(x_pos)
ax.set_xticklabels(list(datasets.keys()), fontsize=9)
ax.legend(fontsize=8)
ax.set_ylim(0, 75)

# Panel 2: Energy vs overlap
ax = axes[0, 1]
plot_data = [
    (et_moe, 'Et MOE', '#2196F3', 'o'),
    (et_rdkit, 'Et RDKit', '#FF5722', '^'),
    (ipr_moe, 'iPr MOE', '#1976D2', 's'),
    (ipr_rdkit, 'iPr RDKit', '#E64A19', 'D'),
]
for df, label, color, marker in plot_data:
    e_col = 'rel_energy_kcal_mol' if 'rel_energy_kcal_mol' in df.columns else 'energy_kcal_mol'
    if e_col in df.columns:
        ax.scatter(df[e_col], df['pi_overlap_pct'], alpha=0.4, s=15,
                   label=label, color=color, marker=marker)
ax.set_xlabel('Relative Energy (kcal/mol)')
ax.set_ylabel('Pi-Overlap (%)')
ax.set_title('Energy vs Pi-Overlap')
ax.legend(fontsize=7, loc='upper right')
ax.axhline(30, color='orange', ls='--', alpha=0.5)
ax.axhline(50, color='green', ls='--', alpha=0.5)

# Panel 3: Cumulative overlap distribution
ax = axes[0, 2]
cum_data = [
    (et_moe, 'Et MOE', '#2196F3', '-'),
    (et_rdkit, 'Et RDKit', '#FF5722', '--'),
    (ipr_moe, 'iPr MOE', '#1976D2', '-'),
    (ipr_rdkit, 'iPr RDKit', '#E64A19', '--'),
]
for df, label, color, ls in cum_data:
    sorted_ov = np.sort(df['pi_overlap_pct'])
    cdf = np.arange(1, len(sorted_ov) + 1) / len(sorted_ov) * 100
    ax.plot(sorted_ov, cdf, label=label, color=color, ls=ls, lw=2)
ax.axvline(30, color='orange', ls=':', alpha=0.5, label='Weak threshold')
ax.axvline(50, color='green', ls=':', alpha=0.5, label='Strong threshold')
ax.set_xlabel('Pi-Overlap (%)')
ax.set_ylabel('Cumulative % of Conformers')
ax.set_title('Cumulative Overlap Distribution')
ax.legend(fontsize=7)
ax.set_xlim(0, 80)

# Panel 4: Distance vs overlap
ax = axes[1, 0]
for df, label, color, marker in plot_data:
    ax.scatter(df['interplane_distance_A'], df['pi_overlap_pct'], alpha=0.4, s=15,
               label=label, color=color, marker=marker)
ax.set_xlabel('Inter-plane Distance (A)')
ax.set_ylabel('Pi-Overlap (%)')
ax.set_title('Distance vs Overlap (Key Diagnostic)')
ax.legend(fontsize=7)
ax.axhline(30, color='orange', ls='--', alpha=0.5)
ax.axvline(4.5, color='red', ls='--', alpha=0.5)
ax.axvspan(3.3, 3.7, alpha=0.08, color='green')

# Panel 5: Angle vs overlap
ax = axes[1, 1]
for df, label, color, marker in plot_data:
    ax.scatter(df['plane_angle_deg'], df['pi_overlap_pct'], alpha=0.4, s=15,
               label=label, color=color, marker=marker)
ax.set_xlabel('Plane Angle (deg)')
ax.set_ylabel('Pi-Overlap (%)')
ax.set_title('Angle vs Overlap')
ax.legend(fontsize=7)
ax.axvline(20, color='green', ls='--', alpha=0.5)
ax.axvline(60, color='red', ls='--', alpha=0.5)

# Panel 6: Excimer population (all 6 datasets)
ax = axes[1, 2]
all_labels = ['Et\nMOE', 'Et\nRDKit', 'iPr\nMOE', 'iPr\nRDKit', 'cHex\nRDKit', 'tBu\nRDKit']
all_dfs = [et_moe, et_rdkit, ipr_moe, ipr_rdkit, chex_rdkit, tbu_rdkit]
excimer_pcts = []
for df in all_dfs:
    n = len(df)
    exc = (df['classification'].isin(['strong_excimer', 'weak_excimer'])).sum()
    excimer_pcts.append(exc / n * 100)

bar_colors = ['#2196F3', '#FF5722', '#1976D2', '#E64A19', '#FF9800', '#795548']
ax.bar(range(len(all_labels)), excimer_pcts, color=bar_colors, edgecolor='white')
ax.set_xticks(range(len(all_labels)))
ax.set_xticklabels(all_labels, fontsize=9)
ax.set_ylabel('Excimer Population (%)')
ax.set_title('Excimer Population (All Datasets)')
ax.axhline(50, color='gray', ls='--', alpha=0.5)
for i, v in enumerate(excimer_pcts):
    ax.text(i, v + 1, f'{v:.0f}%', ha='center', fontsize=9, fontweight='bold')
ax.set_ylim(0, 80)

plt.tight_layout()
plt.savefig('plots/moe_vs_rdkit_classification.png', dpi=200, bbox_inches='tight')
print("Saved: plots/moe_vs_rdkit_classification.png")
plt.close()

# ===== FIGURE 3: Excimer quality analysis =====
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
fig.suptitle('Excimer Quality Analysis: MOE vs RDKit', fontsize=14, fontweight='bold')

# Panel 1: Violin plot of overlap for excimer conformers
ax = axes[0]
exc_data = []
exc_labels = []
vcolors = ['#2196F3', '#FF5722', '#1976D2', '#E64A19']
for df, label in [(et_moe, 'Et MOE'), (et_rdkit, 'Et RDKit'),
                   (ipr_moe, 'iPr MOE'), (ipr_rdkit, 'iPr RDKit')]:
    exc = df[df['classification'].isin(['strong_excimer', 'weak_excimer'])]
    if len(exc) > 0:
        exc_data.append(exc['pi_overlap_pct'].values)
        exc_labels.append(f'{label}\n(n={len(exc)})')

positions = range(1, len(exc_data) + 1)
parts = ax.violinplot(exc_data, positions=positions, showmeans=True, showmedians=True)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(vcolors[i])
    pc.set_alpha(0.6)
ax.set_xticks(list(positions))
ax.set_xticklabels(exc_labels, fontsize=9)
ax.set_ylabel('Pi-Overlap (%)')
ax.set_title('Excimer Overlap Quality')
ax.axhline(50, color='green', ls='--', alpha=0.5, label='Strong threshold')
ax.legend(fontsize=8)

# Panel 2: Box plot of distances for excimer conformers
ax = axes[1]
d_data = []
d_labels = []
for df, label in [(et_moe, 'Et MOE'), (et_rdkit, 'Et RDKit'),
                   (ipr_moe, 'iPr MOE'), (ipr_rdkit, 'iPr RDKit')]:
    exc = df[df['classification'].isin(['strong_excimer', 'weak_excimer'])]
    if len(exc) > 0:
        d_data.append(exc['interplane_distance_A'].values)
        d_labels.append(f'{label}\n(n={len(exc)})')

bp = ax.boxplot(d_data, labels=d_labels, patch_artist=True)
for i, box in enumerate(bp['boxes']):
    box.set_facecolor(vcolors[i])
    box.set_alpha(0.6)
ax.axhspan(3.3, 3.7, alpha=0.15, color='green', label='Ideal (3.3-3.7 A)')
ax.set_ylabel('Inter-plane Distance (A)')
ax.set_title('Excimer Stacking Distance')
ax.legend(fontsize=8)

# Panel 3: Strong excimer detail
ax = axes[2]
for df, label, color in [(et_moe, 'Et MOE', '#2196F3'), (ipr_moe, 'iPr MOE', '#1976D2')]:
    strong = df[df['classification'] == 'strong_excimer']
    weak = df[df['classification'] == 'weak_excimer']
    if len(weak) > 0:
        ax.scatter(weak['plane_angle_deg'], weak['pi_overlap_pct'],
                   s=15, alpha=0.2, color=color, zorder=1)
    if len(strong) > 0:
        ax.scatter(strong['plane_angle_deg'], strong['pi_overlap_pct'],
                   s=50, alpha=0.8, color=color, edgecolors='black', zorder=3,
                   label=f'{label} strong (n={len(strong)})')

ax.set_xlabel('Plane Angle (deg)')
ax.set_ylabel('Pi-Overlap (%)')
ax.set_title('Strong vs Weak Excimers (MOE)')
ax.axhline(50, color='green', ls='--', alpha=0.5)
ax.axvline(20, color='green', ls='--', alpha=0.5)
ax.legend(fontsize=8)
ax.set_xlim(-2, 65)

plt.tight_layout()
plt.savefig('plots/excimer_quality_analysis.png', dpi=200, bbox_inches='tight')
print("Saved: plots/excimer_quality_analysis.png")
plt.close()

print("\nAll 3 comparison figures generated successfully!")
