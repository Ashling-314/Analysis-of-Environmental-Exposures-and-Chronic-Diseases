import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

def sf(v):
    try: return float(v)
    except: return None

exposure_en = {
    'Smoking_initiation': 'Smoking',
    'Alcohol_intake_freq': 'Alcohol',
    'Coffee_intake': 'Coffee',
    'NO2_air_pollution': 'NO2',
    'Fresh_fruit_intake': 'Fresh Fruit'
}
outcome_en = {
    'COPD': 'COPD',
    'LungCancer': 'Lung Cancer',
    'CAD': 'CAD',
    'HeartFailure': 'Heart Failure',
    'T2D': 'T2D',
    'IschStroke': 'Isch. Stroke',
    'Alzheimer': 'Alzheimer',
    'MDD': 'MDD'
}
exp_colors = {
    'Smoking_initiation': '#d62728',
    'Alcohol_intake_freq': '#ff7f0e',
    'Coffee_intake': '#2ca02c',
    'NO2_air_pollution': '#1f77b4',
    'Fresh_fruit_intake': '#9467bd'
}

# ============ Figure 1: Forest Plot ============
with open('results/tables/03_batch_mr_results.csv') as f:
    p2 = list(csv.DictReader(f))
sig = [r for r in p2 if r['method'] == 'IVW' and sf(r['pval']) and float(r['pval']) < 0.05]
sig.sort(key=lambda x: float(x['pval']))

fig, ax = plt.subplots(figsize=(10, 7))
fig.patch.set_facecolor('white')
labels, or_vals, or_lo, or_hi, pvals, colors = [], [], [], [], [], []
for r in sig:
    labels.append(f"{exposure_en.get(r['exposure'], r['exposure'])} -> {outcome_en.get(r['outcome'], r['outcome'])}")
    or_vals.append(float(r['or']))
    or_lo.append(float(r['or_lci']))
    or_hi.append(float(r['or_uci']))
    pvals.append(float(r['pval']))
    colors.append(exp_colors.get(r['exposure'], '#333333'))

y_pos = range(len(labels))
for i, y in enumerate(y_pos):
    ax.errorbar(or_vals[i], y, xerr=[[or_vals[i]-or_lo[i]], [or_hi[i]-or_vals[i]]],
                fmt='o', color=colors[i], markersize=8, capsize=4, linewidth=1.5,
                markeredgecolor='black', markeredgewidth=0.5)
ax.axvline(x=1, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel('Odds Ratio (95% CI)', fontsize=12)
ax.set_title('Figure 1. Mendelian Randomization Causal Estimates (IVW, P < 0.05)', fontsize=13, fontweight='bold')
ax.set_xlim(0.2, 4.2)
ax.invert_yaxis()
for i, y in enumerate(y_pos):
    p = pvals[i]
    pstr = f"P = {p:.2e}" if p < 0.001 else f"P = {p:.3f}"
    ax.text(4.25, y, pstr, fontsize=9, va='center', color='#555555')
legend_patches = [mpatches.Patch(color=c, label=l) for l, c in [
    ('Smoking', '#d62728'), ('Alcohol', '#ff7f0e'),
    ('Coffee', '#2ca02c'), ('NO2', '#1f77b4'), ('Fresh Fruit', '#9467bd')
]]
ax.legend(handles=legend_patches, loc='lower right', fontsize=9, framealpha=0.9)
plt.tight_layout()
plt.savefig('results/figures/Fig1_forest_plot.png', dpi=300, bbox_inches='tight')
plt.close()
print("Figure 1: Forest plot (English)")

# ============ Figure 2: Mediation MR ============
with open('results/tables/06_mediation_results.csv') as f:
    p5 = list(csv.DictReader(f))
sig_med = [r for r in p5 if r['sig_mediation'] == 'TRUE']
sig_med.sort(key=lambda x: abs(float(x['prop_mediated'])), reverse=True)

outcome_en_med = {'COPD': 'COPD', 'LungCancer': 'Lung Cancer', 'CAD': 'CAD', 'T2D': 'T2D', 'MDD': 'MDD'}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
fig.patch.set_facecolor('white')
gene_labels = [f"{r['gene']} -> {outcome_en_med.get(r['outcome'], r['outcome'])}" for r in sig_med]
med_effects = [float(r['mediation_effect']) for r in sig_med]
props = [float(r['prop_mediated']) for r in sig_med]
bar_colors = ['#d62728' if m > 0 else '#1f77b4' for m in med_effects]

y = range(len(gene_labels))
ax1.barh(y, med_effects, color=bar_colors, edgecolor='black', linewidth=0.5)
ax1.set_yticks(y)
ax1.set_yticklabels(gene_labels, fontsize=9)
ax1.set_xlabel('Mediation Effect (beta1 x beta2)', fontsize=10)
ax1.set_title('A. Mediation Effect', fontsize=11, fontweight='bold')
ax1.axvline(x=0, color='gray', linewidth=0.8)
ax1.invert_yaxis()

bar_colors2 = ['#2ca02c' if abs(p) <= 1 else '#ff7f0e' if abs(p) <= 2 else '#d62728' for p in props]
ax2.barh(y, props, color=bar_colors2, edgecolor='black', linewidth=0.5)
ax2.set_yticks(y)
ax2.set_yticklabels(gene_labels, fontsize=9)
ax2.set_xlabel('Proportion Mediated', fontsize=10)
ax2.set_title('B. Proportion Mediated', fontsize=11, fontweight='bold')
ax2.axvline(x=0, color='gray', linewidth=0.8)
ax2.axvline(x=1, color='red', linestyle='--', linewidth=0.8, alpha=0.5)
ax2.axvline(x=-1, color='red', linestyle='--', linewidth=0.8, alpha=0.5)
ax2.invert_yaxis()

plt.suptitle('Figure 2. Significant Mediation Pathways (Smoking -> Gene Expression -> Disease)',
             fontsize=12, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('results/figures/Fig2_mediation.png', dpi=300, bbox_inches='tight')
plt.close()
print("Figure 2: Mediation (English)")

# ============ Figure 3: GTEx ============
with open('results/tables/10_gtex_multitissue_validation.csv') as f:
    gtex = list(csv.DictReader(f))
gtex_data = {}
for r in gtex:
    if r['gene'] not in gtex_data:
        gtex_data[r['gene']] = {}
    gtex_data[r['gene']][r['tissue']] = r

genes_order = list(gtex_data.keys())
tissues = ['Lung', 'Whole_Blood']
tissue_labels = ['Lung', 'Whole Blood']

matrix = []
for g in genes_order:
    row = []
    for t in tissues:
        r = gtex_data[g].get(t, {})
        if r.get('found') == 'True':
            if r.get('is_egenie') == 'True':
                qval = float(r.get('qval', 1))
                row.append(-np.log10(max(qval, 1e-300)))
            else:
                row.append(0)
        else:
            row.append(np.nan)
    matrix.append(row)
matrix = np.array(matrix)

fig, ax = plt.subplots(figsize=(6, 5))
fig.patch.set_facecolor('white')
display_matrix = np.full_like(matrix, -1, dtype=float)
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        if np.isnan(matrix[i, j]): display_matrix[i, j] = -1
        elif matrix[i, j] == 0: display_matrix[i, j] = 0
        else: display_matrix[i, j] = min(matrix[i, j], 10)

from matplotlib.colors import ListedColormap, BoundaryNorm
cmap = ListedColormap(['#e0e0e0', '#ffffff', '#a8d8a8', '#2ca02c'])
bounds = [-1.5, -0.5, 0.5, 2, 10.5]
norm = BoundaryNorm(bounds, cmap.N)

im = ax.imshow(display_matrix, cmap=cmap, norm=norm, aspect='auto')
ax.set_xticks(range(len(tissues)))
ax.set_xticklabels(tissue_labels, fontsize=10)
ax.set_yticks(range(len(genes_order)))
ax.set_yticklabels(genes_order, fontsize=9)
ax.set_title('Figure 3. GTEx v8 Multi-Tissue eQTL Validation\n(-log10 q-value, darker = more significant)',
             fontsize=11, fontweight='bold')
for i in range(len(genes_order)):
    for j in range(len(tissues)):
        if np.isnan(matrix[i, j]): txt = 'NT'
        elif matrix[i, j] == 0: txt = 'ns'
        else: txt = f"{matrix[i,j]:.1f}"
        ax.text(j, i, txt, ha='center', va='center', fontsize=8,
                color='black', fontweight='bold' if matrix[i,j] > 2 else 'normal')

legend_elements = [
    mpatches.Patch(facecolor='#2ca02c', edgecolor='black', label='eGene (q < 0.05)'),
    mpatches.Patch(facecolor='#a8d8a8', edgecolor='black', label='eGene (q < 0.01)'),
    mpatches.Patch(facecolor='#ffffff', edgecolor='black', label='Not eGene'),
    mpatches.Patch(facecolor='#e0e0e0', edgecolor='black', label='Not tested')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=7, framealpha=0.9)
plt.tight_layout()
plt.savefig('results/figures/Fig3_gtex_validation.png', dpi=300, bbox_inches='tight')
plt.close()
print("Figure 3: GTEx (English)")

# ============ Figure 4: Colocalization ============
with open('results/tables/07_coloc_results.csv') as f:
    p7_150 = list(csv.DictReader(f))
with open('results/tables/09_coloc_narrow_50kb.csv') as f:
    p7_50 = list(csv.DictReader(f))

pairs = [f"{r['exposure']} -> {r['outcome']}" for r in p7_150]
h4_150 = [float(r['PP_H4']) for r in p7_150]
p50_map = {f"{r['exposure']} -> {r['outcome']}": float(r['PP_H4']) for r in p7_50}
h4_50 = [p50_map.get(p, 0) for p in pairs]

fig, ax = plt.subplots(figsize=(8, 5))
fig.patch.set_facecolor('white')
x = np.arange(len(pairs))
width = 0.35
bars1 = ax.bar(x - width/2, h4_150, width, label='+-150 kb', color='#58a6ff', edgecolor='black', linewidth=0.5)
bars2 = ax.bar(x + width/2, h4_50, width, label='+-50 kb', color='#f85149', edgecolor='black', linewidth=0.5)
ax.axhline(y=0.5, color='red', linestyle='--', linewidth=1, alpha=0.7, label='PP.H4 = 0.5 threshold')
ax.set_xticks(x)
ax.set_xticklabels(pairs, fontsize=9, rotation=30, ha='right')
ax.set_ylabel('PP.H4 (Posterior Probability)', fontsize=10)
ax.set_title('Figure 4. Colocalization Sensitivity Analysis\n(+-150 kb vs +-50 kb Windows)',
             fontsize=11, fontweight='bold')
ax.legend(fontsize=8)
ax.set_ylim(0, 0.6)
for bar in bars1:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
            f'{bar.get_height():.3f}', ha='center', va='bottom', fontsize=7)
for bar in bars2:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
            f'{bar.get_height():.3f}', ha='center', va='bottom', fontsize=7)
plt.tight_layout()
plt.savefig('results/figures/Fig4_coloc_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("Figure 4: Coloc (English)")

print("\nAll 4 figures regenerated with English labels")
