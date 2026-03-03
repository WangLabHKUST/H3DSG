from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from constants import celltype_order, g1_samples, g2_samples, color_dict
import numpy as np
from matplotlib.patches import Patch
h5ad_path = Path(r'D:\local_storage\Work\SingleCell\Multiplex-IHC\spatial_analysis\20250929\h5ad')

adata_dict = np.load(h5ad_path / 'CODEX_adata_dict.npy', allow_pickle=True).item()

sample_lst = list(adata_dict.keys())
ct_count_df = pd.DataFrame(index=sample_lst, columns=celltype_order)
for sample_name in sample_lst:
    this_stat = adata_dict[sample_name].obs['ct_agg'].value_counts().to_frame().T
    ct_count_df.loc[sample_name] = this_stat[celltype_order].values[0]

ct_count_df['total'] = ct_count_df.sum(axis=1)
## group by group, calculate the total 
ct_count_df.to_csv(f'01_ct_count.csv')
ct_prop_df = ct_count_df.div(ct_count_df['total'], axis=0).drop(columns=['total']).astype(float)
ct_prop_df.to_csv(f'01_ct_prop.csv')

sample_order = g1_samples + g2_samples
n_samples = len(sample_order)
n_group1 = len(g1_samples)

#### stacked bar plot
if True:
    fig, ax = plt.subplots(figsize=(max(14, n_samples * 0.5), 8))
    # Create a blank gap between Group1 (left) and others (right)
    n_left = n_group1
    n_right = n_samples - n_group1
    x_left = np.arange(n_left)
    gap = 0.5  # size of the blank gap in "bar widths"
    x_right = np.arange(n_group1, n_samples) + gap

    bottom_left = np.zeros(n_left)
    bottom_right = np.zeros(n_right)
    for cell_type in celltype_order:
        values = ct_prop_df.loc[sample_order, cell_type].values
        ax.bar(x_left,  values[:n_left],  0.7, bottom=bottom_left,
                color=color_dict[cell_type], alpha=0.7, edgecolor='grey', linewidth=1)
        ax.bar(x_right, values[n_group1:], 0.7, bottom=bottom_right,
                color=color_dict[cell_type], alpha=0.7, edgecolor='grey', linewidth=1)
        bottom_left  += values[:n_left]
        bottom_right += values[n_group1:]

    legend_elements = [Patch(facecolor=color_dict[ct], edgecolor='black', label=ct, alpha=0.7) 
                    for ct in celltype_order]

    ax.set_xlabel('Sample', fontsize=12)
    ax.set_ylabel('Proportion', fontsize=12)
    # ax.set_title('Predicted Cell Type Proportions', fontsize=14, fontweight='bold')
    # Adjust ticks and limits to account for the gap
    x_ticks = np.concatenate([x_left, x_right])
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(sample_order, rotation=45, ha='right')
    ax.set_xlim(-0.5, x_right[-1] + 0.5)
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax.set_ylim(0, 1.05)
    # remove background grid explicitly
    ax.grid(False)
    # hide top and right spines (bounding box edges)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ensure bottom and left spines are black
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    # make bounding box (bottom/left spines) thinner
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    plt.tight_layout()
    plt.savefig('stacked_bar_ct_agg.png', dpi=300, bbox_inches='tight')
    plt.show()
