from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import squidpy as sq
from constants import celltype_order, g1_samples, g2_samples, color_dict
import scipy.stats as stats
h5ad_path = Path(r'D:\local_storage\Work\SingleCell\Multiplex-IHC\spatial_analysis\20250929\h5ad')
adata_dict = np.load(h5ad_path / 'CODEX_adata_dict.npy', allow_pickle=True).item()
def zscore_to_pvalue(zscore, alternative='two-sided'):
    if alternative == 'two-sided':
        return 2 * (1 - stats.norm.cdf(abs(zscore)))
    elif alternative == 'greater':
        return 1 - stats.norm.cdf(zscore)
    elif alternative == 'less':
        return stats.norm.cdf(zscore)

th = 40
df_lst = []
sample_lst = list(adata_dict.keys())
def process_sample(sample_name: str) -> None:
    adata = adata_dict[sample_name]
    adata.obs['ct_agg'] = pd.Categorical(adata.obs['ct_agg'], categories=celltype_order)

    if f'r{th}_distances' not in adata.obsp:
        raise RuntimeError(f"Missing 'r{th}_distances' in obsp")

    
    # neighborhood enrichment to get counts/zscores
    sq.gr.nhood_enrichment(adata, cluster_key='ct_agg', connectivity_key=f"r{th}", copy=False)
    res = adata.uns['ct_agg_nhood_enrichment']
    count = res['count']
    zscore = res['zscore']

    # normalization by summed degrees (exclude self-pairs)
    cnvt_matrix_th = adata.obsp[f"r{th}_connectivities"]
    node_degree = cnvt_matrix_th.sum(axis=1).A1
    adata.obs[f"node_degree_r{th}"] = node_degree
    celltype_degree = adata.obs.groupby('ct_agg')[f"node_degree_r{th}"].sum()
    interaction_matrix = count.copy()
    normalization_matrix = celltype_degree.values.reshape(-1, 1) + celltype_degree.values.reshape(1, -1)
    normalization_matrix = normalization_matrix - interaction_matrix
    interaction_matrix = interaction_matrix / normalization_matrix

    # export
    interaction_df = pd.DataFrame(interaction_matrix, index=adata.obs['ct_agg'].cat.categories, columns=adata.obs['ct_agg'].cat.categories)

    zscore_df = pd.DataFrame(zscore, index=adata.obs['ct_agg'].cat.categories, columns=adata.obs['ct_agg'].cat.categories)
    zscore_df = zscore_df.stack().reset_index().fillna(0)
    zscore_df.columns = ['cell_type1', 'cell_type2', 'zscore']

    interaction_df = interaction_df.stack().reset_index().fillna(0)
    interaction_df.columns = ['cell_type1', 'cell_type2', 'count']
    interaction_df = interaction_df.merge(zscore_df, on=['cell_type1', 'cell_type2'], how='left')
    interaction_df = interaction_df.merge(zscore_df, on=['cell_type1', 'cell_type2'], how='left')
    interaction_df['p_two_sided'] = zscore_to_pvalue(interaction_df['zscore'], alternative='two-sided')
    interaction_df['p_one_sided'] = zscore_to_pvalue(interaction_df['zscore'], alternative='ryy-one_sided')
    interaction_df['greater'] = zscore_to_pvalue(interaction_df['zscore'], alternative='greater')
    interaction_df['less'] = zscore_to_pvalue(interaction_df['zscore'], alternative='less')
    df_lst.append(interaction_df)


######### 1. sample-wise summary for clustering coefficient
all_interaction_df = pd.concat(df_lst)
clustering_coef_df = all_interaction_df.loc[all_interaction_df['cell_type1'] == all_interaction_df['cell_type2']]
(clustering_coef_df['p_two_sided'] > 0.05).sum()
clustering_coef_df['cell_type1'].unique()
clustering_coef_df.loc[(clustering_coef_df['p_two_sided'] > 0.05)]
len(sample_lst)

df = clustering_coef_df.copy()
value_col = 'count'
types = df['cell_type1'].unique()
types.sort() 

results = []

for cell_type in types:
    sub = df[df['cell_type1'] == cell_type]
    
    g1_vals = sub[sub['sample'].isin(g1_samples)][value_col].dropna().values
    g2_vals = sub[sub['sample'].isin(g2_samples)][value_col].dropna().values
    
    n_g1 = len(g1_vals)
    n_g2 = len(g2_vals)
    
    stat, p = stats.mannwhitneyu(g1_vals, g2_vals, alternative='two-sided')
    note = ""
    
    med_g1 = np.median(g1_vals) if n_g1 > 0 else np.nan
    med_g2 = np.median(g2_vals) if n_g2 > 0 else np.nan
    
    results.append({
        'cell_type': cell_type,
        'median_G1': round(med_g1, 4) if not np.isnan(med_g1) else np.nan,
        'median_G2': round(med_g2, 4) if not np.isnan(med_g2) else np.nan,
        'p_value': p if not np.isnan(p) else np.nan,
        'enrich': 'Group1' if med_g1 > med_g2 else 'Group2' if med_g1 < med_g2 else 'Equal'
    })

result_df = pd.DataFrame(results)
result_df = result_df.sort_values('p_value')

print(result_df.round(4))