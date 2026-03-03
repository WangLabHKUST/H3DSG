from pathlib import Path
import numpy as np
from graph_utils import th_conn_matrix
from constants import celltype_order, g1_samples, g2_samples
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

h5ad_path = Path(r'D:\local_storage\Work\SingleCell\Multiplex-IHC\spatial_analysis\20250929\h5ad')
adata_dict = np.load(h5ad_path / 'CODEX_adata_dict.npy', allow_pickle=True).item()


sample_lst = g1_samples + g2_samples
conn_matrix_dict = {}
for sample_name in sample_lst:
    print(f'Processing sample: {sample_name}')
    adata = adata_dict[sample_name]
    adata = th_conn_matrix(adata, 'r50_adj', 5)
    conn_matrix = adata.obsp['r50_adj_th5_connectivities']
    conn_matrix = conn_matrix.astype(int)
    conn_matrix.eliminate_zeros()
    conn_matrix_dict[sample_name] = conn_matrix



res_matrix_dict = {}; exp_matrix_dict = {}; res_matrix_norm_log2_dict = {}; res_matrix_norm_dict = {}
for sample_name in sample_lst:
    print(f'Processing sample: {sample_name}')
  
    conn_matrix = conn_matrix_dict[sample_name]
    celltype_df = adata_dict[sample_name].obs

    ### test if conn_matrix is symmetric
    assert conn_matrix.shape[0] == celltype_df.shape[0]

    pair_df = pd.DataFrame({'i':conn_matrix.nonzero()[0], 'j':conn_matrix.nonzero()[1]})
    pair_df['i_ct_agg'] = celltype_df.iloc[pair_df['i']]['ct_agg'].values
    pair_df['j_ct_agg'] = celltype_df.iloc[pair_df['j']]['ct_agg'].values
    pair_df.shape
    # pair_df = pair_df[pair_df['i'] > pair_df['j']]
    ct_to_idx = {ct: idx for idx, ct in enumerate(celltype_order)}

    row_idx = pair_df['i_ct_agg'].map(ct_to_idx)
    col_idx = pair_df['j_ct_agg'].map(ct_to_idx)

    res_matrix = pd.crosstab(row_idx, col_idx).reindex(
        index=range(len(celltype_order)),
        columns=range(len(celltype_order)),
        fill_value=0
    )

    res_matrix.columns = celltype_order; res_matrix.index = celltype_order
    np.allclose(res_matrix.values, res_matrix.values.T)

    exp_matrix = pd.DataFrame(index=celltype_order, columns=celltype_order, dtype=float)
    for i, ct1 in enumerate(celltype_order):
        for j, ct2 in enumerate(celltype_order):
            n1 = res_matrix.iloc[i, :].sum()
            n2 = res_matrix.iloc[:, j].sum()
            exp_matrix.iloc[i, j] = (n1 * n2 / pair_df.shape[0])
    # exp_matrix = exp_matrix / 2

    res_matrix_norm = (res_matrix / exp_matrix).astype(float)
    res_matrix_norm_log2 = np.log2(res_matrix_norm)


    res_matrix_dict[sample_name] = res_matrix
    exp_matrix_dict[sample_name] = exp_matrix
    res_matrix_norm_log2_dict[sample_name] = res_matrix_norm_log2
    res_matrix_norm_dict[sample_name] = res_matrix_norm



obs_matrix_all = pd.DataFrame(index=celltype_order, columns=celltype_order, dtype=float, data=0)
exp_matrix_all = pd.DataFrame(index=celltype_order, columns=celltype_order, dtype=float, data=0)
for sample_name in sample_lst:
    this_data = res_matrix_dict[sample_name]#; this_data = this_data.set_index('Unnamed: 0')
    obs_matrix_all += this_data
    this_data = exp_matrix_dict[sample_name]#; this_data = this_data.set_index('Unnamed: 0')
    exp_matrix_all += this_data
res_matrix_all = obs_matrix_all / exp_matrix_all
res_matrix_all_log2 = np.log2(obs_matrix_all / exp_matrix_all + 1)

res_matrix_all_log2_masked = res_matrix_all_log2.copy()
res_matrix_all_log2_masked[res_matrix_all_log2_masked < 0] = np.nan
fig, ax = plt.subplots(figsize=(7,6))
sns.heatmap(res_matrix_all_log2_masked, cmap='Reds', ax=ax)
plt.savefig('all_samples_res_matrix_log2_masked.png',dpi=300, bbox_inches='tight');plt.close()
