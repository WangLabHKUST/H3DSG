from pathlib import Path
import pandas as pd
import scanpy as sc
import numpy as np
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from graph_utils import th_conn_matrix
from data_loading import to_sp_id
from path_constants import common_data_path, h5ad_path, csv_path

# samples.csv, the h5ad files and adata_obs_dict are keyed by the original sample keys, so
# sample_lst keeps those keys; to_sp_id() is only used to apply the exclusion by SP ID
sample_info = pd.read_csv(common_data_path / 'samples.csv')
sample_lst = sample_info['sample'].tolist()
sample_lst = [f for f in sample_lst if to_sp_id(f) not in ['SP02', 'SP50']]
# sample_lst = [f for f in sample_lst if to_sp_id(f) not in ['SP08']]

adata_dict = {}
for sample_name in sample_lst:
    print(f'Processing {sample_name}...')
    adata = sc.read_h5ad(h5ad_path / f'{sample_name}_masked_counts_no_other_modified.h5ad')
    adata_dict[sample_name] = adata

ct_dict = {}
for sample_name in sample_lst:
    print(f'Processing {sample_name}...')
    ct_df = pd.read_csv(csv_path / f'{sample_name}.csv', index_col=0)
    ct_dict[sample_name] = ct_df

######### update cell types
ct_col = 'cell_type_2'
to_map = ['T_NK','Myeloid']

celltype_order2 = ['Cycling','OPC_OC','AC_MES', 
                  'Dendritic_cell', 'Myeloid_CD68','Myeloid_CD163','Myeloid_M2_like','Myeloid_other',
                  'Neutrophil',
                  'Tcell_CD4','Tcell_CD8','Treg','NK',
                  'B_cell','Endothelial_cell',
                  'Astrocyte',  'Oligo_pro']

for sample_name in sample_lst:
    print(f'Processing {sample_name}...')
    adata = adata_dict[sample_name]
    # sns.scatterplot(x=adata.obsm['spatial'][:,0], y=adata.obsm['spatial'][:,1], hue=adata.obs['ct_agg'], s=0.005)
    # plt.show()
    ct_df = ct_dict[sample_name]
    len(set(adata.obs_names) - set(ct_df.index))
    ct_df = ct_df.loc[adata.obs_names]
    ct_df['cell_type'].value_counts()
    adata.obs['ct_agg'].value_counts().to_frame().T
    adata.obs['cell_type'] = ct_df['cell_type']
    mapping_mask = adata.obs['ct_agg'].isin(to_map)
    indexs = adata.obs_names[mapping_mask].tolist()
    adata.obs[ct_col] = adata.obs['ct_agg'].copy().astype(str)
    len(set(indexs) - set(ct_df.index))
    len(indexs)
    adata.obs.loc[indexs, ct_col] = ct_df.loc[indexs, 'cell_type']
    ## update all Treg_* to Treg
    adata.obs.loc[adata.obs[ct_col].str.find('Treg') != -1, ct_col] = 'Treg'
    adata.obs[ct_col].value_counts().to_frame()
    adata.obs[ct_col] = adata.obs[ct_col].astype(pd.CategoricalDtype(categories=celltype_order2))
    adata_dict[sample_name] = adata

adata_obs_dict = {}
for sample_name in sample_lst:
    print(f'Processing {sample_name}...')
    adata = adata_dict[sample_name]
    obs_df = adata.obs.copy()[['ct_agg','cell_type',ct_col]]
    adata_obs_dict[sample_name] = obs_df

np.save(h5ad_path / f'adata_obs_dict.npy', adata_obs_dict)
adata_obs_dict = np.load(h5ad_path / f'adata_obs_dict.npy', allow_pickle=True).item()