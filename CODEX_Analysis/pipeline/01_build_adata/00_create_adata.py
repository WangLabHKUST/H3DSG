import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
import squidpy as sq
import seaborn as sns
from utils import all_proteins, binary_features

from path_constants import raw_data_path, h5ad_path, csv_path, common_data_path, cell_loc_path, adata_obs_dict_path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from data_loading import to_raw_key
################# input here ##############
sample_name = 'SP08'                   # SP ID; the files below are named by the original key
file_key = to_raw_key(sample_name)
adata_obs_dict = np.load(adata_obs_dict_path, allow_pickle=True).item()
raw_df = pd.read_csv(raw_data_path / f'{file_key}.csv'); raw_df.index = raw_df['Object_ID']
# raw_df.columns
# sns.scatterplot(x='Centroid_X_µm', y='Centroid_Y_µm', data=raw_df, s=0.05)
# plt.show()
raw_data = raw_df[all_proteins].values.astype(float)
binary_data = raw_df[binary_features].values.astype(int)
filtered_data = raw_data * binary_data
# filtered_data.shape
adata_bak = sc.read_h5ad(h5ad_path / f'{file_key}_full_raw.h5ad')
# sns.scatterplot(x='Centroid_X_µm', y='Centroid_Y_µm', data=adata_bak.obs, s=0.05)
# plt.show()

raw_df['r'] = raw_df['Nucleus_Length_µm'] / 2 / np.pi
raw_df['l'] = (raw_df['Nucleus_Max_diameter_µm'] + raw_df['Nucleus_Min_diameter_µm']) / 2
raw_df['r_cell'] = (raw_df['r'] + raw_df['l']) / 2
ct_df = pd.read_csv(csv_path / f'{file_key}.csv', index_col=0)
raw_df['ct_agg'] = ct_df['ct_agg']
adata = ad.AnnData(X=filtered_data, obs=raw_df[['r', 'l', 'r_cell', 'ct_agg']], obsm={'spatial': raw_df[['Centroid_X_µm', 'Centroid_Y_µm']].values})
adata.obsp = adata_bak.obsp
adata.obs.head()

adata.write_h5ad(h5ad_path / f'{file_key}_masked_counts.h5ad')

adata = sc.read_h5ad(h5ad_path / f'{file_key}_masked_counts.h5ad')

# sample_lst_31 of the external utils module holds the original keys, which is also how
# adata_obs_dict and the h5ad files are named, so the loop below stays on those keys
from utils import sample_lst_31 as sample_lst
adata_obs_dict = np.load(h5ad_path / f'adata_obs_dict.npy', allow_pickle=True).item()
cell_data_path = cell_loc_path
for sample_name in sample_lst:
    print(f'Processing {sample_name}...')
    adata = sc.read_h5ad(h5ad_path / f'{sample_name}_masked_counts_no_other_modified.h5ad')
    adata.obs = adata_obs_dict[sample_name]
    adata.obs.head()
    cell_data = pd.DataFrame(index=adata.obs_names, columns=['x', 'y'], data = adata.obsm['spatial'])
    cell_data[['cell_type_2','ct_agg']] = adata.obs[['cell_type_2','ct_agg']]
    cell_data.to_csv(cell_data_path / f'{sample_name}.csv')

