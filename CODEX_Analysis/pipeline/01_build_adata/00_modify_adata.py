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

from path_constants import raw_data_path, h5ad_path, csv_path, common_data_path, adata_obs_dict_path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from data_loading import to_raw_key
################# input here ##############
sample_name = 'SP08'                   # SP ID; the files below are named by the original key
file_key = to_raw_key(sample_name)
adata = sc.read_h5ad(h5ad_path / f'{file_key}_masked_counts.h5ad')
adata = adata[~adata.obs['ct_agg'].isin(['Other', 'Tumor_other'])]
adata.obs['ct_agg'].value_counts().to_frame().T
print(f'writing {sample_name}_masked_counts_no_other...')
adata.write_h5ad(h5ad_path / f'{file_key}_masked_counts_no_other.h5ad')

from graph_interaction_plotting_utils import save_graph_obsp, clean_uns_dict
print(f'adding spatial neighbors for {sample_name}...')
sq.gr.spatial_neighbors(adata, coord_type='generic',delaunay=True)
adata = save_graph_obsp(adata, key='spatial', new_key = 'de')
for th in [40,50]:
       sq.gr.spatial_neighbors(adata, coord_type='generic',radius=th)
       adata = save_graph_obsp(adata, key='spatial', new_key = f'r{th}')
sq.gr.spatial_neighbors(adata, coord_type='generic',radius=27.5)
adata = save_graph_obsp(adata, key='spatial', new_key = f'r27_5')
adata.uns = clean_uns_dict(adata.uns)
print(adata.uns.keys())
print(adata.obsp.keys())
for key in adata.uns.keys():
       k2 = key.replace('_neighbors', '')
       print(k2,'num edges:',adata.obsp[f'{k2}_connectivities'].nnz, 'max distance:', adata.obsp[f'{k2}_distances'].max())
print(f'writing {sample_name}...')
adata.write_h5ad(h5ad_path / f'{file_key}_masked_counts_no_other_modified.h5ad')


