import sys
sys.path.append(r'D:\local_storage\Work\SingleCell\Multiplex-IHC\spatial_analysis\codes\all')
from utils import add_large_legend, join_comb

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
import warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
import seaborn as sns
import time
from scipy.sparse import csr_matrix
import networkx as nx
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
from itertools import combinations  
from collections import Counter
from itertools import combinations_with_replacement
from matplotlib import rcParams
from squidpy.gr._utils import (_save_data)
from squidpy._constants._constants import CoordType, Transform
from squidpy._constants._pkg_constants import Key
import scipy
from scipy import stats
import math
import re
tumor_states_color_dict = {'Cycling':'#1C86EE','OPC_OC':'#FF7F00','AC_MES':'#b62526','Myeloid':'#57b894','Neutrophil':'#bdbd45',
'T_NK':'#519e3e','B_cell':'#3c75af','Endothelial_cell':'#b0aad1','Astrocyte':'#f7b6d2','Oligo_pro':'#9edae5'}
def create_color_palette_for_pairs(cols, palette='Set3'):
    color_lst = sns.color_palette(palette, len(cols))
    color_dict = {cols[i]: color_lst[i] for i in range(len(cols))}
    return color_dict


plt.switch_backend('Agg')  # Use non-interactive backend for multiprocessing

def _spatial_scatter_without_squidpy(adata_obj, color_key, connectivity_key, ax, size_key=None, na_color="lightgrey", edge_color="grey", edge_width=0.1, cmap="viridis", color_series=None, scatter_size=10):
    # coordinates
    coords = adata_obj.obsm['spatial']
    # color series (align to adata_obj)
    if isinstance(color_series, pd.Series):
        series = color_series.reindex(adata_obj.obs_names)
    else:
        series = adata_obj.obs[color_key]
    # node colors
    if pd.api.types.is_numeric_dtype(series):
        norm = mcolors.Normalize(vmin=np.nanmin(series.values), vmax=np.nanmax(series.values))
        scalar_map = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        node_colors = scalar_map.to_rgba(series.fillna(np.nan).values)
        nan_mask = series.isna().values
        node_colors[nan_mask] = mcolors.to_rgba(na_color)
    else:
        categories = pd.Categorical(series)
        cats = list(categories.categories)
        palette = sns.color_palette("tab20", n_colors=max(3, len(cats)))
        color_map = {cat: palette[i % len(palette)] for i, cat in enumerate(cats)}
        node_colors = [mcolors.to_rgba(color_map.get(val, na_color)) if pd.notna(val) else mcolors.to_rgba(na_color) for val in categories]
    # node sizes
    if scatter_size is not None:
        sizes = scatter_size
    elif size_key is not None and size_key in adata_obj.obs:
        sizes = adata_obj.obs[size_key].to_numpy()
        sizes = np.where(np.asarray(sizes) > 0, sizes, 10)
    else:
        sizes = 10
    # edges via networkx (based on squidpy _plot_edges)
    if connectivity_key in adata_obj.obsp:
        try:
            from networkx import Graph
            from networkx.drawing import draw_networkx_edges
            g = Graph(adata_obj.obsp[connectivity_key])
            if len(g.edges):
                edge_collection = draw_networkx_edges(
                    g,
                    coords,
                    width=edge_width,
                    edge_color=edge_color,
                    arrows=False,
                    ax=ax,
                )
                if edge_collection is not None:
                    ax.add_collection(edge_collection)
        except Exception:
            # Fallback: draw with LineCollection if networkx approach fails
            adj = adata_obj.obsp[connectivity_key]
            rows, cols = adj.nonzero()
            segs = []
            for i, j in zip(rows, cols):
                if i < j:
                    segs.append([coords[i], coords[j]])
            if len(segs) > 0:
                lc = LineCollection(segs, colors=edge_color, linewidths=edge_width, alpha=0.7)
                ax.add_collection(lc)
    # nodes
    ax.scatter(coords[:, 0], coords[:, 1], c=node_colors, s=sizes, edgecolors='none')
    return ax

def get_middle_points_for_cell(cell_idx, sparse_dist_matrix):
    """获取连接到某个细胞的所有 middle points"""
    col = sparse_dist_matrix[:, cell_idx]
    middle_point_indices = col.nonzero()[0]
    distances = col.data
    return middle_point_indices, distances

def get_cell_in_niche(cell_id, conn_matrix):
    col = conn_matrix[:, cell_id]
    cell_in_niche = col.nonzero()[0].tolist()
    if not cell_id in cell_in_niche:
        cell_in_niche.append(cell_id)
    return cell_in_niche
def th_conn_matrix(adata_used, base_conn_key, th):
    adata_used.obsp[base_conn_key + '_distances'].setdiag(0)
    adata_used.obsp[base_conn_key + '_connectivities'].setdiag(0)
    adata_used.obsp[base_conn_key + '_distances'].eliminate_zeros()
    adata_used.obsp[base_conn_key + '_connectivities'].eliminate_zeros()
    dist_matrix = adata_used.obsp[base_conn_key + '_distances']
    cnvt_matrix = adata_used.obsp[base_conn_key + '_connectivities']
    cnvt_matrix_th = dist_matrix.copy()
    cnvt_matrix_th.data[dist_matrix.data <= th] = 1
    cnvt_matrix_th.data[dist_matrix.data > th] = 0
    cnvt_matrix_th = scipy.sparse.csr_matrix(cnvt_matrix_th)
    cnvt_matrix_th.eliminate_zeros()
    dist_matrix_th = dist_matrix.copy().multiply(cnvt_matrix_th)
    dist_matrix_th.eliminate_zeros()
    adata_used.obs[f'{base_conn_key}_th{th}_degree'] = cnvt_matrix_th.sum(axis=1).A1
    adata_used.obsp[f'{base_conn_key}_th{th}_connectivities'] = cnvt_matrix_th
    adata_used.obsp[f'{base_conn_key}_th{th}_distances'] = scipy.sparse.csr_matrix((len(adata_used.obs), len(adata_used.obs)))# dummy, all zeros
    key_added = f'{base_conn_key}_th{th}'
    neighs_key = Key.uns.spatial_neighs(key_added)
    conns_key = Key.obsp.spatial_conn(key_added)
    dists_key = Key.obsp.spatial_dist(key_added)
    neighbors_dict = {
        "connectivities_key": conns_key,
        "distances_key": dists_key,
        "params": {
            "n_neighbors": 6,
            "coord_type": 'generic',
            "radius": 50,
            "transform": Transform.NONE,
        }
    }
    _save_data(adata_used, attr="obsp", key=key_added, data=cnvt_matrix_th)
    _save_data(adata_used, attr="obsp", key=key_added, data=dist_matrix_th, prefix=False)
    _save_data(adata_used, attr="uns", key=neighs_key, data=neighbors_dict, prefix=False, time=None)
    return adata_used


def sample_cells_by_ct(adata, n_cells=10, conn_key='r50_connectivities', k='ct_agg'):
    adata.obs[k] = adata.obs[k].astype('category')
    ct_lst = adata.obs[k].unique()
    ct_stat = adata.obs[k].value_counts()
    sampled_cells = {}
    existing_cells_whole = []
    existing_cells_whole_str = []
    for ct in ct_lst:
        # adata_sub = adata[adata.obs[k] == ct]
        existing_cells, existing_cells_str = sample_cells_with_mask(adata, n_cells=min(n_cells, ct_stat[ct]), conn_key=conn_key, mask=adata.obs[k] == ct)
        existing_cells_whole += existing_cells
        existing_cells_whole_str += existing_cells_str
        sampled_cells[ct] = [existing_cells, existing_cells_str]
    return sampled_cells, existing_cells_whole, existing_cells_whole_str


def sample_cells_with_mask(adata, n_cells=100, conn_key='r50_connectivities', mask=None, existing_cells_whole=[]):
    if mask is None:
        sampled_bak = np.random.choice(adata.n_obs, n_cells * 10, replace=False)
    else:
        keep_index = adata.obs.reset_index().loc[mask.values].index
        sampled_bak = np.random.choice(keep_index, min(n_cells * 10, len(keep_index)), replace=False)
    existing_cells = [sampled_bak[0]]
    existing_cells_str = [adata.obs_names[sampled_bak[0]]]
    for i, cell_id in enumerate(sampled_bak[1:]):
        conn_this_n_exists = adata.obsp[conn_key][cell_id, existing_cells+existing_cells_whole] 
        conn_cells = adata.obsp[conn_key][cell_id].nnz
        if conn_this_n_exists.sum() == 0 and conn_cells > 2:
            existing_cells.append(cell_id)
            existing_cells_str.append(adata.obs_names[cell_id])
        if len(existing_cells) == n_cells:
            print(f'stopping at {i+1} cells')
            break
    return existing_cells, existing_cells_str

def sample_cells(adata, n_cells=100, conn_key='r50_connectivities'):
    sampled_bak = np.random.choice(adata.n_obs, min(n_cells * 10, adata.n_obs), replace=False)
    existing_cells = [sampled_bak[0]]
    existing_cells_str = [adata.obs_names[sampled_bak[0]]]
    for i, cell_id in enumerate(sampled_bak[1:]):
        conn_this_n_exists = adata.obsp[conn_key][cell_id, existing_cells]
        if conn_this_n_exists.sum() == 0:
            existing_cells.append(cell_id)
            existing_cells_str.append(adata.obs_names[cell_id])
        if len(existing_cells) == n_cells:
            print(f'stopping at {i+1} cells')
            break
    return existing_cells, existing_cells_str


def sample_cells_remove(adata_bak, n_cells=3000, conn_key='r50_connectivities'): ## !!! not achieved, still cannot prevent overlapping
    adata = adata_bak.copy()
    sampled_bak = np.random.choice(adata.n_obs, min(n_cells * 10, adata.n_obs), replace=False)
    sampled_cell_id = sampled_bak.tolist()
    sampled_cell_id_str = adata.obs_names[sampled_cell_id].tolist()
    existing_cells = []
    existing_cells_str = []
    for i, (cell_id_str, cell_id) in enumerate(zip(sampled_cell_id_str, sampled_cell_id)):
        if cell_id_str not in adata.obs_names:
            continue
        new_id = adata.obs_names.tolist().index(cell_id_str)
        existing_cells_str.append(cell_id_str)
        existing_cells.append(cell_id)
        connected_cells = adata.obsp[conn_key][new_id].nonzero()[1].tolist() + [new_id]
        # fig, ax = plt.subplots(figsize=(10, 10))
        # plot_niche(cell_id_str, adata[connected_cells+[new_id]], ax=ax)
        # for cell_id in connected_cells:
        #     ax.add_patch(plt.Circle(adata.obsm['spatial'][cell_id], 10, color='red', fill=True, linewidth=1, alpha=0.5))
        # plt.savefig(f'{cell_id_str}_connected_cells.png', dpi=300, bbox_inches='tight')
        # plt.close()
        # break
        ## remove connected_cells
        remove_cells = adata.obs_names[connected_cells].tolist()
        adata = adata[~adata.obs_names.isin(remove_cells)]
        if len(existing_cells_str) % 100 == 0:
            print(f'current adata size: {adata.n_obs}')
        # print(f'removing {len(remove_cells)} cells')
        if len(existing_cells_str) == n_cells:
            print(f'stopping at {i+1} cells')
            break
    return existing_cells, existing_cells_str


def get_adata_sub(cell_id_str, adata, conn_key):
    cell_id = adata.obs_names.tolist().index(cell_id_str)
    cell_in_niche = get_cell_in_niche(cell_id, adata.obsp[conn_key])
    cell_in_niche.remove(cell_id)
    cell_in_niche = [cell_id] + cell_in_niche
    adata_sub = adata[cell_in_niche].copy()
    return adata_sub

def plot_niche(cell_id_str, adata_sub, figsize=(10, 10), ax=None, base_conn_key='r50',th_conn_key='r50_adj_th5', tumor_states_color_dict=tumor_states_color_dict, k='ct_agg', edge_width=1, scatter_size=10, niche_radius=50,out_bound_color='black',out_bound_width=1):
    cell_id = adata_sub.obs_names.tolist().index(cell_id_str)
    radius = int(base_conn_key.strip('r').split('_')[0])
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
    _spatial_scatter_without_squidpy(adata_sub, color_key=k, connectivity_key=th_conn_key + '_connectivities', edge_width=edge_width,ax=ax, size_key='size_show', na_color="lightgrey", scatter_size=scatter_size)
    
    # filled circles for all cells in the bounding box
    for idx in range(adata_sub.n_obs):
        radius = adata_sub.obs["r_cell"][idx] if "r_cell" in adata_sub.obs else 10
        ax.add_patch(
            plt.Circle(
                adata_sub.obsm['spatial'][idx],
                radius,
                facecolor=tumor_states_color_dict[adata_sub.obs[k].values[idx]],
                edgecolor=tumor_states_color_dict[adata_sub.obs[k].values[idx]],
                fill=True,
                linewidth=1,
                alpha=0.5,
            )
        )   
    ## add a big circle for the cell
    ax.add_patch(
        plt.Circle(
            adata_sub.obsm['spatial'][cell_id],
            niche_radius,
            color=out_bound_color,
            fill=False,
            linewidth=out_bound_width,
        )
    )
    ## add large legend
    add_large_legend(ax=ax, hue_order=tumor_states_color_dict.keys(), palette=tumor_states_color_dict, data=adata_sub.obs, hue='ct_agg', title='Cell Type', marker_size=150, loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=1, fontsize=12, title_fontsize=14, frameon=True, edgecolor='black', linewidth=0.8)
    ax.invert_yaxis()
    ax.axis('equal')
    ### add x y axis ticks
    plt.tight_layout()
    ### debug
    # plt.savefig(out_path / f'{sample_name}_cell_id_{cell_id_str}_niche.png', dpi=300, bbox_inches='tight'); plt.close()
    return ax
def get_pairs_in_niche(adata_used, conn_key, k='ct_agg'):
    conn_matrix = adata_used.obsp[conn_key]
    pairs_df = pd.DataFrame({'i':conn_matrix.nonzero()[0], 'j':conn_matrix.nonzero()[1]})
    pairs_df['i_ct_agg'] = adata_used.obs.iloc[pairs_df['i']]['ct_agg'].values
    pairs_df['j_ct_agg'] = adata_used.obs.iloc[pairs_df['j']]['ct_agg'].values
    pairs_df.shape
    pairs_df = join_comb(pairs_df, celltype_order=adata_used.obs[k].cat.categories, class_name='ct_agg')
    return pairs_df

def show_pairs_count(pairs_count_df, celltype_count, ax=None, celltype_color_dict=tumor_states_color_dict, pair_palette='Set3'):
    total_pairs = pairs_count_df.values.sum(); total_celltypes = celltype_count.values.sum()
    all_count = pd.concat([ celltype_count,pairs_count_df], axis=1).unstack().to_frame().reset_index().drop(columns=['level_1'])
    all_count.columns = ['cell_type', 'count']
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    else:
        fig = ax.get_figure()
    ### customized palette
    pair_palette = create_color_palette_for_pairs(pairs_count_df.columns, palette=pair_palette)
    print(pair_palette)
    this_palette = {k:tumor_states_color_dict[k] for k in celltype_count.columns}
    this_palette.update(pair_palette)

    sns.barplot(all_count, x='cell_type', y='count', estimator="sum", errorbar=None, ax=ax, palette=this_palette)
    for i in range(len(ax.containers)):
        ax.bar_label(ax.containers[i], fontsize=10)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_title(f'Total Pairs: {total_pairs}, Total Celltypes: {total_celltypes}')
    plt.tight_layout()
    return ax


def plot_niche_and_pairs(cell_id_str,cell_id_overall, adata_sub, pairs_count_df, celltype_count, figsize=(18, 8), out_path=None, sample_name=None, ax=None, k='ct_agg', base_conn_key='r50', th=5):
    if ax is None:
        fig, axes = plt.subplots(1, 2, figsize=figsize, gridspec_kw={'width_ratios': [1.7, 1]})  
    else:
        fig = ax.get_figure()
        axes = [ax]
    plt.sca(axes[0])  # 设置当前 axes 为左图
    # cell_id = adata_sub.obs_names.tolist().index(cell_id_str)
    ax1 = plot_niche(cell_id_str, adata_sub, figsize=(figsize[0] * (1.7/(1.7+1)), figsize[1]), ax=axes[0])  # 它会在 axes[0] 上画图
    axes[0].set_title(f'Cell ID: {cell_id_overall} {cell_id_str}')  # 可选：加标题


    # 同理处理 show_pairs_count
    plt.sca(axes[1])  # 切换到右图
    ax2 = show_pairs_count(pairs_count_df, celltype_count, ax=axes[1])  # 假设它返回 ax
    axes[1].set_title('Interaction Pair Counts')

    # 整体调整布局和标题
    # fig.suptitle(f'Sample: {sample_name} | Cell ID: {cell_id}', fontsize=16)
    fig.suptitle(f'Total Pairs: {pairs_count_df.values.sum()}, Total Celltypes: {celltype_count.values.sum()}')
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # 留出 suptitle 的空间

    # 保存合并后的图
    if out_path is not None:
        save_path = out_path / f'{sample_name}_cell_id_{cell_id_overall}_niche_and_pairs.png'
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)  # 记得关闭
        print(f"Saved combined figure to: {save_path}")
    return ax
