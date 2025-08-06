import anndata as ad
from pathlib import Path
import pandas as pd
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
from tqdm import tqdm
import scipy.stats as stats
warnings.filterwarnings("ignore")
from matplotlib.colors import ListedColormap
color_palette = plt.colormaps["tab20c"]

# Set up paths
base_path = Path(r'D:\OneDrive - HKUST Connect\Work\SingleCell\N220708002\codes\py\20250628_imc_data')
sample_lst = list(base_path.glob('csv/*.csv'))
fig_path = base_path / 'fig'
fig_path.mkdir(exist_ok=True)
obj_path = base_path / 'obj'
obj_path.mkdir(exist_ok=True)

st_color_dict = {
    'Cycling': '#00a6ff',
    'OPC_OC': '#000000',
    'AC_MES': '#b62526',
    'Myeloid': '#57b894',
    'Lymphocyte': '#2f609f',
    'Neutrophil': '#97d68f',
    'Vascular': '#b0aad1',
    'Astrocyte': '#fb9a99',
    'Oligo_prog': '#c4eaff',
    # 'Other': '#b9b9b9'
}

g1 = ["X158785", "X165114", "X179186", "X179980"]

# Convert dictionary to ListedColormap
st_color_list = list(st_color_dict.values())
st_colormap = ListedColormap(st_color_list)
# Desired order
cell_type_order = list(st_color_dict.keys())

##### 1. perform computation #####
for sample_path in sample_lst:
    sample_id = sample_path.stem
    print(f"Processing sample: {sample_id}")
    # Read the data
    df = pd.read_csv(sample_path)
    # df = df[df['X'] != 'Astrocyte']
    # Display basic info about the data
    print(f"Data shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    print(f"First few rows:")
    print(df.head())
    # Check if required columns exist
    required_cols = ['X', 'coord_x', 'coord_y']
    # Create AnnData object for squidpy
    # Prepare the data
    # Assuming 'X' is the cell type column, we'll use it as observations
    obs_data = df[['X']].copy()
    obs_data.columns = ['cell_type']
    
    # Create spatial coordinates
    spatial_coords = df[['coord_x', 'coord_y']].values
    cell_type_matrix = pd.get_dummies(df['X'])
    
    # Create AnnData object
    adata = ad.AnnData(
        X=cell_type_matrix,  # Dummy expression matrix
        obs=obs_data,
        obsm={'spatial': spatial_coords}
    )
    # Set 'cell_type' as a categorical field with your desired order
    adata.obs['cell_type'] = pd.Categorical(
        adata.obs['cell_type'],
        categories=cell_type_order,
        ordered=True
    )
    
    sq.gr.spatial_neighbors(adata, coord_type='generic', radius=40)
    ### centrality scores
    ## remove cat not in use
    adata.obs['cell_type'] = adata.obs['cell_type'].cat.remove_unused_categories()
 
    print("Calculating interaction matrix")
    sq.gr.interaction_matrix(adata, cluster_key="cell_type", normalized=False)
    print("Calculating nhood enrichment")
    sq.gr.nhood_enrichment(adata, cluster_key="cell_type")
    adata.write_h5ad(obj_path / f'{sample_id}.h5ad')

##### 2. load back all objects #####
obj_lst = {}
for sample_path in tqdm(sample_lst):
    sample_id = sample_path.stem
    adata = ad.read_h5ad(obj_path / f'{sample_id}.h5ad')
    obj_lst[sample_id] = adata

##### 3. get the degree of each cell (node) #####
for sample_id, adata in obj_lst.items():
    print(sample_id)
    cntv_matrix = adata.obsp['spatial_connectivities']
    node_degree1 = cntv_matrix.sum(axis=1)
    node_degree2 = cntv_matrix.sum(axis=0)
    np.array_equal(node_degree1.T, node_degree2)
    node_degree2_lst = node_degree2[0].tolist()[0]
    node_degree2_lst = [int(x) for x in node_degree2_lst]
    adata.obs['node_degree'] = node_degree2_lst
    celltype_degree = adata.obs.groupby('cell_type')['node_degree'].sum()
    obj_lst[sample_id].obs['celltype_degree'] = celltype_degree
    adata.write_h5ad(obj_path / f'{sample_id}.h5ad')

##### 4. filtering and normalize the interaction matrix #####
interaction_df_norm_lst = []
interaction_matrix_A = pd.DataFrame(index=cell_type_order, columns=cell_type_order)
interaction_matrix_B = pd.DataFrame(index=cell_type_order, columns=cell_type_order)
interaction_matrix_A.fillna(0, inplace=True); interaction_matrix_B.fillna(0, inplace=True)
for sample_id in obj_lst.keys():
    adata = obj_lst[sample_id]
    ## remove unused cell type
    adata.obs['cell_type'] = adata.obs['cell_type'].cat.remove_unused_categories()
    adata.obs['cell_type'].cat.categories
    #### normalize 
    interaction_matrix = adata.uns['cell_type_nhood_enrichment']['count']
    celltype_degree = adata.obs.groupby('cell_type')['node_degree'].sum()
    normalization_matrix = celltype_degree.values.reshape(-1, 1) + celltype_degree.values.reshape(1, -1)
    normalization_matrix = normalization_matrix - interaction_matrix  #### union of the edges with A or B = degree of A + degree of B - the number of edges with A and B
    interaction_matrix_norm = interaction_matrix / normalization_matrix
    interaction_df_norm = pd.DataFrame(interaction_matrix_norm, index=adata.obs['cell_type'].cat.categories, columns=adata.obs['cell_type'].cat.categories)
    if sample_id in g1:
        interaction_matrix_A.loc[interaction_df_norm.index, interaction_df_norm.columns] = interaction_df_norm
    else:
        interaction_matrix_B.loc[interaction_df_norm.index, interaction_df_norm.columns] = interaction_df_norm
    interaction_df_norm = interaction_df_norm.stack().reset_index()
    interaction_df_norm = interaction_df_norm.fillna(0)
    interaction_df_norm.columns = ['cell_type1', 'cell_type2', 'interaction_norm']
    interaction_df_norm.index = interaction_df_norm['cell_type1'] + '_' + interaction_df_norm['cell_type2']
    #### get the z-score of the interaction matrix
    z_score_matrix = adata.uns['cell_type_nhood_enrichment']['zscore']
    z_score_df = pd.DataFrame(z_score_matrix, index=adata.obs['cell_type'].cat.categories, columns=adata.obs['cell_type'].cat.categories)
    z_score_df = z_score_df.stack().reset_index()
    z_score_df = z_score_df.fillna(0)
    z_score_df.columns = ['cell_type1', 'cell_type2', 'z_score']
    z_score_df.index = z_score_df['cell_type1'] + '_' + z_score_df['cell_type2']

    interaction_df_norm['z_score'] = z_score_df.loc[interaction_df_norm.index, 'z_score']
    interaction_df_norm['sample_id'] = sample_id
    interaction_df_norm['group'] = 'PL' if sample_id in g1 else 'TL'

    interaction_df_norm_lst.append(interaction_df_norm)

##### 4.1. visualize the interaction matrix #####
interaction_matrix_A = interaction_matrix_A / len(g1)
interaction_matrix_B = interaction_matrix_B / (len(obj_lst) - len(g1))

### define a function to plot heatmap
def plot_heatmap(matrix, fig_path, colors='Reds', color_bar=True, vmax=None, pad_inches=0):
    import matplotlib.pyplot as plt
    import seaborn as sns

    if vmax is None:
        vmax = matrix.max().max()
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    sns.heatmap(matrix, cmap=colors, vmin=0, vmax=vmax, ax=ax,
                xticklabels=False, yticklabels=False, linecolor='#AAAAAA', linewidths=2,
                cbar=color_bar)
    ax.axis('off')
    plt.tight_layout(pad=pad_inches)
    if pad_inches > 0:
        plt.savefig(fig_path, dpi=300)
    else:
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()

max_val = max(interaction_matrix_A.max().max(), interaction_matrix_B.max().max())
plot_heatmap(interaction_matrix_A, fig_path / '00_interaction_matrix_A_tight.png', 
             color_bar=False, vmax=max_val, pad_inches=2)
plot_heatmap(interaction_matrix_B, fig_path / '00_interaction_matrix_B_tight.png',
              colors='Blues', color_bar=False, vmax=max_val, pad_inches=2)

##### 4.2. combine two images #####
from imageio import imread, imwrite
img1 = imread(fig_path / '00_interaction_matrix_A_tight.png')
img2 = imread(fig_path / '00_interaction_matrix_B_tight.png')
combined_image = np.zeros(img1.shape)
for i in range(img1.shape[0]):
    for j in range(img1.shape[1]):
        if i > j:
            combined_image[i, j, :] = img1[i, j, :]
        else:
            combined_image[i, j, :] = img2[i, j, :]

combined_image = combined_image.astype(np.uint8)
imwrite(fig_path / '00_combined_interaction_matrix.png', combined_image)

##### 5. gather all and drop duplicates #####
interaction_df_all = pd.concat(interaction_df_norm_lst)
interaction_df_all.columns = ['cell_type1', 'cell_type2', 'score', 'z_score', 'sample_id', 'group']
interaction_df_all['cell_type1_sorted'], interaction_df_all['cell_type2_sorted'] = zip(
    *interaction_df_all.apply(lambda row: sorted([row['cell_type1'], row['cell_type2']]), axis=1))
interaction_df_all_nodup = interaction_df_all.drop_duplicates(
    subset=['cell_type1_sorted', 'cell_type2_sorted', 'sample_id', 'group'])
interaction_df_all_nodup.shape
interaction_df_all.shape
interaction_df_all_nodup['pair'] = interaction_df_all_nodup['cell_type1_sorted'] + '_' + interaction_df_all_nodup['cell_type2_sorted']

##### 6. keep significant interactions #####
interaction_df_all_nodup['p_value'] = 2 * stats.norm.sf(np.abs(interaction_df_all_nodup['z_score']))
interaction_df_all_nodup = interaction_df_all_nodup.loc[interaction_df_all_nodup['p_value'] < 0.05,:]

test_result_lst = []
for pair in interaction_df_all_nodup['pair'].unique():
    pair_df = interaction_df_all_nodup[interaction_df_all_nodup['pair'] == pair]
    group1 = pair_df[pair_df['group'] == 'PL']['score']
    group2 = pair_df[pair_df['group'] == 'TL']['score']
    t_onesided = stats.ttest_ind(group1, group2, alternative='less')
    t_onesided2 = stats.ttest_ind(group1, group2, alternative='greater')
    p_value = min(t_onesided.pvalue, t_onesided2.pvalue)
    mean_diff = group1.mean() - group2.mean()
    enrich = 'PL' if mean_diff > 0 else 'TL'
    if p_value > 0.05:
        enrich = 'NS'
    df_this = pd.DataFrame({
        'pair': pair,
        'p_value': p_value,
        't_test_stat': t_onesided.statistic,
        'enrich': enrich,
        'group1_mean': group1.mean(),
        'group2_mean': group2.mean()
    }, index=[0])
    test_result_lst.append(df_this)

test_result_all = pd.concat(test_result_lst)
test_result_all.loc[test_result_all['enrich'] != 'NS']
test_result_all.to_csv(base_path / 'test_result_all.csv', index=False)




