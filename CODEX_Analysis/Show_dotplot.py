from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
from constants import celltype_order
h5ad_path = Path(r'D:\local_storage\Work\SingleCell\Multiplex-IHC\spatial_analysis\20250929\h5ad')
adata_dict = np.load(h5ad_path / 'CODEX_adata_dict.npy', allow_pickle=True).item()
adata_full = sc.concat(adata_dict.values(), axis=0)
adata_full.shape



adata_full.obs['ct_agg_ordered'] = pd.Categorical(adata_full.obs['ct_agg'], categories=celltype_order)

fig = sc.pl.dotplot(adata_full, groupby='ct_agg_ordered', var_names=adata_full.var_names, title='Cell type proportions by group',
return_fig=True, standard_scale='var')
fig.savefig(f'01_ct_dotplot_var_red.png', dpi=300, bbox_inches='tight')
plt.close()

