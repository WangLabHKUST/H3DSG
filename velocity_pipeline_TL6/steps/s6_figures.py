"""Step 6: projects the velocity stream lines onto the UMAP of step 5 and draws one figure per colouring: side_name, sec_name, local dorsal fraction, root_cells.

Input   results/s4_velocity/<stem>_velocity.h5ad (only obs and velocity_graph are read, not the ~7 GB of layers)
        results/s5_embedding/<stem>_umap_<tag>.csv
Output  results/s6_figures/<prefix>_flow_on_umap_<tag>_<color>_<stem>.png      (one per colouring)
        results/s6_figures/<prefix>_nn10_prop_side_<tag>_<stem>.csv            (UMAP coordinates + side composition of the 10 neighbours of each cell)

nn10_prop_Dorsal_umap = share of cells of the Dorsal samples among the 10 nearest neighbours of a cell (without itself)
in the 2D coordinates of this UMAP. The denominator is always 10, so the Dorsal, Ventral and edge shares add up to 1.
It depends on the UMAP parameters and is descriptive only.
The projection uses a light AnnData without genes (only obs + X_umap + velocity_graph/_neg); it gives the same result
per cell as the projection on the full object.
"""
from _common import pin_threads

pin_threads()

import warnings

import anndata as ad
import h5py
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scvelo as scv
from sklearn.neighbors import NearestNeighbors

try:
    from anndata.io import read_elem
except ImportError:
    from anndata.experimental import read_elem

from _common import atomic_write, parse_args, step_dir
from s5_umap import tag_of

warnings.filterwarnings('ignore')


def load(path, colors):
    with h5py.File(path, 'r') as f:
        obs = read_elem(f['obs'])
        graphs = {k: read_elem(f[f'uns/{k}']) for k in ('velocity_graph', 'velocity_graph_neg')}
    keys = obs['file_id'].astype(str).values + '|' + obs.index.astype(str).values
    obs_small = obs[[c for c in colors if c in obs.columns]].copy()
    obs_small.index = keys
    return keys, obs_small, graphs


def knn_side_props(Y, side, k):
    _, idx = NearestNeighbors(n_neighbors=k + 1).fit(Y).kneighbors(Y)
    nb = side[idx[:, 1:k + 1]]                       # column 0 is the cell itself
    return pd.DataFrame({f'nn{k}_prop_{c}_umap': (nb == c).mean(1) for c in pd.unique(side)})


def main():
    cfg, args = parse_args(__doc__.split('\n')[0], extra=lambda p: p.add_argument(
        '--colors', nargs='*', help='overrides config.FIG_COLORS'))
    s4, s5, out_dir = step_dir(cfg, 4), step_dir(cfg, 5), step_dir(cfg, 6)
    colors = args.colors or cfg.FIG_COLORS
    u, k = cfg.UMAP, cfg.NN_K
    tag = tag_of(u)

    for stem in cfg.STEMS:
        csv = s5 / f'{stem}_umap_{tag}.csv'
        if not csv.exists():
            raise SystemExit(f'UMAP coordinates {csv} missing; run step 5 first')
        keys, obs_small, graphs = load(s4 / f'{stem}_velocity.h5ad', list(dict.fromkeys(list(colors) + ['side_name'])))
        df = pd.read_csv(csv)
        assert np.array_equal((df['file_id'].astype(str) + '|' + df['barcode'].astype(str)).values, keys), \
            'row order of the UMAP table differs from the obs of the h5ad'
        Y = df[['x', 'y']].to_numpy(float)

        props = knn_side_props(Y, obs_small['side_name'].astype(str).to_numpy(), k)
        for c in props.columns:
            obs_small[c] = props[c].values
        feat_csv = out_dir / f'{cfg.FIG_PREFIX}_nn{k}_prop_side_{tag}_{stem}.csv'
        atomic_write(pd.concat([df, props], axis=1), feat_csv)
        print(f'saved {feat_csv}', flush=True)

        lite = ad.AnnData(obs=obs_small, obsm={'X_umap': Y},
                          layers={'velocity': np.zeros((len(Y), 0), dtype=np.float32)},
                          uns={**graphs, 'velocity_params': {}})
        scv.tl.velocity_embedding(lite, basis='umap')
        print(f'{stem}  {len(Y)} cells  {tag}', flush=True)

        for color in colors:
            kw = dict(cfg.COLOR_SPECS[color])
            pad = 0 if 'palette' in kw else cfg.CBAR_PAD
            fig, ax = plt.subplots(figsize=(cfg.FIGSIZE[0] * (1 + pad), cfg.FIGSIZE[1]))
            scv.pl.velocity_embedding_stream(lite, basis='umap', color=color, ax=ax, show=False,
                                             alpha=cfg.FIG_ALPHA, **kw, **cfg.STREAM_KW)
            if pad:
                x0, x1 = ax.get_xlim()
                ax.set_xlim(x0, x1 + (x1 - x0) * pad)
                ax.set_title(ax.get_title(), x=0.5 / (1 + pad))      # centre the title over the data area
            kw_txt = ', '.join(f'{a}={b}' for a, b in {**cfg.STREAM_KW, 'alpha': cfg.FIG_ALPHA,
                               **({'vmax': kw['vmax']} if 'vmax' in kw else {})}.items())
            fig.text(0.01, 0.005, f"UMAP on X_pca: n_pcs={u['n_pcs']}, metric={u['metric']}, "
                                  f"n_neighbors={u['n_neighbors']}, min_dist={u['min_dist']:g}\nstream: {kw_txt}",
                     fontsize=7, color='0.35', va='bottom')
            out = out_dir / f'{cfg.FIG_PREFIX}_flow_on_umap_{tag}_{color}_{stem}.png'
            fig.savefig(out, dpi=cfg.FIG_DPI, bbox_inches='tight')
            plt.close(fig)
            print(f'saved {out}', flush=True)


if __name__ == '__main__':
    main()
