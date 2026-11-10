"""Step 5: UMAP of the final figures, recomputed on X_pca; only a table of coordinates is written, not the h5ad.

Input   obsm['X_pca'] of results/s4_velocity/<stem>_velocity.h5ad (nothing else is read)
Output  results/s5_embedding/<stem>_umap_npc{n}_{metric}_nn{nn}_md{md}.csv
        columns barcode, file_id, x, y; rows in the obs order of the h5ad. Barcodes repeat across samples; file_id + barcode is unique.

umap-learn with a random_state runs on one thread and is reproducible.
"""
from _common import pin_threads

pin_threads()

import warnings

import h5py
import numpy as np
import pandas as pd

try:
    from anndata.io import read_elem
except ImportError:                       # anndata < 0.11
    from anndata.experimental import read_elem

from _common import atomic_write, parse_args, skip_if_done, step_dir


def tag_of(u):
    return f"npc{u['n_pcs']}_{u['metric']}_nn{u['n_neighbors']}_md{u['min_dist']:g}"


def main():
    cfg, args = parse_args(__doc__.split('\n')[0])
    s4, out_dir = step_dir(cfg, 4), step_dir(cfg, 5)
    u = cfg.UMAP
    for stem in cfg.STEMS:
        out = out_dir / f'{stem}_umap_{tag_of(u)}.csv'
        if skip_if_done(out, args.force, stem):
            continue
        with h5py.File(s4 / f'{stem}_velocity.h5ad', 'r') as f:
            obs = read_elem(f['obs'])
            X = np.asarray(read_elem(f['obsm/X_pca']), dtype=float)
        keys = pd.DataFrame({'barcode': obs.index.astype(str), 'file_id': obs['file_id'].astype(str).values})
        assert not (keys['file_id'] + '|' + keys['barcode']).duplicated().any()
        print(f'{stem}: {X.shape[0]} cells, X_pca {X.shape[1]} dims -> UMAP {tag_of(u)}', flush=True)

        warnings.filterwarnings('ignore')
        import umap
        Y = umap.UMAP(n_neighbors=u['n_neighbors'], min_dist=u['min_dist'], metric=u['metric'],
                      random_state=u['random_state']).fit_transform(X[:, :u['n_pcs']])
        df = keys.copy()
        df['x'], df['y'] = Y[:, 0], Y[:, 1]
        atomic_write(df, out)
        print(f'wrote {out}', flush=True)


if __name__ == '__main__':
    main()
