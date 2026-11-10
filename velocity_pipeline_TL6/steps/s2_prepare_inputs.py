"""Step 2: metadata, PCA and UMAP exported from Seurat -> three tables with the common key (orig.ident, file_id, barcode).

Input   config.META_CSV / PCA_CSV / UMAP_CSV (same cell_id order), samples.tsv (seurat_id -> sample)
Output  results/s2_inputs/{meta,harmony,umap}_df_combined.csv
        harmony_df_combined.csv is the table of PCs of the CCA integration; the file name is historical.

The cell_id column is dropped here: step 3 treats every column other than barcode, file_id and orig.ident as coordinates.
"""
import pandas as pd

from _common import parse_args, read_samples, skip_if_done, step_dir, atomic_write

KEY = ['orig.ident', 'file_id', 'barcode']


def main():
    cfg, args = parse_args(__doc__.split('\n')[0])
    out_dir = step_dir(cfg, 2)
    outs = {n: out_dir / f'{n}_df_combined.csv' for n in ('meta', 'harmony', 'umap')}
    if not cfg.UMAP_CSV:
        outs.pop('umap')
    if all(skip_if_done(o, args.force, n) for n, o in outs.items()):
        return

    samples = read_samples(cfg)
    meta = pd.read_csv(cfg.META_CSV)
    pca = pd.read_csv(cfg.PCA_CSV)
    umap = pd.read_csv(cfg.UMAP_CSV) if cfg.UMAP_CSV else None
    for name, df in (('pca', pca), ('umap', umap)):
        if df is not None:
            assert (meta['cell_id'] == df['cell_id']).all(), f'cell_id order of meta and {name} differs'

    file_id = meta['orig.ident'].map(dict(zip(samples['seurat_id'], samples['sample'])))
    missing = meta.loc[file_id.isna(), 'orig.ident'].unique()
    assert len(missing) == 0, f'orig.ident of meta not in the seurat_id column of samples.tsv: {missing}'

    for new_col, src_col in cfg.OBS_ALIASES.items():
        if new_col not in meta.columns and src_col in meta.columns:
            meta[new_col] = meta[src_col]

    def to_schema(df, coord_cols=None):
        out = df.drop(columns=[c for c in ('cell_id', 'file_id') if c in df.columns])
        out.insert(1, 'file_id', file_id.values)
        if coord_cols is not None:
            return out[KEY + coord_cols]
        return out[KEY + [c for c in out.columns if c not in KEY]]

    pc_cols = [c for c in pca.columns if c.startswith('PC_')]
    assert pc_cols, 'no PC_* columns in PCA_CSV'
    tables = {'meta': to_schema(meta), 'harmony': to_schema(pca, pc_cols)}
    if umap is not None:
        tables['umap'] = to_schema(umap, list(cfg.SEURAT_UMAP_COLS))
    for name, df in tables.items():
        atomic_write(df, outs[name])
        print(f'wrote {outs[name]}  shape={df.shape}', flush=True)
    print('\ncells per sample:')
    print(tables['meta'].groupby(['file_id', 'orig.ident']).size(), flush=True)


if __name__ == '__main__':
    main()
