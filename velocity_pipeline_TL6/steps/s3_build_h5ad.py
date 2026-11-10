"""Step 3: loom + Seurat tables -> AnnData (h5ad), one per sample, then merged per analysis set.

Input   results/s1_loom/<sample>/*.loom, results/s2_inputs/*_df_combined.csv
Output  results/s3_h5ad/<sample>.h5ad                (spliced/unspliced/ambiguous + obs + X_pca [+ X_umap])
        results/s3_h5ad/<stem>_combined.h5ad         (stem = ZDK-inc-edge / ZDK-exc-edge)

Barcodes: a CellID of the loom looks like "<name>:AAACCCAAGCGTCTCGx"; the bare barcode is what follows the colon,
without the final x (and -1). Every barcode of the Seurat metadata must be found in the loom (assert).
"""
import re

import pandas as pd
import scanpy as sc

from _common import (parse_args, read_samples, skip_if_done, step_dir, stem_samples)


def find_loom(s1_dir, sample):
    hits = sorted((s1_dir / sample).glob('*.loom'))
    if len(hits) != 1:
        raise FileNotFoundError(f'[{sample}] exactly one loom needed in {s1_dir / sample}, found: {hits}')
    return hits[0]


def bare_barcode(cell_id):
    return re.sub(r'(-1)?x?$', '', cell_id.split(':', 1)[1])


def build_sample(sample, loom, meta, umap, pca, out):
    ld = sc.read_loom(loom)
    if ':' in ld.obs_names[0]:
        ld.obs_names = [bare_barcode(b) for b in ld.obs_names]
    ld.obs_names_make_unique()
    ld.var_names_make_unique()

    meta_s = meta[meta['file_id'] == sample].set_index('barcode')
    pca_s = pca[pca['file_id'] == sample].drop(columns=['file_id', 'orig.ident']).set_index('barcode').astype(float)
    common = meta_s.index.intersection(ld.obs_names)
    print(f'[{sample}] common barcodes: {len(common)} / meta {len(meta_s)} / loom {ld.n_obs}', flush=True)
    assert len(common) == len(meta_s), f'[{sample}] barcodes of meta not found in the loom'
    ad = ld[common].copy()
    ad.obs = meta_s.loc[ad.obs_names, :]
    ad.obsm['X_pca'] = pca_s.loc[ad.obs_names, :].values
    if umap is not None:
        umap_s = umap[umap['file_id'] == sample].drop(columns=['file_id', 'orig.ident']).set_index('barcode')
        ad.obsm['X_umap'] = umap_s.loc[ad.obs_names, :].values
    ad.write_h5ad(out)


def main():
    cfg, args = parse_args(__doc__.split('\n')[0])
    s1, s2, out_dir = step_dir(cfg, 1), step_dir(cfg, 2), step_dir(cfg, 3)
    samples = read_samples(cfg)
    meta = pd.read_csv(s2 / 'meta_df_combined.csv')
    pca = pd.read_csv(s2 / 'harmony_df_combined.csv')
    umap = pd.read_csv(s2 / 'umap_df_combined.csv') if (s2 / 'umap_df_combined.csv').exists() else None

    in_meta = set(meta['file_id'].unique())
    for sample in samples['sample']:
        if sample not in in_meta:
            print(f'[skip] {sample}: no cells of this sample in meta', flush=True)
            continue
        out = out_dir / f'{sample}.h5ad'
        if skip_if_done(out, args.force, sample):
            continue
        print(f'Processing {sample}...', flush=True)
        build_sample(sample, find_loom(s1, sample), meta, umap, pca, out)

    for stem in cfg.STEMS:
        out = out_dir / f'{stem}_combined.h5ad'
        if skip_if_done(out, args.force, stem):
            continue
        members = [s for s in stem_samples(cfg, stem, samples) if s in in_meta]
        print(stem, members, flush=True)
        ad_lst = []
        for s in members:
            ad = sc.read_h5ad(out_dir / f'{s}.h5ad')
            assert ad.obs['orig.ident'].nunique() == 1
            ad.obs['orig.ident'] = s
            ad_lst.append(ad)
        combined = sc.concat(ad_lst)
        print(combined.obs['orig.ident'].value_counts(), flush=True)
        combined.write_h5ad(out)


if __name__ == '__main__':
    main()
