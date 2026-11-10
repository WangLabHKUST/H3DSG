"""Step 4: scVelo RNA velocity (velocity graph in PCA space, root and end point probabilities).

Input   results/s3_h5ad/<stem>_combined.h5ad
Output  results/s4_velocity/<stem>_velocity.h5ad
        adds to obs: velocity_length / velocity_confidence / velocity_pseudotime / root_cells / end_points
        adds to uns: velocity_graph / velocity_graph_neg (the stream lines of the final figures are projected from these)

X_pca = PC_1..PC_30 exported from Seurat; scv.pp.moments(n_pcs=30) builds the neighbour graph on it, the UMAP plays no part.
"""
import scanpy as sc
import scvelo as scv

from _common import parse_args, skip_if_done, step_dir


def main():
    cfg, args = parse_args(__doc__.split('\n')[0])
    s3, out_dir = step_dir(cfg, 3), step_dir(cfg, 4)
    for stem in cfg.STEMS:
        out = out_dir / f'{stem}_velocity.h5ad'
        if skip_if_done(out, args.force, stem):
            continue
        print(f'Processing {stem}...', flush=True)
        ad = sc.read_h5ad(s3 / f'{stem}_combined.h5ad')
        print(ad.shape, flush=True)
        scv.pp.moments(ad, **cfg.MOMENTS)
        scv.tl.velocity(ad, mode=cfg.VELOCITY_MODE, n_jobs=cfg.N_JOBS)
        scv.tl.velocity_graph(ad, n_neighbors=cfg.VELOCITY_GRAPH_NEIGHBORS, n_jobs=cfg.N_JOBS)
        scv.tl.velocity_confidence(ad)
        scv.tl.velocity_pseudotime(ad)
        tmp = out.with_name(out.name + '.tmp')
        ad.write_h5ad(tmp)
        tmp.replace(out)
        print(f'Wrote {out}', flush=True)


if __name__ == '__main__':
    main()
