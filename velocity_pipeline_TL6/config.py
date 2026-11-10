"""The only configuration file of the pipeline (plain Python, comments allowed).

For another machine or data set, change the "Raw inputs" section and samples.tsv. The other parameters are the
ones of the final figures. All paths are relative to this folder.
"""
from pathlib import Path

PIPE_DIR = Path(__file__).resolve().parent

# =====================================================================
# 1. Raw inputs (never modified by the pipeline; not in the repository, see README)
# =====================================================================
RAW_ROOT = PIPE_DIR / 'data'

# Sample table: one row per sample = BAM + Cell Ranger barcodes + orig.ident in Seurat
SAMPLE_SHEET = PIPE_DIR / 'samples.tsv'

# Gene annotation for velocyto (GENCODE v44)
GTF = RAW_ROOT / 'ref' / 'gencode.v44.annotation.gtf'

# Three tables exported from Seurat; rows = malignant cells, same cell_id order:
#   META_CSV : cell metadata (side_name / sec_name / State2025 / nn10_prop_* ...)
#   PCA_CSV  : PC_1..PC_30 used for moments and the velocity graph (= X_pca of the h5ad)
#   UMAP_CSV : the original Seurat UMAP (stored as X_umap of the h5ad, not used by the final figures; None to skip)
SEURAT_DIR = RAW_ROOT / 'seurat_tables'
SEURAT_SUFFIX = 'excl_Immune_Stromal_OPC-like'
META_CSV = SEURAT_DIR / f'meta_subset_{SEURAT_SUFFIX}.csv'
PCA_CSV = SEURAT_DIR / f'coords_pca_{SEURAT_SUFFIX}.csv'
UMAP_CSV = SEURAT_DIR / f'coords_umap_subset_{SEURAT_SUFFIX}.csv'
SEURAT_UMAP_COLS = ['UMAPsubset_1', 'UMAPsubset_2']

# Alias columns of the metadata: {new column: column in Seurat}
OBS_ALIASES = {'Filbin_MP_max': 'State2025', 'tumor_clusters_res0.5': 'integrated_snn_res.0.5'}

# =====================================================================
# 2. Output: one subfolder per step, see README
# =====================================================================
OUT_DIR = PIPE_DIR / 'results'

# =====================================================================
# 3. Analysis sets
#    inc-edge = all samples of samples.tsv; exc-edge = without the samples with margin=1 in samples.tsv
#    The final figures use inc-edge only; exc-edge adds another velocity h5ad of about 7 GB
# =====================================================================
STEMS = ['ZDK-inc-edge']
STEM_DEFS = {'ZDK-inc-edge': dict(exclude_margin=False),
             'ZDK-exc-edge': dict(exclude_margin=True)}

# =====================================================================
# 4. Step 1: velocyto run (default logic, no mask)
# =====================================================================
VELOCYTO_PARALLEL = 5          # samples run at the same time
VELOCYTO_SAMTOOLS_THREADS = 8  # samtools threads per sample, if the BAM is not yet sorted by cell
VELOCYTO_SAMTOOLS_MEM_MB = 2000

# =====================================================================
# 5. Step 4: scVelo
# =====================================================================
MOMENTS = dict(n_pcs=30, n_neighbors=30)
VELOCITY_MODE = 'deterministic'
VELOCITY_GRAPH_NEIGHBORS = 10
N_JOBS = 10

# =====================================================================
# 6. Step 5: UMAP of the final figures (recomputed on X_pca, not the Seurat UMAP)
# =====================================================================
UMAP = dict(n_pcs=20, metric='correlation', n_neighbors=100, min_dist=0.6, random_state=0)

# =====================================================================
# 7. Step 6: figures
# =====================================================================
FIG_PREFIX = '20260916'        # date prefix of the file names; fixed, so that a re-run keeps the names
FIG_DPI = 300
FIG_ALPHA = 0.3
NN_K = 10                      # number of neighbours (without the cell itself); denominator of the local dorsal fraction
STREAM_KW = dict(smooth=1.0, min_mass=4, linewidth=1.2, density=2.5, arrow_size=2.0, size=80)
FIG_COLORS = ['side_name', 'sec_name', 'nn10_prop_Dorsal_umap', 'root_cells']
COLOR_SPECS = {
    # categories: Dorsal / Ventral / edge = red / blue / green
    'side_name': dict(palette=['red', 'blue', 'green'], legend_loc='right margin'),
    'sec_name': dict(palette=['red', 'blue', 'green'], legend_loc='right margin'),
    # share of Dorsal cells among the 10 nearest neighbours on the UMAP of step 5; colour scale cut at 0.8
    'nn10_prop_Dorsal_umap': dict(cmap='coolwarm', vmin=0, vmax=0.8),
    # root probability from velocity_pseudotime; depends only on the velocity graph in PCA space
    'root_cells': dict(cmap='gnuplot', vmin=0, vmax=1),
    'end_points': dict(cmap='gnuplot', vmin=0, vmax=1),
}
CBAR_PAD = 0.15                # scVelo draws the colour bar inside the axes, bottom right; blank space on the right takes it
FIGSIZE = (8, 8)
