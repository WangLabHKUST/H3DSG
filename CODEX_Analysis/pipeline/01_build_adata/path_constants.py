from pathlib import Path
import pandas as pd

# Every input of the scripts in this directory is expected under raw_data/ of this module
# (the archive version only; raw_data/ is in .gitignore). Only common_data/ is shipped
# there; the other directories are far too large to copy (QuPath tables 22 GB, h5ad 357 GB,
# cell type tables 2.8 GB) and have to be placed or linked there to run these scripts.
REPO = Path(__file__).resolve().parent.parent.parent
RAW = REPO / 'raw_data'

raw_data_path = RAW / 'qupath_csv'                 # QuPath measurement tables, one per sample
common_data_path = RAW / 'common_data'             # samples.csv, SP_ID_pathology_ID_mapping.xlsx
base_path = RAW / 'detailed_types'                 # interacting-pairs tables by detailed cell type
h5ad_path = RAW / 'h5ad'                           # {sample}_*.h5ad, adata_obs_dict.npy
csv_path = RAW / 'cell_type_csv'                   # cell type tables, one per sample
cell_loc_path = RAW / 'cell_loc'                   # per-sample cell coordinates + cell types
