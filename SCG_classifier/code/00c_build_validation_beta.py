"""(optional) Rebuild data/val_data/beta_{5,14,55}_filtered.csv from the raw
per-batch validation beta matrices (~1.2 GB, not shipped), using
data/sample_manifest.csv for the column -> current SP ID mapping.

Usage: python code/00c_build_validation_beta.py [RAW_DIR]
RAW_DIR defaults to <project>/val_data (two levels above this archive) and
must contain beta_5.csv, beta_14.csv, beta_55.csv.

Per batch (same processing as codes/20250613_test.py, minus its 2025-06 IDs):
probe names are cut at the first '_' (EPIC v2 suffixes such as _BC11),
duplicated probes keep the first occurrence, the top-10000 CpGs of
data/beta_score_full.csv that are present are kept in ranking order, only
this study's samples are kept, and columns are renamed to current SP IDs.

If RAW_DIR also holds the legacy beta_{n}_renamed_filtered.csv (columns =
2025-06 SP IDs), its values are compared after translating those IDs through
the manifest's SP_ID_2025_06 -> SP_ID.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from patient_id import check_manifest, load_manifest, load_table_s1

BASE = Path(__file__).resolve().parent.parent
data_dir = BASE / 'data'
out_dir = data_dir / 'val_data'
out_dir.mkdir(exist_ok=True)
raw_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else BASE.parent.parent / 'val_data'

table = load_table_s1(data_dir / 'Table S1.xlsx')
manifest = load_manifest(data_dir / 'sample_manifest.csv')
check_manifest(manifest, table)
raw_columns = pd.read_csv(data_dir / 'id_mapping' / 'raw_val_beta_columns.csv', dtype=str)
beta_score = pd.read_csv(data_dir / 'beta_score_full.csv', index_col=0)
top10k = beta_score.sort_values(by='score', ascending=False).index[:10000]

for n in (5, 14, 55):
    batch = f'beta_{n}'
    out_path = out_dir / f'{batch}_filtered.csv'
    if out_path.exists():
        raise FileExistsError(f'{out_path} already exists, not overwriting')
    raw_path = raw_dir / f'{batch}.csv'
    header = pd.read_csv(raw_path, index_col=0, nrows=0).columns.astype(str).tolist()
    if header != raw_columns.loc[raw_columns['batch'] == batch, 'raw_column'].tolist():
        raise ValueError(f'{raw_path}: header differs from data/id_mapping/raw_val_beta_columns.csv -- not the same raw file')

    m = manifest[manifest['batch'] == batch]
    df = pd.read_csv(raw_path, index_col=0, usecols=[0] + [header.index(c) + 1 for c in m['raw_column']])
    df.columns = df.columns.astype(str)
    df.index = df.index.astype(str).str.split('_').str[0]
    df = df[~df.index.duplicated(keep='first')]
    df = df.loc[[p for p in top10k if p in df.index], m['raw_column'].tolist()]
    df.columns = m['SP_ID'].tolist()
    df.to_csv(out_path)
    print(f'{batch}: {df.shape[0]} CpGs x {df.shape[1]} samples -> {out_path}')

    legacy_path = raw_dir / f'{batch}_renamed_filtered.csv'
    if legacy_path.exists():
        legacy = pd.read_csv(legacy_path, index_col=0)
        legacy.columns = m.set_index('SP_ID_2025_06').loc[legacy.columns, 'SP_ID'].tolist()
        same_cpgs = set(legacy.index) == set(df.index)
        max_diff = (legacy.loc[df.index, df.columns] - df).abs().max().max()
        renamed = {a: b for a, b in zip(m['SP_ID_2025_06'], m['SP_ID']) if a != b}
        print(f'  vs legacy {legacy_path.name}: same CpG set {same_cpgs}, max|diff| = {max_diff:.2e}; '
              f'2025-06 -> current IDs changed: {renamed or "none"}')
