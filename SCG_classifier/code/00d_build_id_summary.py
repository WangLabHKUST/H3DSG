"""Write results/id_summary_all_samples.csv: one row per Table S1 sample with
every ID that sample is known by, so any table from any stage can be joined
back to the current SP ID.

  current SP ID        Table S1 (the ID everything is reported under)
  SP_ID_2025_06        Table S1_0606 ID -- what the methylation matrices were
                       linked to before SP70-SP73 were inserted
  ID_in_network        P01-P36, discovery only ('S' + this == SP ID)
  Multiomics_analysis_ID, PID   the other two discovery-cohort spellings
  Original_ID          hospital number; beta.csv column (discovery) and raw
                       beta_55.csv column
  methy_ID             Sentrix ID; raw beta_5.csv column
  GSM                  GEO accession (GSE270993); raw beta_14.csv column
  batch / beta_column / raw_column   where the methylation values come from

Sources: data/Table S1.xlsx, data/sample_manifest.csv (00b) and
data/ID_relationship_36cases_NF1_TP53.csv. Joins are checked with
patient_id.check_manifest() and by requiring SP ID == 'S' + ID_in_network.
Samples without methylation data have no batch/beta_column/methy_ID/GSM.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from patient_id import FINGERPRINT, check_manifest, load_manifest, load_table_s1

BASE = Path(__file__).resolve().parent.parent
data_dir = BASE / 'data'
out_path = BASE / 'results' / 'id_summary_all_samples.csv'
out_path.parent.mkdir(exist_ok=True)
if out_path.exists():
    raise FileExistsError(f'{out_path} already exists, not overwriting')

table = load_table_s1(data_dir / 'Table S1.xlsx')
manifest = load_manifest(data_dir / 'sample_manifest.csv')
check_manifest(manifest, table)

### mapping only: load the four ID columns and nothing else, so this script cannot
### reach ID_relationship's Group/NF1/TP53 even by accident
rel = pd.read_csv(data_dir / 'ID_relationship_36cases_NF1_TP53.csv',
                  usecols=['ID_in_network', 'Multiomics_analysis_ID', 'Original_ID', 'PID'],
                  dtype={'Original_ID': str})
discovery = manifest[manifest['batch'] == 'discovery'].set_index('Original_ID')['SP_ID']
rel['SP_ID'] = discovery.loc[rel['Original_ID']].values
if not (rel['SP_ID'] == 'S' + rel['ID_in_network']).all():
    raise ValueError('discovery: SP ID is not "S" + ID_in_network:\n'
                     + rel.loc[rel['SP_ID'] != 'S' + rel['ID_in_network'], ['SP_ID', 'ID_in_network']].to_string())

m = manifest.set_index('SP_ID')
r = rel.set_index('SP_ID')
out = table[['Discovery/Validation', 'Resource', 'DNA methyaltion']].copy()
out['SP_ID_2025_06'] = m['SP_ID_2025_06']
for c in ['ID_in_network', 'Multiomics_analysis_ID', 'PID']:
    out[c] = r[c]
out['Original_ID'] = m['Original_ID']
for c in ['methy_ID', 'GSM', 'batch', 'beta_column', 'raw_column']:
    out[c] = m[c]
out[FINGERPRINT] = table[FINGERPRINT]

is_discovery = out['Discovery/Validation'] == 'Discovery'
assert out.loc[is_discovery, 'ID_in_network'].notna().all() and out.loc[~is_discovery, 'ID_in_network'].isna().all()
assert (out['batch'].notna() == (out['DNA methyaltion'].astype(str).str.strip().str.lower() == 'yes')).all()
assert out.loc[is_discovery, 'Original_ID'].equals(r.loc[out.index[is_discovery], 'Original_ID'])
for c in ['Original_ID', 'methy_ID', 'GSM', 'SP_ID_2025_06', 'PID']:
    dup = out[c].dropna()
    assert not dup.duplicated().any(), f'{c} used by more than one sample: {dup[dup.duplicated(keep=False)].to_dict()}'

out.to_csv(out_path)
renumbered = out.index[out['SP_ID_2025_06'].notna() & (out['SP_ID_2025_06'] != out.index)]
print(f'{len(out)} samples: {is_discovery.sum()} discovery, {out["batch"].notna().sum()} with methylation data')
print(f'batches: {out["batch"].value_counts().to_dict()}')
print(f'renumbered since Table S1_0606: {len(renumbered)} '
      f'({out.loc[renumbered, "SP_ID_2025_06"].iloc[0]}->{renumbered[0]} .. '
      f'{out.loc[renumbered, "SP_ID_2025_06"].iloc[-1]}->{renumbered[-1]})' if len(renumbered) else 'no renumbering')
print(f'saved -> {out_path}')
