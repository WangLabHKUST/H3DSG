"""Build data/sample_manifest.csv: one row per Table S1 sample, all 86, so that
anything keyed by SP ID can be joined against it directly.

The 73 methylation samples (36 discovery + 37 validation) link the column name
in each beta matrix to the current Table S1 SP ID and carry the patient
fingerprint (Age, Sex, NF1, TP53) that 01/02/03 check against Table S1 on every
run. The other 13 have no methylation data: only SP_ID and cohort are filled,
every other column is empty, and consumers decide "has methylation data" by
'batch' being set rather than by a row being present.

ID chain (only IDs that never change are used to reach the 2025-06 SP ID):
  discovery : beta.csv column = Original_ID -> ID_relationship P-number -> SP ID (2025-06)
  beta_55   : raw column = Original_ID      -> Methy_Cell_ID_map_to_yy 'Cell ID' (2025-06 SP ID)
  beta_5    : raw column = Sentrix position -> same map
  beta_14   : raw column = GEO GSM          -> GSE270993 sample_info (Sentrix) -> same map
Fingerprints are taken from Table S1_0606.xlsx, the table those 2025-06 IDs
refer to.

2025-06 SP ID -> current SP ID: Table S1 must equal Table S1_0606 with rows
inserted only. Each old row is matched, in order, to a Table S1 row with the
same Age/Sex/NF1/TP53/Resource/Discovery-Validation. The match is accepted
only if the left-most and right-most order-preserving matches coincide, i.e.
it is unique -- so an inserted row can never be taken for an old one, even
where fingerprints repeat. Inserted rows are reported and must not carry
methylation data.

Sources (data/id_mapping/): Methy_Cell_ID_map_to_yy.xlsx, Table S1_0606.xlsx,
GSE270993_sample_info.csv, raw_val_beta_columns.csv (column headers of the
raw val_data/beta_{5,14,55}.csv, which are too large to ship).
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from patient_id import CORRECTIONS, FINGERPRINT, check_manifest, fingerprint, load_manifest, load_table_s1

BASE = Path(__file__).resolve().parent.parent
data_dir = BASE / 'data'
src_dir = data_dir / 'id_mapping'
out_path = data_dir / 'sample_manifest.csv'
if out_path.exists():
    raise FileExistsError(f'{out_path} already exists, not overwriting')

table = load_table_s1(data_dir / 'Table S1.xlsx')
old = pd.read_excel(src_dir / 'Table S1_0606.xlsx', dtype={'Sample ID': str}).set_index('Sample ID').sort_index()
assert old.index.is_unique


### ==================== 1. 2025-06 SP ID -> current Table S1 SP ID ====================
def align_keys(df, fields=FINGERPRINT):
    k = fingerprint(df)[list(fields)]
    for c in ('Resource', 'Discovery/Validation'):
        k[c] = df[c].astype(str).str.strip()
    return k.agg('|'.join, axis=1).tolist()


def leftmost_match(a, b):
    """Positions in b of the left-most order-preserving match of every element of a, or None."""
    pos, j = [], 0
    for key in a:
        while j < len(b) and b[j] != key:
            j += 1
        if j == len(b):
            return None
        pos.append(j)
        j += 1
    return pos


def try_align(fields):
    """Unique order-preserving match of old onto table on `fields`, or None."""
    old_keys, new_keys = align_keys(old, fields), align_keys(table, fields)
    left = leftmost_match(old_keys, new_keys)
    right = leftmost_match(old_keys[::-1], new_keys[::-1])
    if left is None or right is None:
        return None
    right = [len(new_keys) - 1 - p for p in reversed(right)]
    return left if left == right else None


### A clinical value corrected after 2025-06 (patient_id.CORRECTIONS) is not a
### changed patient, but it does change the fingerprint, so the full-field match
### fails. Fall back to a key without the corrected fields; the fallback is only
### accepted if it is still unique AND the full-field differences it implies are
### exactly CORRECTIONS -- i.e. nothing else moved.
left = try_align(FINGERPRINT)
if left is None and CORRECTIONS:
    corrected_fields = {f for _, f in CORRECTIONS}
    reduced = [f for f in FINGERPRINT if f not in corrected_fields]
    print(f'full-field match failed; retrying without {sorted(corrected_fields)} '
          f'(corrected after 2025-06): matching on {reduced} + Resource/cohort')
    left = try_align(reduced)
    if left is not None:
        f_old, f_new = fingerprint(old), fingerprint(table)
        implied = {(table.index[l], f) for o, l in zip(old.index, left)
                   for f in FINGERPRINT if f_old.loc[o, f] != f_new.loc[table.index[l], f]}
        if implied != CORRECTIONS:
            raise ValueError('the reduced-key match implies fingerprint changes other than the known '
                             f'corrections: implied {sorted(implied)}, CORRECTIONS {sorted(CORRECTIONS)}')
        print(f'  accepted: implied fingerprint changes are exactly {sorted(implied)}')

if left is None:
    raise ValueError("Table S1_0606 cannot be mapped onto Table S1 by inserted rows only: a patient's "
                     'Age/Sex/NF1/TP53/Resource changed, or rows were removed/reordered -- resolve by hand')
renumber = pd.Series(table.index[left], index=old.index)
inserted = [i for i in table.index if i not in set(renumber)]
moved = renumber[renumber.index != renumber.values]
print(f'Table S1_0606 ({len(old)} rows) -> Table S1 ({len(table)} rows): unique order-preserving match')
print(f'  renumbered ({len(moved)}): ' + ', '.join(f'{a}->{b}' for a, b in moved.items()))
print(f'  inserted in Table S1 ({len(inserted)}): {inserted}, '
      f'DNA methylation = {table.loc[inserted, "DNA methyaltion"].tolist()}')

### ==================== 2. discovery cohort (data/beta.csv) ====================
rel = pd.read_csv(data_dir / 'ID_relationship_36cases_NF1_TP53.csv', usecols=['ID_in_network', 'Original_ID'], dtype={'Original_ID': str})
beta_cols = pd.read_csv(data_dir / 'beta.csv', index_col=0, nrows=0).columns.astype(str)
assert set(beta_cols) == set(rel['Original_ID']) and len(beta_cols) == len(rel)
rows = [{'SP_2025': 'S' + r['ID_in_network'], 'cohort': 'Discovery', 'batch': 'discovery',
         'raw_column': r['Original_ID'], 'Original_ID': r['Original_ID'], 'methy_ID': None, 'GSM': None}
        for _, r in rel.iterrows()]

### ==================== 3. validation cohort (raw val_data/beta_{5,14,55}.csv) ====================
link = pd.read_excel(src_dir / 'Methy_Cell_ID_map_to_yy.xlsx', dtype={'Original_ID': str, 'methy_ID': str})
by_sentrix = link.dropna(subset=['methy_ID']).set_index('methy_ID')
by_original = link.set_index('Original_ID')
assert by_sentrix.index.is_unique and by_original.index.is_unique
gse = pd.read_csv(src_dir / 'GSE270993_sample_info.csv', dtype=str)
sentrix_of_gsm = dict(zip(gse['Sample_Name'], gse['Sentrix_ID'] + '_' + gse['Sentrix_Position']))
raw_columns = pd.read_csv(src_dir / 'raw_val_beta_columns.csv', dtype=str)

for _, r in raw_columns.iterrows():
    batch, col = r['batch'], r['raw_column']
    if batch == 'beta_55':
        ### pathology numbers; this batch also holds the 36 discovery samples, which are not in the map
        hit = by_original.loc[col] if col in by_original.index else None
        original_id, gsm = col, None
        sentrix = None if hit is None else hit['methy_ID']
    else:
        gsm = col if batch == 'beta_14' else None
        sentrix = sentrix_of_gsm[col] if gsm else col
        hit = by_sentrix.loc[sentrix] if sentrix in by_sentrix.index else None
        original_id = None if hit is None else hit['Original_ID']
    if hit is None:
        continue
    rows.append({'SP_2025': hit['Cell ID'], 'cohort': 'Validation', 'batch': batch, 'raw_column': col,
                 'Original_ID': original_id, 'methy_ID': sentrix, 'GSM': gsm})

man = pd.DataFrame(rows)
if man['SP_2025'].duplicated().any():
    raise ValueError(f'linked to more than one methylation column: {man.loc[man.SP_2025.duplicated(), "SP_2025"].tolist()}')
linked = set(link.dropna(subset=['methy_ID'])['Cell ID'])
found = set(man.loc[man['cohort'] == 'Validation', 'SP_2025'])
if linked != found:
    raise ValueError(f'Methy_Cell_ID_map_to_yy: Sentrix ID but no raw column for {sorted(linked - found)}; '
                     f'raw column but no Sentrix ID for {sorted(found - linked)}')

man['SP_ID'] = renumber.loc[man['SP_2025']].values
man['beta_column'] = man['raw_column'].where(man['cohort'] == 'Discovery', man['SP_ID'])
manifest = pd.concat([
    man[['SP_ID', 'cohort', 'batch', 'beta_column', 'raw_column', 'Original_ID', 'methy_ID', 'GSM']],
    man['SP_2025'].rename('SP_ID_2025_06'),
    old.loc[man['SP_2025'], FINGERPRINT].reset_index(drop=True),
], axis=1).sort_values('SP_ID').reset_index(drop=True)

if set(manifest['SP_ID']) & set(inserted):
    raise ValueError(f'rows inserted after 2025-06 received methylation data: {sorted(set(manifest.SP_ID) & set(inserted))}')

### ==================== 4. the samples with no methylation data ====================
### The manifest lists every Table S1 sample, so that anything keyed by SP ID can
### be joined against it directly. The rows added here have no methylation data,
### so every column describing that data stays empty -- batch, beta_column,
### raw_column, Original_ID, methy_ID, GSM -- and so does SP_ID_2025_06 (these
### samples were never linked to a methylation matrix) and the fingerprint (it is
### recorded at linking time, and there was no linking for them). Only SP_ID and
### cohort are filled, both read straight off Table S1.
###
### check_manifest() therefore decides "has methylation data" by batch being set,
### never by a row being present.
missing = [sp for sp in table.index if sp not in set(manifest['SP_ID'])]
no_data = pd.DataFrame({'SP_ID': missing,
                        'cohort': table.loc[missing, 'Discovery/Validation'].values})
manifest = pd.concat([manifest, no_data.reindex(columns=manifest.columns)], ignore_index=True)
manifest = manifest.sort_values('SP_ID').reset_index(drop=True)

### ==================== 5. the PL/TL subtype ====================
### Two columns, both read off Table S1:
###   Group/Predict Group (Table S1)   verbatim. 'PL'/'TL' for the discovery
###                                    cohort (the training label), and
###                                    'PL-SCG classifier'/'TL-SCG classifier'
###                                    for the validation cohort (the call).
###                                    Named as in results/subtype_predictions_all_samples.csv,
###                                    where it holds exactly the same strings.
###   subtype                          the same thing as plain PL/TL, with the
###                                    '-SCG classifier' marker dropped (the
###                                    cohort column already says which is which).
###
### Taken from Table S1 rather than from the model's output so that this file
### stays a pure input: 02_predict_all_samples.py reads the manifest, so reading
### its results back in here would make an input depend on an output. The two
### agree for all 86 samples, which is exactly what results/table_s1_group_update.csv
### (0 rows) and the n_changed column of results/group_counts_by_cohort.csv report.
### Prob and the per-feature probabilities live in the results file, not here.
SUBTYPE_SUFFIX = '-SCG classifier'
verbatim = table.loc[manifest['SP_ID'], 'Group/Predict Group'].astype(str).str.strip().values
manifest['Group/Predict Group (Table S1)'] = verbatim
manifest['subtype'] = [v[:-len(SUBTYPE_SUFFIX)] if v.endswith(SUBTYPE_SUFFIX) else v for v in verbatim]

assert manifest['SP_ID'].tolist() == table.index.tolist(), 'manifest is not exactly Table S1, in order'
assert manifest['batch'].notna().sum() == len(man), f'{len(man)} methylation samples expected to carry a batch'
assert manifest.loc[manifest['batch'].isna(), ['beta_column', 'raw_column', 'Original_ID', 'methy_ID',
                                               'GSM', 'SP_ID_2025_06'] + FINGERPRINT].isna().all().all(), \
    'a sample without methylation data carries a value in a methylation/linking column'
bad = set(manifest['subtype']) - {'PL', 'TL'}
assert not bad, f'Table S1 "Group/Predict Group" has values that are not PL/TL even after stripping: {sorted(bad)}'
assert (manifest.loc[manifest['cohort'] == 'Discovery', 'Group/Predict Group (Table S1)']
        == manifest.loc[manifest['cohort'] == 'Discovery', 'subtype']).all(), \
    'a discovery sample carries the "-SCG classifier" marker'

manifest.to_csv(out_path, index=False)
print(f'{len(manifest)} rows = all of Table S1; {manifest["batch"].notna().sum()} with methylation data, '
      f'{manifest["batch"].isna().sum()} without ({", ".join(missing)})')
print('\nsubtype from Table S1:')
print(manifest.groupby(['cohort', 'subtype']).size().rename('n').to_string())
print(f'  all 86: {manifest["subtype"].value_counts().to_dict()}')
print(manifest.groupby(['cohort', 'batch']).size().to_string())
print(f'saved -> {out_path}')
check_manifest(load_manifest(out_path), table)
