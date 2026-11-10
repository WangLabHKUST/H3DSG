"""Write results/group_counts_by_cohort.csv: how many TL / PL / NA there are in
each cohort, both as Table S1 currently says and as the current model says.
Meant as a checkpoint when drawing figures -- a plot of the validation cohort
should show exactly these numbers.

Rows: All samples; Discovery; Validation; Validation split by Resource.
Columns: n, TL/PL/NA per source ('model' = 'formatted cls', i.e. the training
label for discovery and the classifier call for validation; 'Table S1' = the
current 'Group/Predict Group' cell), and n_changed = rows where they differ.
Input: results/subtype_predictions_all_samples.csv (02_predict_all_samples.py).
"""
import sys
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
pred_path = BASE / 'results' / 'subtype_predictions_all_samples.csv'
out_path = BASE / 'results' / 'group_counts_by_cohort.csv'
if out_path.exists():
    raise FileExistsError(f'{out_path} already exists, not overwriting')

pred = pd.read_csv(pred_path, index_col=0, keep_default_na=False)
strip = lambda s: s.str.replace('-SCG classifier', '', regex=False).str.strip()
model, table_s1 = strip(pred['formatted cls']), strip(pred['Group/Predict Group (Table S1)'])
groups = ['TL', 'PL', 'NA']
unexpected = sorted((set(model) | set(table_s1)) - set(groups))
if unexpected:
    raise ValueError(f'unexpected group values: {unexpected}')

is_val = pred['Discovery/Validation'] == 'Validation'
subsets = {'All samples': pred.index,
           'Discovery': pred.index[~is_val],
           'Validation': pred.index[is_val]}
for res in sorted(pred.loc[is_val, 'Resource'].unique()):
    subsets[f'Validation / {res}'] = pred.index[is_val & (pred['Resource'] == res)]

rows = []
for name, idx in subsets.items():
    row = {'subset': name, 'n': len(idx)}
    row.update({f'{g} (model)': int((model[idx] == g).sum()) for g in groups})
    row.update({f'{g} (Table S1)': int((table_s1[idx] == g).sum()) for g in groups})
    row['n_changed'] = int((model[idx] != table_s1[idx]).sum())
    rows.append(row)
out = pd.DataFrame(rows).set_index('subset')

assert (out.loc['Discovery', 'n'] + out.loc['Validation', 'n']) == out.loc['All samples', 'n']
val_rows = [k for k in subsets if k.startswith('Validation / ')]
assert out.loc[val_rows, 'n'].sum() == out.loc['Validation', 'n']
for src in ('model', 'Table S1'):
    assert (out[[f'{g} ({src})' for g in groups]].sum(axis=1) == out['n']).all()

out.to_csv(out_path)
print(out.to_string())
print(f'\nsamples where the two differ: {sorted(pred.index[model != table_s1])}')
print(f'saved -> {out_path}')
