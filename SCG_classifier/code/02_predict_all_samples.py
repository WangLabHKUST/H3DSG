"""Apply the models trained by 01_train_model.py to all 86 samples in
Table S1, produce the final PL/TL call for each one, and list every Table S1
'Group/Predict Group' cell that differs from it.

Methy source per cohort:
  - 36 discovery samples : LOOV value from training_methy_loov.csv (each
                            sample scored by models that never saw it)
  - 37 validation samples: model_dict_cpg_alldata.npy applied to
                            data/val_data/beta_{5,14,55}_filtered.csv
                            (columns = current SP IDs, built by 00c from
                            data/sample_manifest.csv)
Every SP-ID-keyed join is guarded by patient_id.check_manifest() (patient
fingerprint Age/Sex/NF1/TP53 + "DNA methyaltion = Yes" set).

Prob is the mean of the available feature probabilities among Age, TP53,
NF1, Ki67 and Methy: a sample's methylation is used whenever it has
methylation data. A sample with none of the 5 has no Prob and its call is NA.

results/subtype_predictions_all_samples.csv, one row per sample:
  Discovery/Validation, Resource, Age, Sex, NF1, TP53, Ki67   as in Table S1 (the values used)
  Methylation availability                                   Table S1 'DNA methyaltion'
  Age_prob, TP53_prob, NF1_prob, Ki67_prob, Methy_prob       P(TL) per feature
  n_features, Prob, Prob_bi_str                              combined probability and call
  Group/Predict Group (Table S1)                             current Table S1 value
  formatted cls                                              what that cell should say:
      discovery -> training label (PL / TL); validation -> call + '-SCG classifier', or NA
results/subtype_predictions_all_samples_by_feature.csv: the same columns reordered so each
  value sits next to its probability (Age, Age_prob, TP53, TP53_prob, ..., Methylation
  availability, Methy_prob, n_features, Prob, Prob_bi_str).
results/table_s1_group_update.csv: the rows where the last two columns differ.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from patient_id import FINGERPRINT, check_batch_columns, check_manifest, load_manifest, load_table_s1

BASE = Path(__file__).resolve().parent.parent
data_dir = BASE / 'data'
models_dir = BASE / 'models'
results_dir = BASE / 'results'
results_dir.mkdir(exist_ok=True)

out_path = results_dir / 'subtype_predictions_all_samples.csv'
reordered_path = results_dir / 'subtype_predictions_all_samples_by_feature.csv'
update_path = results_dir / 'table_s1_group_update.csv'
for p in (out_path, reordered_path, update_path):
    if p.exists():
        raise FileExistsError(f'{p} already exists, not overwriting')

meta_features = ['Age', 'TP53', 'NF1', 'Ki67']
prob_features = meta_features + ['Methy']
model_dict_meta = np.load(models_dir / 'model_dict_meta.npy', allow_pickle=True).item()
model_dict_cpg = np.load(models_dir / 'model_dict_cpg_alldata.npy', allow_pickle=True).item()

overall_df = load_table_s1(data_dir / 'Table S1.xlsx')
manifest = load_manifest(data_dir / 'sample_manifest.csv')
check_manifest(manifest, overall_df)

used_df = overall_df[meta_features].copy()
used_df['TP53'] = used_df['TP53'].str.upper().replace({'WT': 0, 'MUT': 1}).astype(float)
used_df['NF1'] = used_df['NF1'].str.upper().replace({'WT': 0, 'MUT': 1}).astype(float)

### ---- meta-feature probabilities, all 86 samples ----
pred_df = pd.DataFrame(index=overall_df.index, columns=meta_features, dtype=float)
for f in meta_features:
    x = used_df[f].values.reshape(-1, 1)
    valid = ~np.isnan(x).ravel()
    pred_df.loc[valid, f] = model_dict_meta[f].predict_proba(x[valid])[:, 1]

### ---- Methy: LOOV for the 36 discovery samples ----
training_methy = pd.read_csv(models_dir / 'training_methy_loov.csv', index_col=0)['Methy']
assert set(training_methy.index) == set(manifest.loc[manifest['batch'] == 'discovery', 'SP_ID'])
pred_df.loc[training_methy.index, 'Methy'] = training_methy

### ---- Methy: alldata model applied to the 37 validation samples ----
top10k = set(model_dict_cpg)
val_methy = {}
for n in (5, 14, 55):
    batch = f'beta_{n}'
    df = pd.read_csv(data_dir / 'val_data' / f'{batch}_filtered.csv', index_col=0)
    check_batch_columns(df.columns, manifest, batch)
    assert set(df.index) <= top10k, f'{batch}: CpGs outside the trained top-10000 found'
    assert not set(df.columns) & set(val_methy), f'{batch}: sample already seen in another batch'
    probs = pd.DataFrame(index=df.index, columns=df.columns, dtype=float)
    for cpg in df.index:
        probs.loc[cpg] = model_dict_cpg[cpg].predict_proba(df.loc[cpg].values.reshape(-1, 1))[:, 1]
    val_methy.update(probs.mean(axis=0).to_dict())
val_methy = pd.Series(val_methy, name='Methy')
pred_df.loc[val_methy.index, 'Methy'] = val_methy

### ---- combine: every available feature, methylation included whenever the sample has it ----
pred_df['n_features'] = pred_df[prob_features].notna().sum(axis=1)
pred_df['Prob'] = pred_df[prob_features].mean(axis=1, skipna=True)
pred_df['Prob_bi_str'] = np.where(pred_df['Prob'].isna(), 'NA', np.where(pred_df['Prob'] > 0.5, 'TL', 'PL'))

is_discovery = overall_df['Discovery/Validation'] == 'Discovery'
### the discovery cohort's PL/TL label is Table S1's own value for those rows
train_label = overall_df.loc[is_discovery, 'Group/Predict Group'].astype(str).str.strip()

formatted = pd.Series('NA', index=pred_df.index, dtype=object)
formatted[is_discovery] = train_label.loc[pred_df.index[is_discovery]].values
has_call = ~is_discovery & (pred_df['Prob_bi_str'] != 'NA')
formatted[has_call] = pred_df.loc[has_call, 'Prob_bi_str'] + '-SCG classifier'
pred_df['formatted cls'] = formatted
n_disc_diff = (pred_df.loc[is_discovery, 'Prob_bi_str'] != formatted[is_discovery]).sum()
print(f'discovery: call (in-sample meta + LOOV Methy) differs from the training label for {n_disc_diff} / {is_discovery.sum()}')

### ---- output: Table S1 values used, methylation availability, probabilities, call ----
methy_available = overall_df['DNA methyaltion'].astype(str).str.strip()
assert ((methy_available.str.lower() == 'yes') == pred_df['Methy'].notna()).all()
current = overall_df['Group/Predict Group'].fillna('NA').astype(str).str.strip()
out_df = pd.concat([
    overall_df[['Discovery/Validation', 'Resource'] + FINGERPRINT + ['Ki67']],
    methy_available.rename('Methylation availability'),
    pred_df[prob_features].add_suffix('_prob'),
    pred_df[['n_features', 'Prob', 'Prob_bi_str']],
    current.rename('Group/Predict Group (Table S1)'),
    pred_df['formatted cls'],
], axis=1)
out_df.to_csv(out_path)
print(f'saved -> {out_path}')

# same content, each value placed next to its probability
by_feature_cols = (['Discovery/Validation', 'Resource', 'Sex']
                   + [c for f in meta_features for c in (f, f'{f}_prob')]
                   + ['Methylation availability', 'Methy_prob', 'n_features', 'Prob', 'Prob_bi_str',
                      'Group/Predict Group (Table S1)', 'formatted cls'])
assert sorted(by_feature_cols) == sorted(out_df.columns)
out_df[by_feature_cols].to_csv(reordered_path)
print(f'saved -> {reordered_path}')

### ---- Table S1 'Group/Predict Group' cells to update ----
changed = out_df.index[out_df['Group/Predict Group (Table S1)'] != out_df['formatted cls']]
update = (out_df.loc[changed].drop(columns='Prob_bi_str')
          .rename(columns={'Group/Predict Group (Table S1)': 'Table S1 (current)', 'formatted cls': 'proposed'}))
update['reason'] = np.where(update['proposed'] == 'NA', 'none of Age/TP53/NF1/Ki67/Methy available -> no call',
                            np.where(update['Discovery/Validation'] == 'Discovery', 'differs from training label',
                                     'differs from classifier call'))
update.to_csv(update_path)
print(f'\n{len(update)} Table S1 "Group/Predict Group" cells to update:')
print(update.to_string() if len(update) else '  none')
print(f'saved -> {update_path}')
