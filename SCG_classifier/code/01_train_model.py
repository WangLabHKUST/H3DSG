"""Train, from raw data only, every model behind the current PL/TL subtype
classifier (Table S1 + BernoulliNB/GaussianNB meta features + GaussianNB
per-CpG methylation model), so the whole pipeline can be reproduced by
someone who only has this `archive/` folder.

Inputs (all under ../data):
  - Table S1.xlsx                         clinical metadata, 86 samples
  - beta.csv                               141824 CpGs x 36 discovery-cohort
                                            beta values
  - beta_score_full.csv                    per-CpG ranking (score, auc) --
                                            fixed feature-selection input,
                                            independent of the metadata table
  - ID_relationship_36cases_NF1_TP53.csv   Original_ID <-> P-number mapping
                                            (mapping only: its Group/NF1/TP53
                                            columns are not model inputs)
  - sample_manifest.csv                    methylation column -> current SP
                                            ID + patient fingerprint (00b)

Outputs (under ../models):
  - model_dict_meta.npy            {feature: fitted BernoulliNB/GaussianNB}
                                    for Age, TP53, NF1, Ki67
  - model_dict_cpg_alldata.npy     {cpg: fitted GaussianNB}, one fit per top-
                                    10000 CpG on ALL 36 discovery samples --
                                    this is what gets applied to new/
                                    validation samples
  - training_methy_loov.csv        leave-one-out out-of-sample P(TL) for
                                    each of the 36 discovery samples,
                                    averaged over the top-10000 CpGs (the
                                    "Methy" value reported for the training
                                    cohort in the final table -- LOOV is used
                                    here, instead of the alldata model,
                                    specifically so these samples aren't
                                    scored with a model that was fit on them)

Before anything is fit: patient_id.check_manifest() confirms every
methylation sample's recorded Age/Sex/NF1/TP53 matches its Table S1 row.
Every value the models are fitted on -- Age, TP53, NF1, Ki67 and the PL/TL
label -- is read from Table S1.

Model choice per meta feature: BernoulliNB if the feature only takes 2
distinct values in the training data (TP53, NF1 -- WT=0/Mut=1), GaussianNB
otherwise (Age, Ki67). This matches
codes/20260514_updated_model_ki67/01_train_updated_data_n_model_TableR1.py.

GaussianNB fitting is deterministic (closed-form mean/variance, no random
state), so re-running this script always reproduces the same models.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.naive_bayes import GaussianNB, BernoulliNB
from sklearn.model_selection import LeaveOneOut

sys.path.insert(0, str(Path(__file__).resolve().parent))
from patient_id import N_DISCOVERY, check_manifest, load_manifest, load_table_s1

BASE = Path(__file__).resolve().parent.parent
data_dir = BASE / 'data'
models_dir = BASE / 'models'
models_dir.mkdir(exist_ok=True)

for name in ('model_dict_meta.npy', 'model_dict_cpg_alldata.npy', 'training_methy_loov.csv'):
    if (models_dir / name).exists():
        raise FileExistsError(f'{models_dir / name} already exists, not overwriting')

### ==================== 0. patient identity ====================
### sorted by SP ID; fails unless IDs are SP01..SP86 and SP01-SP36 are exactly the Discovery rows
overall_df = load_table_s1(data_dir / 'Table S1.xlsx')
manifest = load_manifest(data_dir / 'sample_manifest.csv')
check_manifest(manifest, overall_df)

### ==================== 1. meta-feature models (Age/TP53/NF1/Ki67) ====================
meta_features = ['Age', 'TP53', 'NF1', 'Ki67']
training_df = overall_df.iloc[:N_DISCOVERY, :].copy()
training_df[['TP53', 'NF1']] = training_df[['TP53', 'NF1']].replace({'WT': 0, 'Mut': 1}).astype(float)
y_meta = training_df['Group/Predict Group'].map({'PL': 0, 'TL': 1}).values
assert not np.isnan(y_meta).any()

model_dict_meta = {}
for f in meta_features:
    x = training_df[f].values.astype(float)
    valid = ~np.isnan(x)
    x = x[valid].reshape(-1, 1)
    y_f = y_meta[valid]
    model = BernoulliNB() if np.unique(x).shape[0] == 2 else GaussianNB()
    model.fit(x, y_f)
    model_dict_meta[f] = model
    print(f'{f}: {type(model).__name__} fit on {len(x)} discovery samples')

np.save(models_dir / 'model_dict_meta.npy', model_dict_meta)
print(f'saved -> {models_dir / "model_dict_meta.npy"}')

### ==================== 2. CpG-level methylation models ====================
beta = pd.read_csv(data_dir / 'beta.csv', index_col=0)
beta.columns = beta.columns.astype(str)

### ID_relationship is read for ID mapping only -- its Group/NF1/TP53 columns are
### not model inputs, so only the two mapping columns are loaded at all.
sample_info = pd.read_csv(data_dir / 'ID_relationship_36cases_NF1_TP53.csv', usecols=['ID_in_network', 'Original_ID'], dtype={'Original_ID': str})
sample_info = sample_info.set_index('Original_ID')
assert set(beta.columns) == set(sample_info.index)

### SP ID of each beta.csv column comes from the fingerprint-checked manifest,
### and must agree with ID_relationship's own P-number
discovery_manifest = manifest[manifest['batch'] == 'discovery'].set_index('beta_column')
assert set(beta.columns) == set(discovery_manifest.index)
sp_id_of = discovery_manifest.loc[beta.columns, 'SP_ID']
assert (sp_id_of.values == ('S' + sample_info.loc[beta.columns, 'ID_in_network']).values).all()

### the CpG models' PL/TL labels come from Table S1, the single source of truth
y_cpg = overall_df.loc[sp_id_of, 'Group/Predict Group'].str.strip().map({'PL': 0, 'TL': 1}).values
assert not np.isnan(y_cpg).any()
print(f'CpG labels from Table S1: {int((y_cpg == 0).sum())} PL, {int((y_cpg == 1).sum())} TL')

beta_score = pd.read_csv(data_dir / 'beta_score_full.csv', index_col=0)
top10k = beta_score.sort_values(by='score', ascending=False).index.tolist()[:10000]
beta_top = beta.loc[top10k]

print(f'training LOOV + alldata GaussianNB on {len(top10k)} CpGs (this takes a few minutes) ...')
loo = LeaveOneOut()
loov_probs = pd.DataFrame(index=top10k, columns=beta.columns, dtype=float)
model_dict_cpg_alldata = {}
for i, cpg in enumerate(top10k):
    x = beta_top.loc[cpg].values.reshape(-1, 1)
    ### leave-one-out: used only to report an honest out-of-sample Methy
    ### value for the discovery cohort itself
    for train_idx, test_idx in loo.split(x):
        m = GaussianNB()
        m.fit(x[train_idx], y_cpg[train_idx])
        loov_probs.loc[cpg, beta.columns[test_idx[0]]] = m.predict_proba(x[test_idx])[0, 1]
    ### alldata: the model actually applied to new/validation samples
    model_dict_cpg_alldata[cpg] = GaussianNB().fit(x, y_cpg)
    if (i + 1) % 2000 == 0:
        print(f'  {i + 1} / {len(top10k)} CpGs done')

np.save(models_dir / 'model_dict_cpg_alldata.npy', model_dict_cpg_alldata)
print(f'saved -> {models_dir / "model_dict_cpg_alldata.npy"}')

training_methy = loov_probs.mean(axis=0)
training_methy.index = sp_id_of.loc[training_methy.index].values
training_methy = training_methy.rename('Methy').to_frame()
training_methy.to_csv(models_dir / 'training_methy_loov.csv')
print(f'saved -> {models_dir / "training_methy_loov.csv"}')
print('\ndone.')
