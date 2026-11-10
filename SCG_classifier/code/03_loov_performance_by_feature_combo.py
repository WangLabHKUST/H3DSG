"""Leave-one-out classification performance (Accuracy/F1/Sensitivity/
Specificity) on the 36 discovery-cohort patients, for each single feature
and feature combination in FEATURE_SETS. Computed from the current data
only -- not aligned to any earlier table (res_df_0623.csv / Table S4).

Each feature is scored strictly out-of-sample:
  - Age/TP53/NF1/Ki67: refit per LOO fold here (36 folds), same model choice
    as 01_train_model.py (BernoulliNB for the 2-valued TP53/NF1, GaussianNB
    otherwise).
  - Methy: models/training_methy_loov.csv, produced by 01_train_model.py's
    own 36-fold LOOV over the top-10000 CpGs.
A combination's probability is the mean of its features' LOOV probabilities,
thresholded at 0.5; samples missing any feature of the combination are
dropped for that row (see n).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import LeaveOneOut
from sklearn.naive_bayes import BernoulliNB, GaussianNB

sys.path.insert(0, str(Path(__file__).resolve().parent))
from patient_id import N_DISCOVERY, check_manifest, load_manifest, load_table_s1

BASE = Path(__file__).resolve().parent.parent
data_dir = BASE / 'data'
models_dir = BASE / 'models'
results_dir = BASE / 'results'
results_dir.mkdir(exist_ok=True)

out_path = results_dir / 'loov_performance_by_feature_combo.csv'
if out_path.exists():
    raise FileExistsError(f'{out_path} already exists, not overwriting')

meta_features = ['Age', 'TP53', 'NF1', 'Ki67']
FEATURE_SETS = [
    ('Methy',), ('Age',), ('TP53',), ('Ki67',), ('NF1',),
    ('Methy', 'Age', 'Ki67'), ('Methy', 'Age', 'TP53'), ('Methy', 'Age', 'NF1'),
    ('Methy', 'Ki67', 'TP53'), ('Methy', 'Ki67', 'NF1'), ('Methy', 'NF1', 'TP53'),
    ('Age', 'Ki67', 'TP53'), ('Age', 'Ki67', 'NF1'), ('Ki67', 'TP53', 'NF1'),
    ('Methy', 'Age', 'TP53', 'NF1', 'Ki67'),
]

overall_df = load_table_s1(data_dir / 'Table S1.xlsx')
manifest = load_manifest(data_dir / 'sample_manifest.csv')
check_manifest(manifest, overall_df)

### ---- 1. LOOV over the 36 discovery patients for Age/TP53/NF1/Ki67 ----
training_df = overall_df.iloc[:N_DISCOVERY, :].copy()
training_df['TP53'] = training_df['TP53'].str.upper().map({'WT': 0, 'MUT': 1}).astype(float)
training_df['NF1'] = training_df['NF1'].str.upper().map({'WT': 0, 'MUT': 1}).astype(float)
y = training_df['Group/Predict Group'].map({'PL': 0, 'TL': 1})
assert y.notna().all()

loo = LeaveOneOut()
meta_loov_df = pd.DataFrame(index=training_df.index, columns=meta_features, dtype=float)
for f in meta_features:
    vals = training_df[f].astype(float)
    valid_idx = np.where(vals.notna().values)[0]
    valid_ids = training_df.index[valid_idx]
    x = vals.iloc[valid_idx].values.reshape(-1, 1)
    y_f = y.iloc[valid_idx].values
    model_cls = BernoulliNB if np.unique(x).shape[0] == 2 else GaussianNB
    for train_idx, test_idx in loo.split(x):
        model = model_cls()
        model.fit(x[train_idx], y_f[train_idx])
        meta_loov_df.loc[valid_ids[test_idx[0]], f] = model.predict_proba(x[test_idx])[0, 1]

### ---- 2. bring in the Methy LOOV from 01_train_model.py ----
methy_loov = pd.read_csv(models_dir / 'training_methy_loov.csv', index_col=0)['Methy']
assert set(methy_loov.index) == set(training_df.index) == set(manifest.loc[manifest['batch'] == 'discovery', 'SP_ID'])

inference_df = pd.merge(meta_loov_df, methy_loov.to_frame(), left_index=True, right_index=True, how='inner')
inference_df['label'] = y.loc[inference_df.index].values

### ---- 3. metrics per feature combination ----
rows = []
for features in FEATURE_SETS:
    this_df = inference_df[list(features)].dropna()
    this_label = inference_df.loc[this_df.index, 'label']
    this_pred = (this_df.mean(axis=1).values > 0.5).astype(int)

    report = classification_report(this_label, this_pred, output_dict=True, zero_division=0)
    cm = confusion_matrix(this_label, this_pred, labels=[0, 1])
    rows.append({
        'Features': '+'.join(features),
        'Accuracy': report['accuracy'],
        'F1 score': report['macro avg']['f1-score'],
        'Sensitivity': cm[1, 1] / (cm[1, 1] + cm[1, 0]) if (cm[1, 1] + cm[1, 0]) else np.nan,
        'Specificity': cm[0, 0] / (cm[0, 0] + cm[0, 1]) if (cm[0, 0] + cm[0, 1]) else np.nan,
        'n': len(this_df),
    })
res_df = pd.DataFrame(rows)

res_df.to_csv(out_path)
print(f'\n{res_df.to_string()}')
print(f'\nsaved -> {out_path}')
