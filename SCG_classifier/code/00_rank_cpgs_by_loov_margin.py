"""Rank every CpG in data/beta.csv by how well it alone classifies the 36
discovery-cohort patients, producing data/beta_score_full.csv -- the fixed
ranking that 01_train_model.py's top-10000 feature selection is based on.

For each CpG: 36-fold leave-one-out GaussianNB (that CpG's beta value as the
sole feature), scored by
  - score = mean LOOV margin, P(correct class) - P(incorrect class)
  - auc   = ROC-AUC of the LOOV P(TL) probabilities against the true label
Mirrors codes/20250609_try.py (score) + codes/20250612_training_models.py
(auc). The choice of keeping the top 10000 CpGs by 'score' (rather than some
other cutoff) was made by comparing single- and combined-feature margin/AUC
as a function of how many top-ranked CpGs are averaged together
(codes/20250612_comb_after_ranking_prog_cont_v2.py) and picking the point
where combined performance plateaus.

This is the single expensive step in the whole pipeline (141824 CpGs x 36
LOOV folds, ~1 hour) and is NOT required to reproduce 01_train_model.py /
02_predict_all_samples.py -- beta_score_full.csv already ships in data/.
Run this only to verify the ranking itself is reproducible from beta.csv.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut
from sklearn.naive_bayes import GaussianNB

sys.path.insert(0, str(Path(__file__).resolve().parent))
from patient_id import check_manifest, load_manifest, load_table_s1

BASE = Path(__file__).resolve().parent.parent
data_dir = BASE / 'data'

out_path = data_dir / 'beta_score_full.csv'
if out_path.exists():
    raise FileExistsError(f'{out_path} already exists, not overwriting')

beta = pd.read_csv(data_dir / 'beta.csv', index_col=0)
beta.columns = beta.columns.astype(str)

### ID_relationship is read for ID mapping only -- the PL/TL labels this ranking
### is built on come from Table S1, like every other value the model is fitted on.
overall_df = load_table_s1(data_dir / 'Table S1.xlsx')
manifest = load_manifest(data_dir / 'sample_manifest.csv')
check_manifest(manifest, overall_df)

sample_info = pd.read_csv(data_dir / 'ID_relationship_36cases_NF1_TP53.csv',
                          usecols=['ID_in_network', 'Original_ID'], dtype={'Original_ID': str})
sample_info = sample_info.set_index('Original_ID')
assert set(beta.columns) == set(sample_info.index)
discovery_manifest = manifest[manifest['batch'] == 'discovery'].set_index('beta_column')
sp_id_of = discovery_manifest.loc[beta.columns, 'SP_ID']
assert (sp_id_of.values == ('S' + sample_info.loc[beta.columns, 'ID_in_network']).values).all()

y = overall_df.loc[sp_id_of, 'Group/Predict Group'].str.strip().map({'PL': 0, 'TL': 1}).values
assert not np.isnan(y).any()

loo = LeaveOneOut()
scores = np.empty(len(beta))
aucs = np.empty(len(beta))
for i, cpg in enumerate(beta.index):
    x = beta.loc[cpg].values.reshape(-1, 1)
    margins, probs = [], []
    for train_idx, test_idx in loo.split(x):
        model = GaussianNB().fit(x[train_idx], y[train_idx])
        prob = model.predict_proba(x[test_idx])[0]
        label = y[test_idx][0]
        margins.append(prob[label] - prob[1 - label])
        probs.append(prob[1])
    scores[i] = np.mean(margins)
    aucs[i] = roc_auc_score(y, probs)
    if (i + 1) % 10000 == 0:
        print(f'{i + 1} / {len(beta)} CpGs scored')

beta_score = pd.DataFrame({'score': scores, 'auc': aucs}, index=beta.index)
beta_score.to_csv(out_path)
print(f'saved -> {out_path}')
