"""Independent end-to-end re-check of the archive's inputs and outputs.

Written without importing anything from code/ and without scikit-learn:
Naive Bayes is re-implemented from its closed form (sklearn's defaults:
GaussianNB var_smoothing=1e-9, BernoulliNB alpha=1), the 2025-06 -> Table S1
renumbering is re-derived from different table fields than 00b uses, and the
validation Methy values are recomputed from the raw beta matrices.

Milestones (each prints PASS/FAIL with its evidence):
  M0  inputs   : archive data/Table S1.xlsx == the collaborator's latest
                 clinical table (since 2026-09-16 that file is Table S1.xlsx)
  M1  identity : which patient every methylation column belongs to
                 (ID_relationship is used for ID mapping only)
  M2  features : meta models == closed-form fit on Table S1 SP01-SP36 values
  M3  modality : CpG ranking, CpG models and Methy recomputed; every sample
                 uses exactly the features Table S1 provides
  M4  calls    : Prob / PL-TL recomputed; Table S1 cells to change; the
                 results table carries the Table S1 values it used

Usage: python code/04_independent_check.py [TABLE_S1] [RAW_VAL_DIR]
Report: results/independent_check_report.txt (overwritten on every run).
Exit code 1 if any milestone fails.
"""
import datetime
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 30)

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / 'data'
TABLE_S1 = Path(sys.argv[1]) if len(sys.argv) > 1 else DATA / 'Table S1.xlsx'
RAW = Path(sys.argv[2]) if len(sys.argv) > 2 else BASE.parent.parent / 'val_data'
OTHER_DATA = BASE.parent.parent / 'data' / 'other_data.csv'
REPORT = BASE / 'results' / 'independent_check_report.txt'
FEATURES = ['Age', 'TP53', 'NF1', 'Ki67', 'Methy']

out, status = [], {}


def log(s=''):
    print(s)
    out.append(str(s))


def milestone(name, ok, evidence):
    status[name] = bool(ok)
    log(f'[{"PASS" if ok else "FAIL"}] {name}\n        {evidence}')


def md5(p):
    return hashlib.md5(p.read_bytes()).hexdigest()


def mtime(p):
    return datetime.datetime.fromtimestamp(p.stat().st_mtime).strftime('%Y-%m-%d %H:%M')


def norm(v, numeric=False):
    if pd.isna(v):
        return 'NA'
    return f'{round(float(v), 6):g}' if numeric else str(v).strip().lower()


### ---- closed-form Naive Bayes ----
def gnb_fit(x, y):
    """GaussianNB(var_smoothing=1e-9) fitted independently on every row of x (rows = features, columns = samples)."""
    eps = 1e-9 * x.var(axis=1)
    return [((y == c).sum() / len(y), x[:, y == c].mean(axis=1), x[:, y == c].var(axis=1) + eps) for c in (0, 1)]


def gnb_p_tl(params, x):
    jll = [np.log(prior) - 0.5 * np.log(2 * np.pi * var)[:, None] - 0.5 * (x - mu[:, None]) ** 2 / var[:, None]
           for prior, mu, var in params]
    return np.exp(jll[1] - np.logaddexp(jll[0], jll[1]))


def bnb_p_tl(x_train, y, x):
    """BernoulliNB(alpha=1, binarize=0) on one 0/1 feature."""
    joint = []
    for c in (0, 1):
        xc = x_train[y == c]
        p_mut = ((xc > 0).sum() + 1) / (len(xc) + 2)
        joint.append(len(xc) / len(x_train) * np.where(x > 0, p_mut, 1 - p_mut))
    return joint[1] / (joint[0] + joint[1])


log(f'Independent check, run {datetime.datetime.now():%Y-%m-%d %H:%M}')
log(f'Table S1: {TABLE_S1}\nraw validation beta: {RAW}')

### ==================== M0 inputs ====================
log('\n== M0  inputs: does the archive use the latest Table S1?')
siblings = sorted((p for p in TABLE_S1.parent.iterdir() if p.is_file()), key=lambda p: p.stat().st_mtime)
newer = [p.name for p in siblings if p.stat().st_mtime > TABLE_S1.stat().st_mtime]
log('  files in the collaborator folder (oldest -> newest): ' + '; '.join(f'{p.name} [{mtime(p)}]' for p in siblings))
# content, not bytes: re-saving in Excel (filter, frozen panes, column widths) changes the md5 only
same_bytes = md5(TABLE_S1) == md5(DATA / 'Table S1.xlsx')
milestone('M0.1 archive data/Table S1.xlsx has the same content as the collaborator file (header, row order, every cell)',
          pd.read_excel(TABLE_S1, header=None).equals(pd.read_excel(DATA / 'Table S1.xlsx', header=None)),
          f'md5 {md5(TABLE_S1)} (collaborator, {mtime(TABLE_S1)}) vs {md5(DATA / "Table S1.xlsx")} (archive), '
          f'{"byte-identical" if same_bytes else "bytes differ"}; files in that folder modified later: {newer or "none"}')
### the published Table S1 carries its title in the first row: locate the header
_head = pd.read_excel(TABLE_S1, header=None, nrows=2)
_hrow = next((i for i in (0, 1) if str(_head.iloc[i, 0]).strip() == 'Sample ID'), None)
if _hrow is None:
    raise ValueError(f'{TABLE_S1}: no "Sample ID" header found in the first two rows')
s1 = pd.read_excel(TABLE_S1, index_col=0, header=_hrow)
s1.index = s1.index.astype(str).str.strip()
s1 = s1.sort_index()
layout_ok = (s1.index.tolist() == [f'SP{i:02d}' for i in range(1, 87)]
             and s1.index[s1['Discovery/Validation'] == 'Discovery'].tolist() == [f'SP{i:02d}' for i in range(1, 37)])
milestone('M0.2 Table S1 layout', layout_ok,
          f'{len(s1)} rows, IDs SP01..SP86 gap-free, Discovery == SP01-SP36, '
          f'{int((s1["Discovery/Validation"] == "Validation").sum())} Validation')

### ==================== M1 identity ====================
log('\n== M1  identity: which patient is each methylation column?')
### mapping only: Group/NF1/TP53 are not read, so they cannot reach the model
rel = pd.read_csv(DATA / 'ID_relationship_36cases_NF1_TP53.csv',
                  usecols=['ID_in_network', 'Multiomics_analysis_ID', 'Original_ID'],
                  dtype={'Original_ID': str})
sp_of_orig = dict(zip(rel['Original_ID'], 'S' + rel['ID_in_network']))
beta = pd.read_csv(DATA / 'beta.csv', index_col=0)
beta.columns = beta.columns.astype(str)


milestone('M1.1 discovery: beta.csv columns -> SP01-SP36 through ID_relationship (mapping only)',
          set(beta.columns) == set(rel['Original_ID']) and sorted(sp_of_orig.values()) == s1.index[:36].tolist(),
          f'36/36 columns mapped to SP01-SP36 via Original_ID -> P-number; ID_relationship supplies no model '
          f'input (its Group/NF1/TP53 columns are not read) -- labels and mutation calls come from Table S1')

raw55 = pd.read_csv(RAW / 'beta_55.csv', index_col=0)
raw55.columns = raw55.columns.astype(str)
raw55.index = raw55.index.astype(str).str.split('_').str[0]
raw55 = raw55[~raw55.index.duplicated(keep='first')]
if all(c in raw55.columns for c in beta.columns):
    common = beta.index.intersection(raw55.index)
    B, Rv = beta.loc[common].values, raw55.loc[common].values
    keep = ~np.isnan(Rv).any(axis=1) & ~np.isnan(B).any(axis=1)
    zs = [(m[keep] - m[keep].mean(axis=0)) / m[keep].std(axis=0) for m in (B, Rv)]
    C = zs[0].T @ zs[1] / keep.sum()
    self_pos = np.array([raw55.columns.get_loc(c) for c in beta.columns])
    self_r = C[np.arange(36), self_pos]
    C_other = C.copy()
    C_other[np.arange(36), self_pos] = -np.inf
    n_self = int((C.argmax(axis=1) == self_pos).sum())
    milestone('M1.2 discovery, data-intrinsic: each beta.csv column correlates best with the same Original_ID in raw beta_55.csv',
              n_self == 36, f'{n_self}/36 over {int(keep.sum())} CpGs; same-sample r >= {self_r.min():.4f}, '
                            f'any other sample r <= {C_other.max():.4f}')
else:
    milestone('M1.2 discovery, data-intrinsic check', False, 'not all discovery Original_IDs are columns of raw beta_55.csv')

old = pd.read_excel(DATA / 'id_mapping' / 'Table S1_0606.xlsx', dtype={'Sample ID': str}).set_index('Sample ID')


def key(df, ki67):
    return df.apply(lambda r: '|'.join([norm(r['Age'], True), norm(r['Sex']), norm(r[ki67], True), norm(r['WES']),
                                        norm(r['Target_seq']), norm(r['RNA-seq']), norm(r['Resource'])]), axis=1)


def tiebreak(df):
    return df.apply(lambda r: f"{norm(r['OS (Months)'], True)}|{norm(r['Sensor'], True)}", axis=1)


k_old, k_new, t_old, t_new = key(old, 'MKI67'), key(s1, 'Ki67'), tiebreak(old), tiebreak(s1)
renum, via_tiebreak, unresolved = {}, [], {}
for o in old.index:
    cand = k_new.index[k_new == k_old[o]].tolist()
    if len(cand) > 1:
        via_tiebreak.append(o)
        cand = [c for c in cand if t_new[c] == t_old[o]]
    if len(cand) == 1:
        renum[o] = cand[0]
    else:
        unresolved[o] = cand

link = pd.read_excel(DATA / 'id_mapping' / 'Methy_Cell_ID_map_to_yy.xlsx', dtype=str)
by_sentrix = link.dropna(subset=['methy_ID']).set_index('methy_ID')['Cell ID']
by_orig = link.set_index('Original_ID')['Cell ID']
gse = pd.read_csv(DATA / 'id_mapping' / 'GSE270993_sample_info.csv', dtype=str)
gsm2sentrix = dict(zip(gse['Sample_Name'], gse['Sentrix_ID'] + '_' + gse['Sentrix_Position']))
headers = {'beta_55': raw55.columns.tolist()}
for n in (5, 14):
    headers[f'beta_{n}'] = pd.read_csv(RAW / f'beta_{n}.csv', index_col=0, nrows=0).columns.astype(str).tolist()
old_of_col = {}
for b, cols in headers.items():
    for c in cols:
        o = by_orig.get(c) if b == 'beta_55' else by_sentrix.get(gsm2sentrix.get(c) if b == 'beta_14' else c)
        if isinstance(o, str):
            old_of_col[(b, c)] = o
needed = set(old_of_col.values()) | set(s1.index[:36])
bad = {o: unresolved[o] for o in needed if o in unresolved}
renumbered = {o: n for o, n in renum.items() if o != n}
milestone('M1.3 2025-06 -> Table S1 renumbering, re-derived from Age/Sex/Ki67/WES/Target_seq/RNA-seq/Resource '
          '(OS/Sensor only to break ties)',
          not bad and len(set(renum.values())) == len(renum),
          f'{len(renum)}/{len(old)} old rows matched to exactly one Table S1 row ({len(via_tiebreak)} needed OS/Sensor); '
          f'unresolved among methylation samples: {bad or "none"} (other unresolved: {sorted(set(unresolved) - set(bad)) or "none"}); '
          f'renumbered: {renumbered}; Table S1 rows with no 2025-06 row: {sorted(set(s1.index) - set(renum.values()))}')

val_cols = {b: {c: renum.get(o, f'?{o}') for (bb, c), o in old_of_col.items() if bb == b} for b in headers}
val_ids = [sp for m in val_cols.values() for sp in m.values()]
have = set(val_ids) | set(s1.index[:36])
yes = set(s1.index[s1['DNA methyaltion'].astype(str).str.strip().str.lower() == 'yes'])
### the manifest lists all 86 Table S1 samples; the 13 without methylation data
### have no batch and are not part of this comparison
man = pd.read_csv(DATA / 'sample_manifest.csv', dtype=str)
man = man[man['batch'].notna()]
man_pairs = set(zip(man['batch'], man['raw_column'], man['SP_ID']))
my_pairs = {(b, c, sp) for b, m in val_cols.items() for c, sp in m.items()} | {('discovery', o, sp) for o, sp in sp_of_orig.items()}
milestone('M1.4 methylation samples: set and SP IDs',
          len(val_ids) == len(set(val_ids)) == 37 and have == yes and my_pairs == man_pairs,
          f'{len(val_ids)} validation columns -> {len(set(val_ids))} distinct SP IDs '
          f'(beta_5 {len(val_cols["beta_5"])}, beta_14 {len(val_cols["beta_14"])}, beta_55 {len(val_cols["beta_55"])}); '
          f'36+37 SP IDs with data == Table S1 "DNA methyaltion = Yes" ({len(yes)}): {have == yes}; '
          f'identical to data/sample_manifest.csv (batch, raw column, SP ID), 73 rows: {my_pairs == man_pairs}')

### ==================== M2 training features ====================
log('\n== M2  training features: are the meta models fitted on the current Table S1 values?')
enc = {'WT': 0.0, 'MUT': 1.0}


def fvals(df, f):
    if f in ('TP53', 'NF1'):
        return df[f].map(lambda v: np.nan if pd.isna(v) else enc[str(v).strip().upper()]).astype(float)
    return df[f].astype(float)


disc = s1.iloc[:36]
y = disc['Group/Predict Group'].map({'PL': 0, 'TL': 1}).values
shipped = np.load(BASE / 'models' / 'model_dict_meta.npy', allow_pickle=True).item()
res = pd.read_csv(BASE / 'results' / 'subtype_predictions_all_samples.csv', index_col=0,
                  keep_default_na=False, na_values=[''])
indep = pd.DataFrame(index=s1.index, columns=FEATURES, dtype=float)
param_diff = {}
for f in ('Age', 'TP53', 'NF1', 'Ki67'):
    xd, xa = fvals(disc, f).values, fvals(s1, f).values
    ok, oka = ~np.isnan(xd), ~np.isnan(xa)
    if f in ('Age', 'Ki67'):
        params = gnb_fit(xd[ok][None, :], y[ok])
        theta = np.array([params[0][1][0], params[1][1][0]])
        var = np.array([params[0][2][0], params[1][2][0]])
        diff = max(np.abs(theta - shipped[f].theta_.ravel()).max(), np.abs(var - shipped[f].var_.ravel()).max())
        indep.loc[oka, f] = gnb_p_tl(params, xa[oka][None, :])[0]
        desc = f'n={int(ok.sum())}, class mean PL {theta[0]:.4g} / TL {theta[1]:.4g}'
    else:
        counts = np.array([(xd[ok][y[ok] == c] > 0).sum() for c in (0, 1)])
        diff = np.abs(counts - shipped[f].feature_count_.ravel()).max()
        indep.loc[oka, f] = bnb_p_tl(xd[ok], y[ok], xa[oka])
        desc = f'Mut in PL {counts[0]}/{int((y[ok] == 0).sum())}, in TL {counts[1]}/{int((y[ok] == 1).sum())}'
    class_ok = np.array_equal(shipped[f].class_count_, [(y[ok] == 0).sum(), (y[ok] == 1).sum()])
    param_diff[f] = diff if class_ok else np.inf
    log(f'  {f}: Table S1 SP01-SP36 {desc}; shipped model parameters max|diff| {diff:.1e}, class counts equal {class_ok}')
milestone('M2.1 shipped meta models == closed-form fit on Table S1 SP01-SP36',
          max(param_diff.values()) < 1e-9, f'max parameter difference {max(param_diff.values()):.1e}')
d_meta = {f: float(np.nanmax(np.abs(indep[f] - res[f'{f}_prob']))) for f in ('Age', 'TP53', 'NF1', 'Ki67')}
nan_same = all((indep[f].isna() == res[f'{f}_prob'].isna()).all() for f in d_meta)
milestone('M2.2 Age/TP53/NF1/Ki67 probabilities in results == Table S1 values through those models, all 86 samples',
          nan_same and max(d_meta.values()) < 1e-9,
          f'max|diff| {max(d_meta.values()):.1e}; missing exactly where Table S1 is empty: {nan_same}')
if OTHER_DATA.exists():
    od = pd.read_csv(OTHER_DATA, index_col=0)
    num_to_sp = {int(m.split('_')[1]): 'S' + p for m, p in zip(rel['Multiomics_analysis_ID'], rel['ID_in_network'])}
    od.index = [num_to_sp[int(str(i).split('_')[1])] for i in od.index]
    od_means = od['KI67'].groupby(disc.loc[od.index, 'Group/Predict Group']).mean()
    log(f'  contrast: Ki67 class means from the inaccurate data/other_data.csv would be PL {od_means["PL"]:.4g} / '
        f'TL {od_means["TL"]:.4g}; shipped model has {shipped["Ki67"].theta_.ravel()[0]:.4g} / '
        f'{shipped["Ki67"].theta_.ravel()[1]:.4g} -> other_data.csv is not what the model uses')

### ==================== M3 methylation modality ====================
log('\n== M3  modality: methylation recomputed, features used vs Table S1')
ylab = s1.loc[[sp_of_orig[c] for c in beta.columns], 'Group/Predict Group'].map({'PL': 0, 'TL': 1}).values
X = beta.values
p_loov = np.empty_like(X)
for j in range(X.shape[1]):
    tr = np.arange(X.shape[1]) != j
    p_loov[:, j] = gnb_p_tl(gnb_fit(X[:, tr], ylab[tr]), X[:, [j]])[:, 0]
margin = np.where(ylab == 1, 2 * p_loov - 1, 1 - 2 * p_loov).mean(axis=1)
score_file = pd.read_csv(DATA / 'beta_score_full.csv', index_col=0)['score'].reindex(beta.index)
top_file = set(score_file.sort_values(ascending=False).index[:10000])
top_mine = set(beta.index[np.argsort(-margin, kind='stable')[:10000]])
cpg_models = np.load(BASE / 'models' / 'model_dict_cpg_alldata.npy', allow_pickle=True).item()
d_score = float(np.abs(margin - score_file.values).max())
milestone('M3.1 CpG ranking recomputed from beta.csv (141824 CpGs x 36-fold LOOV)',
          d_score < 1e-6 and top_mine == top_file == set(cpg_models),
          f'score max|diff| vs beta_score_full.csv {d_score:.1e}; top-10000 identical to the file: {top_mine == top_file}; '
          f'== CpGs of the shipped models: {top_file == set(cpg_models)}')

top = [c for c in beta.index if c in top_file]
ti = beta.index.get_indexer(top)
params = gnb_fit(X[ti], ylab)
th_ship = np.array([cpg_models[c].theta_.ravel() for c in top])
va_ship = np.array([cpg_models[c].var_.ravel() for c in top])
d_cpg = max(np.abs(np.stack([params[0][1], params[1][1]], axis=1) - th_ship).max(),
            np.abs(np.stack([params[0][2], params[1][2]], axis=1) - va_ship).max())
methy_disc = pd.Series(p_loov[ti].mean(axis=0), index=[sp_of_orig[c] for c in beta.columns])
ship_loov = pd.read_csv(BASE / 'models' / 'training_methy_loov.csv', index_col=0)['Methy']
d_loov = float(np.abs(methy_disc - ship_loov.reindex(methy_disc.index)).max())
d_res_disc = float(np.abs(methy_disc - res.loc[methy_disc.index, 'Methy_prob']).max())
milestone('M3.2 per-CpG models and discovery Methy (LOOV) recomputed',
          d_cpg < 1e-9 and d_loov < 1e-9 and d_res_disc < 1e-9,
          f'10000 CpG models max|param diff| {d_cpg:.1e}; discovery Methy vs training_methy_loov.csv {d_loov:.1e}, '
          f'vs results {d_res_disc:.1e}')

top_pos = pd.Index(top)
methy_val, ncpg = {}, {}
for b, m in val_cols.items():
    cols = list(m)
    if b == 'beta_55':
        df = raw55[cols]
    else:
        hdr = headers[b]
        df = pd.read_csv(RAW / f'{b}.csv', index_col=0, usecols=[0] + [hdr.index(c) + 1 for c in cols])
        df.columns = df.columns.astype(str)
        df.index = df.index.astype(str).str.split('_').str[0]
        df = df[~df.index.duplicated(keep='first')][cols]
    df = df.loc[df.index.isin(top_file)]
    ncpg[b] = len(df)
    idx = top_pos.get_indexer(df.index)
    p = gnb_p_tl([(prior, mu[idx], var[idx]) for prior, mu, var in params], df.values)
    for k, c in enumerate(cols):
        methy_val[m[c]] = float(np.nanmean(p[:, k]))
methy_val = pd.Series(methy_val)
d_val = float(np.abs(methy_val - res.loc[methy_val.index, 'Methy_prob']).max())
no_methy = sorted(set(s1.index) - have)
milestone('M3.3 validation Methy recomputed from the raw beta matrices',
          d_val < 1e-9 and res.loc[no_methy, 'Methy_prob'].isna().all(),
          f'{len(methy_val)} samples, max|diff| vs results {d_val:.1e}; CpGs used per batch {ncpg}; '
          f'no Methy in results for the {len(no_methy)} samples without methylation data: {", ".join(no_methy)}')


def yes_col(c):
    return s1[c].astype(str).str.strip().str.lower().str.startswith('yes')


expected = pd.DataFrame({'Age': s1['Age'].notna(), 'TP53': s1['TP53'].notna(), 'NF1': s1['NF1'].notna(),
                         'Ki67': s1['Ki67'].notna(), 'Methy': yes_col('DNA methyaltion')})
used = res[[f'{f}_prob' for f in FEATURES]].notna().set_axis(FEATURES, axis=1)
mism = expected.index[(expected != used).any(axis=1)].tolist()
combo = used.apply(lambda r: '+'.join(f for f in FEATURES if r[f]) or 'none', axis=1)
milestone('M3.4 every sample uses exactly the features Table S1 provides (value present; Methy iff DNA methylation = Yes)',
          not mism, f'mismatching samples: {mism or "none"}; samples using each feature: {used.sum().to_dict()}')
log('  feature combinations: ' + '; '.join(f'{k}: {v}' for k, v in combo.value_counts().items()))
seq = yes_col('WES') | yes_col('Target_seq')
log(f'  data-entry notes: TP53/NF1 filled without WES or Target_seq: '
    f'{s1.index[(s1["TP53"].notna() | s1["NF1"].notna()) & ~seq].tolist()}; '
    f'WES/Target_seq = Yes but TP53 or NF1 missing: {s1.index[seq & (s1["TP53"].isna() | s1["NF1"].isna())].tolist()}; '
    'RNA-seq / scRNAseq / CODEX are not model inputs')

### ==================== M4 calls ====================
log('\n== M4  calls and Table S1 cells to change')
indep['Methy'] = pd.concat([methy_disc, methy_val])
prob = indep[FEATURES].mean(axis=1, skipna=True)
call = pd.Series(np.where(prob.isna(), 'NA', np.where(prob > 0.5, 'TL', 'PL')), index=s1.index)
d_prob = float(np.nanmax(np.abs(prob - res['Prob'])))
milestone('M4.1 Prob (mean of available Age/TP53/NF1/Ki67/Methy) and PL/TL call recomputed == results',
          d_prob < 1e-9 and (call == res['Prob_bi_str']).all() and (prob.isna() == res['Prob'].isna()).all(),
          f'Prob max|diff| {d_prob:.1e}; calls identical {int((call == res["Prob_bi_str"]).sum())}/86')

is_disc = s1['Discovery/Validation'] == 'Discovery'
train_label = s1.loc[is_disc, 'Group/Predict Group'].astype(str).str.strip()
proposed = pd.Series(np.where(call == 'NA', 'NA', call + '-SCG classifier'), index=s1.index)
proposed[is_disc] = train_label.reindex(s1.index[is_disc]).values
current = s1['Group/Predict Group'].fillna('NA').astype(str).str.strip()
change = s1.index[current != proposed]
upd = pd.read_csv(BASE / 'results' / 'table_s1_group_update.csv', index_col=0, keep_default_na=False, na_values=[''])
same_list = change.tolist() == upd.index.tolist() and (proposed[change].values == upd['proposed'].values).all()
milestone('M4.2 Table S1 "Group/Predict Group" cells to change', same_list,
          f'{len(change)} cells; identical to results/table_s1_group_update.csv: {same_list}')
log(pd.DataFrame({'Age': s1['Age'], 'Sex': s1['Sex'], 'NF1': s1['NF1'], 'TP53': s1['TP53'], 'features': combo,
                  'Prob': prob.round(3), 'Table S1 now': current, 'change to': proposed}).loc[change].to_string())

s1_cols = ['Discovery/Validation', 'Resource', 'Age', 'Sex', 'NF1', 'TP53', 'Ki67']
same_vals = all(((res[c].isna() & s1[c].isna()) | (res[c].astype(str) == s1[c].astype(str))).all() for c in s1_cols)
avail = res['Methylation availability'].astype(str)
avail_ok = ((avail == s1['DNA methyaltion'].astype(str).str.strip()).all()
            and ((avail.str.lower() == 'yes') == res['Methy_prob'].notna()).all())
grp_ok = (res['Group/Predict Group (Table S1)'] == current).all() and (res['formatted cls'] == proposed).all()
milestone('M4.3 results table carries the Table S1 values it used',
          same_vals and avail_ok and grp_ok,
          f'{", ".join(s1_cols)} identical to Table S1 for 86/86: {same_vals}; "Methylation availability" == Table S1 '
          f'"DNA methyaltion" and == Methy_prob present: {avail_ok} (Yes {int((avail.str.lower() == "yes").sum())}, '
          f'No {int((avail.str.lower() == "no").sum())}); Table S1 group and formatted cls columns correct: {grp_ok}')

log(f'  discovery: model call differs from the training label for {s1.index[is_disc & (call != train_label.reindex(s1.index))].tolist() or "none"}')
border = s1.index[~is_disc & (prob - 0.5).abs().lt(0.05)]
log('  validation calls within 0.05 of 0.5: ' + ', '.join(f'{i} {prob[i]:.3f} ({combo[i]})' for i in border))
p79 = indep.loc['SP79', ['Age', 'NF1', 'Ki67', 'Methy']].mean()
log(f'  SP79 without its unsourced TP53: Prob {p79:.3f} -> {"TL" if p79 > 0.5 else "PL"} (with TP53: {prob["SP79"]:.3f} -> {call["SP79"]})')

log('\n== summary')
for k, v in status.items():
    log(f'  {"PASS" if v else "FAIL"}  {k}')
REPORT.write_text('\n'.join(out) + '\n', encoding='utf-8')
print(f'report -> {REPORT}')
sys.exit(0 if all(status.values()) else 1)
