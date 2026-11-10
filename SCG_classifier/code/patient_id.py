"""Patient-identity checks shared by the scripts in code/.

SP IDs are not stable: between the 2025-06 table (Table S1_0606) and the
2026-09 Table S1, four newly collected cases were inserted as SP70-SP73 and
the 13 GSE270993 cases moved from SP70-SP82 to SP74-SP86. The methylation
matrices were linked to patients in 2025-06, so every SP-ID-keyed join is
checked against a patient fingerprint

    Age, Sex, NF1, TP53

recorded in data/sample_manifest.csv (taken from Table S1_0606 at linking
time) and compared with the current Table S1 row of the same SP ID.

sample_manifest.csv lists every Table S1 sample, so it can be joined against
anything keyed by SP ID. The 13 samples with no methylation data have only
SP_ID and cohort filled; every column describing methylation data, the 2025-06
linking ID and the fingerprint are empty for them. A sample therefore "has
methylation data" iff its 'batch' is set, and that is what every consumer
filters on.

The fingerprint is a check, not a key: several Table S1 patients share one
(e.g. SP41/SP44, SP76/SP81), so a swap inside such a group would pass;
check_manifest() prints those groups on every run. The second guard is that
the set of SP IDs with methylation data must equal Table S1's
'DNA methyaltion' == 'Yes' -- the check that exposes a shifted block of IDs.
"""
import pandas as pd

FINGERPRINT = ['Age', 'Sex', 'NF1', 'TP53']

### Clinical values corrected in the current table after the manifest was built.
### These differ from the 2025-06 fingerprint on purpose, so they must not be
### read as "this SP ID is a different patient now".
###   SP36 NF1 WT -> Mut   (Table S1, 2026-09-16)
CORRECTIONS = {('SP36', 'NF1')}
N_DISCOVERY = 36
MANIFEST_DTYPES = {c: str for c in ('SP_ID', 'cohort', 'batch', 'beta_column', 'raw_column', 'Original_ID',
                                    'methy_ID', 'GSM', 'SP_ID_2025_06', 'Sex', 'NF1', 'TP53',
                                    'Group/Predict Group (Table S1)', 'subtype')}


def _norm(field, v):
    if pd.isna(v):
        return 'NA'
    if field == 'Age':
        return f'{float(v):g}'
    s = str(v).strip()
    if field == 'Sex':
        return s.lower()
    return {'WT': 'WT', 'MUT': 'Mut'}.get(s.upper(), s)


def fingerprint(df):
    """Normalized Age/Sex/NF1/TP53 strings for every row of df (missing -> 'NA')."""
    return pd.DataFrame({f: [_norm(f, v) for v in df[f]] for f in FINGERPRINT}, index=df.index)


def load_table_s1(path):
    """Read the clinical table sorted by SP ID; fail if the ID layout the pipeline relies on has changed.

    The file is kept as a verbatim copy of the collaborator's workbook, so the
    header is located rather than assumed: the published Table S1 puts its
    title in the first row and the header in the second.
    """
    head = pd.read_excel(path, header=None, nrows=2)
    header_row = next((i for i in (0, 1) if str(head.iloc[i, 0]).strip() == 'Sample ID'), None)
    if header_row is None:
        raise ValueError(f'{path}: no "Sample ID" header found in the first two rows')
    table = pd.read_excel(path, index_col=0, header=header_row)
    table.index = table.index.astype(str).str.strip()
    if table.index.duplicated().any():
        raise ValueError(f'Table S1: duplicated sample IDs {table.index[table.index.duplicated()].tolist()}')
    table = table.sort_index()
    expected = [f'SP{i:02d}' for i in range(1, len(table) + 1)]
    if table.index.tolist() != expected:
        raise ValueError(f'Table S1: sample IDs are not exactly SP01..SP{len(table):02d} -- rows were added, removed '
                         'or renumbered; rebuild data/sample_manifest.csv with 00b_build_sample_manifest.py and '
                         'review its renumbering report first')
    missing = [c for c in FINGERPRINT + ['Ki67', 'Resource', 'DNA methyaltion', 'Discovery/Validation',
                                         'Group/Predict Group'] if c not in table.columns]
    if missing:
        raise ValueError(f'Table S1: missing columns {missing}')
    discovery = table.index[table['Discovery/Validation'] == 'Discovery'].tolist()
    if discovery != expected[:N_DISCOVERY]:
        raise ValueError(f'Table S1: discovery cohort is no longer exactly SP01-SP{N_DISCOVERY}: {discovery}')
    return table


def load_manifest(path):
    return pd.read_csv(path, dtype=MANIFEST_DTYPES)


def check_manifest(manifest, table):
    """Raise unless every methylation sample in the manifest matches its Table S1 row on the fingerprint,
    and the SP IDs carrying methylation data are exactly Table S1's 'DNA methyaltion' == 'Yes'.

    The manifest lists every Table S1 sample, including those with no methylation data, so a sample
    "has methylation data" iff its 'batch' is set -- never merely because it has a row. The rows
    without data carry no fingerprint either (it is recorded when a methylation matrix is linked to a
    patient, and nothing was linked for them), so they are left out of both checks below.
    """
    m = manifest.set_index('SP_ID')
    if m.index.duplicated().any():
        raise ValueError(f'sample_manifest: SP_ID used twice: {m.index[m.index.duplicated()].tolist()}')
    unknown = sorted(set(m.index) - set(table.index))
    if unknown:
        raise ValueError(f'sample_manifest: SP IDs not in Table S1: {unknown}')

    m = m[m['batch'].notna()]
    recorded, current = fingerprint(m), fingerprint(table.loc[m.index])
    differ = [(sp, f, recorded.loc[sp, f], current.loc[sp, f])
              for sp in m.index for f in FINGERPRINT if recorded.loc[sp, f] != current.loc[sp, f]]
    bad = [d for d in differ if (d[0], d[1]) not in CORRECTIONS]
    for sp, f, was, now in differ:
        if (sp, f) in CORRECTIONS:
            print(f'  known correction allowed: {sp} {f} {was} -> {now}')
    if bad:
        raise ValueError('patient fingerprint mismatch -- the methylation sample recorded under this SP ID is not '
                         'the Table S1 patient with that ID (or Table S1 was edited):\n'
                         + pd.DataFrame(bad, columns=['SP_ID', 'field', 'manifest', 'Table S1']).to_string(index=False))

    with_data = set(m.index)
    says_yes = set(table.index[table['DNA methyaltion'].astype(str).str.strip().str.lower() == 'yes'])
    if with_data != says_yes:
        raise ValueError('methylation data vs Table S1 "DNA methyaltion": data but not "Yes" for '
                         f'{sorted(with_data - says_yes)}; "Yes" but no data for {sorted(says_yes - with_data)}')

    key = fingerprint(table).agg('|'.join, axis=1)
    shared = key[key.duplicated(keep=False)]
    groups = sorted('/'.join(ids) for ids in shared.groupby(shared).groups.values() if set(ids) & with_data)
    print(f'patient check OK: {len(m)} methylation samples match Table S1 on {"/".join(FINGERPRINT)}; '
          'samples with data == Table S1 "DNA methyaltion = Yes"')
    if groups:
        print(f'  fingerprint shared within {", ".join(groups)} -- a swap inside one of these would not be detected')


def check_batch_columns(columns, manifest, batch):
    expected = set(manifest.loc[manifest['batch'] == batch, 'beta_column'])
    got = set(map(str, columns))
    if got != expected:
        raise ValueError(f'{batch}: columns differ from sample_manifest.csv -- extra {sorted(got - expected)}, '
                         f'missing {sorted(expected - got)}')
