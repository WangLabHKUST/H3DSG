# -*- coding: utf-8 -*-
"""Build ../../sample_metadata.tsv -- the single sample manifest of this repository -- and
../../raw_data/sample_key.tsv, the table that ties raw_data/ and data_full/ to it.

Input (all under ../../raw_data/common_data/; archive version only, because they hold the
hospital record numbers, which never leave raw_data/):
  Table S1_0921.xlsx                        Table S1 of the paper, **the authoritative
                                            source of the PL/TL call**
  SP_ID_pathology_ID_mapping.xlsx           record number <-> SP number <-> original key
  samples.csv                               the list of 35 CODEX samples, by original key
  g1_samples/g2_samples in ../../constants.py
                                            the grouping that actually produced the
                                            existing outputs (frozen)

Output
  ../../sample_metadata.tsv   35 rows, published. Columns:

    CODEX_ID          the SP ID: the key of adata_dict and what every file name and
                      variable name of this repository uses
    Original_ID       the SP ID again. The record number the hospital assigned is not
                      published; it is read from the mapping table here and goes no
                      further than the key table below
    SP_ID             the Sample ID of Table S1
    group             PL / TL, taken from Group/Predict Group of **Table S1**
    group_as_run      PL / TL, taken from g1/g2 in constants.py, i.e. the grouping that
                      produced the existing outputs
    in_analysis       whether the sample enters the analysis (31 True, 4 False)
    exclusion_reason  low_quality / recurrent / empty

  ../../raw_data/sample_key.tsv   35 rows, columns raw_key and SP_ID, in ascending order
                      of raw_key. raw_data/ and data_full/ are keyed by the original
                      sample key; data_loading.py reads this table to translate the keys
                      to SP IDs, and its row order is the order the samples had when the
                      tables in data/ were first written (data_loading.legacy_sample_order).
                      Not published.

Three things to know:

1. Two original keys do not match the record number in the mapping table (a digit was
   mistyped when the key was made), so the record number cannot be read off the key. This
   script always joins on the key column of the mapping table and then takes the Sample ID
   from it; it reports how many keys differ, never which.

2. **PL/TL always comes from Table S1; do not use the Group column of the mapping
   table.** Measured, that column agrees with Table S1 for 34 of the 35 CODEX samples,
   but records SP41 as PL, whereas Table S1 has `TL-SCG classifier`, the SCG classifier
   gives Prob = 0.5173 > 0.5, and samples.csv on the CODEX side also puts it in Group2
   (TL). In other words that one cell of the mapping table is stale. The script asserts
   that the two agree below, so if this ever happens again it fails loudly instead of
   quietly carrying on.

3. Four samples are in samples.csv but do not enter the analysis; the exclusion happens
   two steps upstream:
       SP02, SP50   low_quality  (sample_lst_33 in codex/utils.py)
       SP42, SP43   recurrent    (sample_lst_31 in codex/utils.py)
"""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from constants import g1_samples, g2_samples

COMMON_DATA = REPO / 'raw_data' / 'common_data'
TABLE_S1 = COMMON_DATA / 'Table S1_0921.xlsx'
OUT_PATH = REPO / 'sample_metadata.tsv'
KEY_PATH = REPO / 'raw_data' / 'sample_key.tsv'

### Upstream exclusion: sample_lst_33 / sample_lst_31 in codes/cleaned_202604/codex/utils.py
EXCLUSION = {'SP02': 'low_quality', 'SP50': 'low_quality',
             'SP42': 'recurrent', 'SP43': 'recurrent'}

### The validation cohort in Table S1 carries this suffix; strip it to get a bare PL/TL
SUBTYPE_SUFFIX = '-SCG classifier'
### The Group column of the mapping table writes TL as TML
GROUP_FROM_MAPPING = {'PL': 'PL', 'TML': 'TL'}


def load_table_s1_subtype():
    """Read Table S1 and return SP number -> bare PL/TL.

    Row 1 of Table S1 is the table title and row 2 is the header, so the header row has
    to be located rather than assumed.
    """
    head = pd.read_excel(TABLE_S1, header=None, nrows=2)
    header_row = next((i for i in (0, 1) if str(head.iloc[i, 0]).strip() == 'Sample ID'), None)
    if header_row is None:
        raise ValueError(f'{TABLE_S1}: no "Sample ID" header found in the first two rows')
    t = pd.read_excel(TABLE_S1, index_col=0, header=header_row)
    t.index = t.index.astype(str).str.strip()

    raw = t['Group/Predict Group'].astype(str).str.strip()
    subtype = raw.map(lambda v: v[:-len(SUBTYPE_SUFFIX)] if v.endswith(SUBTYPE_SUFFIX) else v)
    bad = set(subtype) - {'PL', 'TL'}
    if bad:
        raise ValueError(
            f'Group/Predict Group in Table S1 still has non-PL/TL values after stripping '
            f'the suffix: {sorted(bad)}')
    return subtype


def main():
    samples = pd.read_csv(COMMON_DATA / 'samples.csv', dtype=str)
    mapping = pd.read_excel(COMMON_DATA / 'SP_ID_pathology_ID_mapping.xlsx', dtype=str)
    s1_subtype = load_table_s1_subtype()

    mapping = mapping.dropna(subset=['CODEX_ID']).set_index('CODEX_ID')
    if not mapping.index.is_unique:
        raise ValueError(
            f'duplicate CODEX_ID in the mapping table: '
            f'{mapping.index[mapping.index.duplicated()].tolist()}')

    rows, key_rows, mapping_conflicts, key_differs = [], [], [], 0
    for raw_key in samples['sample']:
        digits = raw_key.lstrip('X')
        if digits not in mapping.index:
            raise ValueError('a sample of samples.csv is not in SP_ID_pathology_ID_mapping.xlsx: '
                             f'{raw_key!r}')
        hit = mapping.loc[digits]
        sp_id = hit['ID_Cell']
        if sp_id not in s1_subtype.index:
            raise ValueError(f'{sp_id} is not in Table S1')
        key_differs += int(hit['Original_ID'] != digits)

        group = s1_subtype[sp_id]                       # authoritative: Table S1
        from_mapping = GROUP_FROM_MAPPING.get(str(hit['Group']).strip())
        if from_mapping != group:
            mapping_conflicts.append((sp_id, from_mapping, group))

        if sp_id in g1_samples:
            as_run = 'PL'
        elif sp_id in g2_samples:
            as_run = 'TL'
        else:
            as_run = ''                                 # the 4 excluded ones have no as-run group

        key_rows.append({'raw_key': raw_key, 'SP_ID': sp_id})
        rows.append({
            'CODEX_ID': sp_id,
            'Original_ID': sp_id,
            'SP_ID': sp_id,
            'group': group,
            'group_as_run': as_run,
            'in_analysis': sp_id not in EXCLUSION,
            'exclusion_reason': EXCLUSION.get(sp_id, ''),
        })

    out = pd.DataFrame(rows).sort_values(
        'SP_ID', key=lambda s: s.str.extract(r'(\d+)', expand=False).astype(int))
    key = pd.DataFrame(key_rows).sort_values('raw_key').reset_index(drop=True)

    ### ---------------- self-checks ----------------
    assert len(out) == 35, f'expected 35 rows, got {len(out)}'
    assert out['in_analysis'].sum() == 31, \
        f'expected 31 samples in the analysis, got {out["in_analysis"].sum()}'
    for col in ('SP_ID',):
        dup = out[col][out[col].duplicated()]
        assert dup.empty, f'duplicate {col}: {dup.tolist()}'
    assert key['raw_key'].is_unique and key['SP_ID'].is_unique, 'sample_key.tsv must be one to one'
    assert (out.loc[~out['in_analysis'], 'group_as_run'] == '').all(), \
        'an excluded sample must not have an as-run group'
    assert (out.loc[out['in_analysis'], 'group_as_run'] != '').all(), \
        'every sample in the analysis must have an as-run group'
    assert set(out.loc[out['in_analysis'] & (out['group_as_run'] == 'PL'), 'CODEX_ID']) == set(g1_samples)
    assert set(out.loc[out['in_analysis'] & (out['group_as_run'] == 'TL'), 'CODEX_ID']) == set(g2_samples)

    ### The existing outputs were produced under group_as_run; it must agree with
    ### Table S1, or the results do not match the paper
    inc = out[out['in_analysis']]
    drift = inc[inc['group'] != inc['group_as_run']]
    assert drift.empty, ('the grouping of the existing outputs disagrees with Table S1, '
                         'so the results need to be rerun:\n'
                         + drift[['SP_ID', 'group', 'group_as_run']].to_string(index=False))

    out.to_csv(OUT_PATH, sep='\t', index=False)
    key.to_csv(KEY_PATH, sep='\t', index=False, lineterminator='\n')

    ### ---------------- report ----------------
    print(f'{len(out)} samples, {len(inc)} in the analysis, {len(out) - len(inc)} excluded')
    print(f"  group (Table S1, authoritative) : {inc['group'].value_counts().to_dict()}")
    print(f"  group_as_run (as the outputs ran): {inc['group_as_run'].value_counts().to_dict()}")
    print('  the two agree -- the existing results need no rerun')

    if mapping_conflicts:
        print(f'\n  WARNING: the Group column of SP_ID_pathology_ID_mapping.xlsx disagrees '
              f'with Table S1 on {len(mapping_conflicts)} sample(s) (Table S1 was used):')
        for sp_id, got, want in mapping_conflicts:
            print(f'      {sp_id}: mapping table={got}, Table S1={want}')
    else:
        print('  the Group column of the mapping table agrees with Table S1 throughout')

    print(f'\n  {key_differs} original key(s) differ from the record number in the mapping table')

    for reason, grp in out[~out['in_analysis']].groupby('exclusion_reason'):
        print(f"  excluded ({reason}): " + ', '.join(grp['SP_ID']))

    print(f'\nsaved -> {OUT_PATH}')
    print(f'saved -> {KEY_PATH}')


if __name__ == '__main__':
    main()
