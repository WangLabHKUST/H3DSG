"""Write results/Table S1 - updated.xlsx: a copy of Table S1 with every
'Group/Predict Group' cell listed in results/table_s1_group_update.csv (from
02_predict_all_samples.py) replaced by its proposed value. Nothing else is
touched -- the workbook is edited cell by cell with openpyxl, so other values,
formatting and row order stay as in the collaborator's file.

Before a cell is changed, the xlsx row found by Sample ID must match the
update row on the patient fingerprint Age/Sex/NF1/TP53 and on the current
'Group/Predict Group' value. After saving, the file is read back and checked:
only the listed cells differ from the source, and 'Group/Predict Group' equals
'formatted cls' of results/subtype_predictions_all_samples.csv for all samples.

usage: python code/02b_write_table_s1_updated.py [SOURCE_XLSX]   (default data/Table S1.xlsx)
"""
import hashlib
import sys
from pathlib import Path

import openpyxl
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from patient_id import FINGERPRINT, fingerprint, load_table_s1

BASE = Path(__file__).resolve().parent.parent
source_path = Path(sys.argv[1]) if len(sys.argv) > 1 else BASE / 'data' / 'Table S1.xlsx'
update_path = BASE / 'results' / 'table_s1_group_update.csv'
pred_path = BASE / 'results' / 'subtype_predictions_all_samples.csv'
out_path = BASE / 'results' / 'Table S1 - updated.xlsx'
if out_path.exists():
    raise FileExistsError(f'{out_path} already exists, not overwriting')

GROUP = 'Group/Predict Group'
print(f'source: {source_path} (md5 {hashlib.md5(source_path.read_bytes()).hexdigest()})')

update = pd.read_csv(update_path, index_col=0)
wb = openpyxl.load_workbook(source_path)
if len(wb.worksheets) != 1:
    raise ValueError(f'expected one sheet, found {wb.sheetnames}')
ws = wb.active
### the collaborator's workbook may carry its title in the first row (Table S1 does)
header_row = next((i for i in (1, 2) if str(ws.cell(i, 1).value).strip() == 'Sample ID'), None)
if header_row is None:
    raise ValueError(f'{source_path}: no "Sample ID" header found in the first two rows')
header = [c.value for c in ws[header_row]]
col = {name: i for i, name in enumerate(header)}
rows = {str(r[0].value).strip(): r for r in ws.iter_rows(min_row=header_row + 1) if r[0].value is not None}

sheet_values = pd.DataFrame({f: [rows[sp][col[f]].value for sp in update.index] for f in FINGERPRINT},
                            index=update.index)
fp_sheet, fp_update = fingerprint(sheet_values), fingerprint(update)
changes = []
for sp, u in update.iterrows():
    cell = rows[sp][col[GROUP]]
    if not fp_sheet.loc[sp].equals(fp_update.loc[sp]):
        raise ValueError(f'{sp}: fingerprint in {source_path.name} {fp_sheet.loc[sp].tolist()} != '
                         f'update list {fp_update.loc[sp].tolist()}')
    if str(cell.value).strip() != u['Table S1 (current)']:
        raise ValueError(f'{sp}: {GROUP} is {cell.value!r} in {source_path.name}, '
                         f'update list expects {u["Table S1 (current)"]!r}')
    changes.append((sp, cell.coordinate, cell.value, u['proposed']))
    cell.value = u['proposed']
wb.save(out_path)

### ---- read back and check ----
src = list(openpyxl.load_workbook(source_path).active.iter_rows(values_only=True))
dst = list(openpyxl.load_workbook(out_path).active.iter_rows(values_only=True))
assert len(src) == len(dst) and all(len(a) == len(b) for a, b in zip(src, dst))
differ = {(str(a[0]).strip(), header[j], a[j], b[j])
          for a, b in zip(src, dst) for j in range(len(a)) if a[j] != b[j]}
assert differ == {(sp, GROUP, old, new) for sp, _, old, new in changes}, differ

table = load_table_s1(out_path)
pred = pd.read_csv(pred_path, index_col=0, keep_default_na=False)
group = table[GROUP].fillna('NA').astype(str).str.strip()
mismatch = group.index[group != pred.loc[group.index, 'formatted cls']].tolist()
assert not mismatch, f'{GROUP} still differs from formatted cls for {mismatch}'

print(f'{len(changes)} cells changed:')
for sp, coord, old, new in changes:
    print(f'  {sp} ({coord}, {"/".join(fp_update.loc[sp])}): {old} -> {new}')
print(f'read back: only these cells differ from the source; {GROUP} == formatted cls for all {len(group)} samples')
print(f'saved -> {out_path}')
