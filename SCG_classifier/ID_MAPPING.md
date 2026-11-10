# Identifiers and how to convert between them

A patient in this study appears under up to eight different identifiers. This file
describes each one, which files are authoritative for it, and how to get from a
methylation matrix column back to the patient. Every count below was recomputed from
the files in this folder.

## The key

**The primary key is the current SP ID**, the `Sample ID` column of Table S1, taking
the values `SP01`–`SP86`. Everything in `results/` is reported under it. Every other
identifier is an alias, to be converted by lookup and **not by string manipulation** —
the two exceptions are called out explicitly below.

The cohort split is fixed: `SP01`–`SP36` are the discovery cohort (36 patients),
`SP37`–`SP86` the validation cohort (50). `patient_id.load_table_s1()` enforces this,
so "the discovery cohort is the first 36 rows" is a safe assumption.

73 patients have methylation data (`DNA methyaltion = Yes` in Table S1); 13 do not
(SP45, SP47, SP61, SP63–SP65, SP67–SP69, SP70–SP73).

## Identifier inventory

| Identifier | Example | Covers | Authoritative source | Derivable? |
|---|---|---|---|---|
| `Sample ID` (current SP ID) | `SP01` … `SP86` | all 86 | `data/Table S1.xlsx` | — (this is the key) |
| `SP_ID_2025_06` | `SP01` … `SP82` | 73 (those with methylation data) | `data/id_mapping/Table S1_0606.xlsx`; written to `data/sample_manifest.csv` | **yes**, see "Renumbering" |
| `ID_in_network` | `P01` … `P36` | discovery only | `data/ID_relationship_36cases_NF1_TP53.csv` | **yes**: current SP ID = `'S' + ID_in_network` |
| `Multiomics_analysis_ID` | `SPH3_37` | discovery only | same | **no**, must be looked up |
| `PID` | `NNNNNN` | discovery only | same | **yes**: `'P' + Original_ID` |
| `Original_ID` | a 6-digit hospital number, or `S01`–`S25` for GSE270993 | 73 | `ID_relationship…csv` (discovery), `data/id_mapping/Methy_Cell_ID_map_to_yy.xlsx` (validation) | **no** |
| `methy_ID` (Sentrix) | `203789410095_R05C01` | 37 validation | `Methy_Cell_ID_map_to_yy.xlsx` | **no** |
| `GSM` | `GSM8366794` | the 13 GSE270993 patients | `data/id_mapping/GSE270993_sample_info.csv` | **no** |

Each identifier is one-to-one within its coverage; `00d_build_id_summary.py` asserts
uniqueness for `Original_ID`, `methy_ID`, `GSM`, `SP_ID_2025_06` and `PID`.

### `Multiomics_analysis_ID` must be looked up

It has the form `SPH3_<n>` with 36 values, `<n>` is not contiguous (it runs up to
`SPH3_45`), and **the order does not follow the P-numbers**: `P01 = SPH3_37`,
`P06 = SPH3_7`, `P36 = SPH3_34`. Any assumption that the two run in parallel is wrong.

### `Original_ID` spans two namespaces

Newly collected patients (60 of the 73) have a 6-digit hospital number. The 13
GSE270993 patients instead carry the GEO-side sample names `S01`–`S25`, which are not
contiguous. **`S01` is not `SP01`** — `S01` is `SP74`. Do not treat the two as the
same namespace, and do not sort them together.

A further 8 hospital numbers exist in `Methy_Cell_ID_map_to_yy.xlsx` for patients
without methylation data (SP45, SP47, SP61, SP63–SP65, SP67, SP68), so they are absent
from `sample_manifest.csv` and from `results/id_summary_all_samples.csv`, both of
which cover only the 73.

## From a methylation matrix column back to the patient

Four chains, one per `batch` value in `data/sample_manifest.csv`:

| `batch` | n | Column name in the raw matrix (`raw_column`) | Via | Gives |
|---|---|---|---|---|
| `discovery` | 36 | `Original_ID` (hospital number) | `ID_relationship…csv` → `ID_in_network` | current SP ID = `'S' + Pxx` |
| `beta_55` | 19 | `Original_ID` (hospital number) | `Methy_Cell_ID_map_to_yy.xlsx`, `Cell ID` column | 2025-06 SP ID |
| `beta_5` | 5 | Sentrix, e.g. `208776820119_R05C01` | same, matched on `methy_ID` | 2025-06 SP ID |
| `beta_14` | 13 | `GSM` | `GSE270993_sample_info.csv` gives `Sentrix_ID + '_' + Sentrix_Position`, then the table above | 2025-06 SP ID |

The last step, 2025-06 SP ID → current SP ID, is below. In practice only the 13
`beta_14` patients changed number; discovery, `beta_55` and `beta_5` did not.

### What the delivered matrices are actually keyed by

Easy to get wrong, so note the difference:

| File | Columns are |
|---|---|
| `data/beta.csv` (discovery, 141824 CpGs × 36) | **`Original_ID`**, e.g. `NNNNNN` |
| `data/val_data/beta_5_filtered.csv` (5) | **current SP ID**: SP46, SP48, SP60, SP62, SP66 |
| `data/val_data/beta_14_filtered.csv` (13) | **current SP ID**: SP74–SP86 |
| `data/val_data/beta_55_filtered.csv` (19) | **current SP ID**: SP37–SP44, SP49–SP59 |

So the validation matrices are already converted and can be used by SP ID directly;
`beta.csv` is not and has to go through `ID_relationship`. In `sample_manifest.csv`,
`beta_column` is the column name in the *delivered* matrix and `raw_column` the name
in the raw matrix.

The raw validation matrices are not delivered (too large). Their column headers are
kept in `data/id_mapping/raw_val_beta_columns.csv`: raw `beta_55` has 55 columns, 36
of which are the discovery hospital numbers and 19 the validation patients; raw
`beta_14` has 28 columns of which only 13 belong to this study; raw `beta_5` has 5,
all used.

## Renumbering: 2025-06 SP ID → current SP ID

The methylation data was linked to patients in 2025-06 against `Table S1_0606.xlsx`
(82 rows). The current Table S1 has 86: **four newly collected validation patients
were inserted as `SP70`–`SP73`** — none of them has methylation data, so no linkage is
affected — and the 13 GSE270993 patients that followed them shifted up by four.

```
current = 2025-06                 for SP01 .. SP69
current = 2025-06 + 4             for 2025-06 SP70 .. SP82
SP70 .. SP73                      new in the current table, absent in 2025-06
```

`00b_build_sample_manifest.py` derives this by matching the old table onto the current
one order-preservingly on Age/Sex/NF1/TP53/Resource/cohort, accepting the result only
when the leftmost and rightmost order-preserving matches coincide — that is, when the
match is unique, so an inserted row can never be taken for an old one even where
fingerprints repeat.

> **A blank `SP_ID_2025_06` does not mean "newly inserted".** The column is filled only
> for the 73 patients with methylation data, so its 13 blanks are the 13 patients
> without data. Of those, only SP70–SP73 are actually new; the other nine (SP45, SP47,
> SP61, SP63–SP65, SP67–SP69) existed in 2025-06 under the same number.

## The machine-readable tables

Two files carry the mappings themselves:

- **`data/sample_manifest.csv`** — 86 rows, one per Table S1 patient, 15 columns. The
  join entry point for anything keyed by SP ID. The 13 patients without methylation
  data have only `SP_ID`, `cohort` and the two subtype columns filled; every column
  describing methylation data, the 2025-06 linking ID and the fingerprint are empty for
  them. **A patient "has methylation data" iff `batch` is set.** Not in the repository
  (see `.gitignore`); rebuild with `00b_build_sample_manifest.py`.
- **`results/id_summary_all_samples.csv`** — 86 rows, 18 columns: every identifier side
  by side, plus cohort, resource, `DNA methyaltion` and the fingerprint. Not in the
  repository either, for the same reason; rebuild with `00d_build_id_summary.py`.

Blank always means "this identifier does not exist for this patient", nothing else.

The last two columns of `sample_manifest.csv` carry the subtype, both read off
Table S1 rather than off `results/` — `02_predict_all_samples.py` reads the manifest,
so reading its own output back in would make an input depend on an output:

| Column | Values |
|---|---|
| `Group/Predict Group (Table S1)` | `PL`/`TL` for discovery (the training label), `PL-SCG classifier`/`TL-SCG classifier` for validation. Same name and same strings as in `results/subtype_predictions_all_samples.csv`, so the two join directly |
| `subtype` | the same thing as plain `PL`/`TL`. 48 TL / 38 PL over all 86 |

`Prob` and the per-feature probabilities live in the results table, not in the
manifest.

## Rules that are easy to break

1. **When joining on SP ID, carry and check the fingerprint `Age`/`Sex`/`NF1`/`TP53`.**
   SP IDs have been renumbered, so matching the number alone cannot detect a shift.
   `patient_id.check_manifest()` does this on every run of 00/00b/00c/00d/01/02/03.
2. **The fingerprint is a check, not a key.** Twelve groups of patients share one, and
   a swap inside a group would pass: SP01/SP12/SP45, SP06/SP70, SP08/SP11, SP24/SP31,
   SP26/SP59, SP27/SP57, SP29/SP68, SP33/SP61, SP41/SP44, SP48/SP75, SP62/SP73,
   SP76/SP81. What uniquely identifies a methylation sample is `Original_ID`,
   `methy_ID` or `GSM`. `check_manifest()` prints these groups every run.
3. **`ID_relationship_36cases_NF1_TP53.csv` is for ID mapping only.** Its `Group`,
   `NF1` and `TP53` columns are hg19-era calls — SP36 is recorded there as `CNV_Loss` —
   and must not enter anything. Every script restricts the read with `usecols`. Labels
   and mutation calls come from Table S1's `Group/Predict Group`, `NF1` and `TP53`.
4. **Clinical values come from Table S1 only.** No other copy of Age, Ki-67, NF1 or
   TP53 is authoritative.
5. **Locate Table S1's header row; do not assume it.** Row 1 is the table title and
   row 2 is the header. `patient_id.load_table_s1()` looks for `Sample ID` in the first
   two rows.
6. **Copy column spellings exactly**: `DNA methyaltion` (the typo is in the published
   table), `Group/Predict Group`, and `Ki67` — which the 2025-06 table called `MKI67`.
7. **The subtype column in the results tables is `Group/Predict Group (Table S1)`**,
   not `Group/Predict Group`.
