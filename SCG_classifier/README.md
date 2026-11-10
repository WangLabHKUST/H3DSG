# SCG classifier

Naive Bayes classifier that assigns a spinal cord glioma to the PL (H3-DSG-PL) or
TL (H3-DSG-TL) subtype from clinical, mutational and DNA methylation features.

The model is trained on the 36 discovery patients and applied to the 50 validation
patients. Each feature is scored by its own model and the per-feature probabilities
of being TL are averaged, skipping features a patient does not have:

| Feature | Model | Source |
|---|---|---|
| Age | GaussianNB | Table S1 |
| Ki-67 | GaussianNB | Table S1 |
| TP53 | BernoulliNB (WT = 0, Mut = 1) | Table S1 |
| NF1 | BernoulliNB (WT = 0, Mut = 1) | Table S1 |
| Methy | GaussianNB per CpG, averaged over the top 10000 CpGs | 450K/EPIC beta values |

`Prob` is the mean of the available probabilities; `Prob > 0.5` is called TL. A
patient with none of the five features is left NA. Discovery patients get their
methylation probability from leave-one-out fits, so they are never scored by a model
fitted on themselves.

Current result: 48 TL / 38 PL over all 86 patients — discovery 20 TL / 16 PL
(the training labels), validation 28 TL / 22 PL. The call agrees with Table S1's
`Group/Predict Group` for all 86 patients, which is why `table_s1_group_update.csv`
is empty.

## Layout

```
code/        11 scripts, see below
data/        all inputs (not in the repository, see "Data")
models/      the three fitted model files (not in the repository)
results/     the eight output tables
```

Every script resolves its own paths from its location — `BASE` is the parent of
`code/` — so the tree runs as it stands, with one exception documented under
`04_independent_check.py` below.

## Scripts

Run in this order. Steps 00, 00b and 00c only need to be run to rebuild an input
that already ships in `data/`; the pipeline proper is 01 → 04.

| Script | What it does | Run time |
|---|---|---|
| `00_rank_cpgs_by_loov_margin.py` | ranks all 141824 CpGs by their single-CpG leave-one-out margin; the top 10000 become the methylation feature set. Writes `data/beta_score_full.csv` | ~1 h |
| `00b_build_sample_manifest.py` | builds `data/sample_manifest.csv`, one row per Table S1 patient, linking each methylation column to a current SP ID | seconds |
| `00c_build_validation_beta.py` | extracts the top 10000 CpGs from the three raw validation batches and renames the columns to SP IDs | minutes |
| `00d_build_id_summary.py` | writes `results/id_summary_all_samples.csv`, every identifier each patient is known by, side by side | seconds |
| `01_train_model.py` | fits the four clinical models, the per-CpG models and the discovery leave-one-out methylation values | ~10 min |
| `02_predict_all_samples.py` | scores all 86 patients; writes the prediction tables and the list of Table S1 cells to change | seconds |
| `02b_write_table_s1_updated.py` | writes `results/Table S1 - updated.xlsx`, a copy of Table S1 with those cells applied | seconds |
| `02c_summarize_group_counts.py` | PL/TL counts per cohort, model call against Table S1 | seconds |
| `03_loov_performance_by_feature_combo.py` | leave-one-out accuracy, F1, sensitivity and specificity for 15 feature combinations | ~1 min |
| `04_independent_check.py` | re-derives every number above in closed form, without scikit-learn, and reports 15 pass/fail milestones | ~1 h |
| `patient_id.py` | identity checks shared by the scripts above | — |

Scripts refuse to overwrite an output that already exists; delete it first to
rerun.

### `04_independent_check.py` needs the original source tree

This is the one script that does not run from this folder alone, deliberately: it
checks the release against the *original* sources rather than against the copies
in `data/`. Two paths point outside, both overridable on the command line:

```
argv[1]  the collaborator's Table S1 workbook   (milestone M0.1 compares it cell by cell)
argv[2]  the directory holding the raw, undelivered beta_{5,14,55}.csv
```

Without them it runs as far as milestone M1.1 and then stops on a missing
`beta_55.csv`. The report it produced on the reference machine ships in
`results/independent_check_report.txt`: 15 milestones, all PASS.

## Data

Unlike the earlier version of this module, `data/` and `models/` are present on
disk — the folder is self-contained — but they are listed in `.gitignore` and are
not part of the repository. `data/` holds hospital numbers, Sentrix array IDs and
the raw beta matrices; `models/` holds binary artifacts that `01_train_model.py`
reproduces exactly from `data/`.

```
data/Table S1.xlsx                     the clinical table, 86 patients; the single
                                       source of Age, Ki-67, NF1, TP53 and the label
data/beta.csv                          discovery beta matrix, 141824 CpGs x 36
data/beta_score_full.csv               the fixed CpG ranking from step 00
data/sample_manifest.csv               methylation column -> current SP ID + fingerprint
data/ID_relationship_36cases_NF1_TP53.csv   Original_ID <-> P-number, ID mapping only
data/id_mapping/                       the four tables the manifest is built from
data/val_data/beta_{5,14,55}_filtered.csv   validation betas, top 10000 CpGs, SP ID columns
models/model_dict_meta.npy             the four clinical models
models/model_dict_cpg_alldata.npy      the 10000 per-CpG models
models/training_methy_loov.csv         discovery leave-one-out P(TL)
```

Clinical values are read from Table S1 only — never from any other copy.
`ID_relationship_36cases_NF1_TP53.csv` is used for ID mapping and nothing else:
every script that opens it restricts the read to its ID columns, because its
`Group`, `NF1` and `TP53` columns hold hg19-era calls that disagree with Table S1.

`results/` is published except for `id_summary_all_samples.csv`, which carries
hospital numbers and array IDs; rerun `00d_build_id_summary.py` to regenerate it.

## Sample identity

SP IDs are not stable across table versions: four cases were inserted as SP70-SP73
and the 13 GSE270993 cases moved from SP70-SP82 to SP74-SP86. Every join on SP ID is
therefore checked against a patient fingerprint (Age, Sex, NF1, TP53) recorded in
`data/sample_manifest.csv` when the methylation data was linked, and against the
requirement that the set of patients with methylation data equals the
`DNA methyaltion = Yes` column of Table S1. The fingerprint is a check rather than a
key, since 12 groups of patients share one; `patient_id.check_manifest()` prints
those groups on every run.

`sample_manifest.csv` has one row per Table S1 patient, all 86. The 73 with
methylation data carry the full set of columns; the other 13 carry only `SP_ID`,
`cohort` and the subtype. **A patient "has methylation data" iff `batch` is set** —
never merely because a row exists.

`ID_MAPPING.md` describes every identifier in the study and how to convert between
them.

## One clinical value was corrected after the methylation data was linked

Table S1 records SP36's NF1 as `Mut`; the 2025-06 table recorded `WT`. The
correction follows the hg38 SAVI call in Table S2 — an NF1 frameshift at
chr17:31337558 A>AGAAG, exon 43, 31% VAF in tumour and 0/77 reads in blood.

Two consequences, both handled in the code:

- `patient_id.CORRECTIONS = {('SP36', 'NF1')}` lets `check_manifest()` accept the
  one fingerprint cell that differs, instead of reading it as a patient mismatch.
- `00b_build_sample_manifest.py` matches the 2025-06 table onto Table S1 on the
  fingerprint, so the correction breaks the full-field match. It falls back to a key
  without the corrected field and accepts the fallback only if the match is still
  unique **and** the full-field differences it implies are exactly `CORRECTIONS`.

The correction moves the NF1 model only — P(TL | WT) 0.338 → 0.309, P(TL | Mut)
0.877 → 0.885 — and flips no call among the 86 patients. It does change the NF1
feature's reported performance: leave-one-out accuracy 0.778 → 0.806, F1
0.777 → 0.805, sensitivity 0.65 → 0.70, specificity unchanged at 0.9375; in-sample
AUC 0.79 → 0.82. The methylation models are untouched.

## Environment

Python 3.12, scikit-learn 1.8.0, pandas 3.0.5, numpy 1.26.4, scipy 1.13.1,
openpyxl 3.1.5.
