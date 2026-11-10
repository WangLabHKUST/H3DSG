celltype_order = ['Cycling','OPC_OC_like','AC_MES_like','Myeloid','Neutrophil','T_NK','B_cell','Endothelial_cell','Astrocyte','OPC_OC']

# ---------------------------------------------------------------------------
# g1_samples / g2_samples are the grouping that ACTUALLY produced data/ and figures/.
#
# Keep these two lists frozen: every result already in the repository was computed with
# them. New code should use load_samples() below, which reads sample_metadata.tsv and
# keeps the two notions in separate columns:
#
#   group_as_run  Group1/Group2, equivalent to g1/g2 here -- use this column to
#                 reproduce the existing outputs
#   group         PL/TL, the authoritative grouping from Table S1 of the paper
#
# The two columns currently **agree completely** (PL 14 / TL 17 over the 31 samples), so
# the existing results match the paper. The one sample that ever looked contradictory is
# SP41: the Group column of common_data/SP_ID_pathology_ID_mapping.xlsx records it as PL,
# while Table S1 has TL-SCG classifier, the classifier gives Prob > 0.5, and samples.csv
# on the CODEX side also puts it in Group2. The conclusion is that the cell in that
# spreadsheet is stale and SP41 is TL; sample_metadata.tsv follows Table S1. See
# "Samples and groups" in the top-level README.
# ---------------------------------------------------------------------------
g1_samples = ['SP11', 'SP06', 'SP15', 'SP08', 'SP03', 'SP04', 'SP16', 'SP05', 'SP12', 'SP13', 'SP10', 'SP01', 'SP09',  'SP14']
g2_samples = ['SP31', 'SP35', 'SP28', 'SP24', 'SP18', 'SP17', 'SP19', 'SP21', 'SP22', 'SP32', 'SP30', 'SP23', 'SP34', 'SP25', 'SP41', 'SP26', 'SP36']

color_dict = {'Cycling':'#1C86EE','OPC_OC_like':'#FF7F00','AC_MES_like':'#b62526','Myeloid':'#57b894','Neutrophil':'#bdbd45',
'T_NK':'#519e3e','B_cell':'#3c75af','Endothelial_cell':'#b0aad1','Astrocyte':'#f7b6d2','OPC_OC':'#9edae5'}

# ---------------------------------------------------------------------------
# Naming convention
#
# Development (codes/all, codes/cleaned_202604) used one set of short names; the
# submitted manuscript uses celltype_order above. The difference is a _like suffix on
# the malignant states, which separates them from the normal cell types of the same name:
#
#   old OPC_OC    (malignant state) -> new OPC_OC_like
#   old AC_MES    (malignant state) -> new AC_MES_like
#   old Oligo_pro (oligodendrocyte precursor, a normal cell) -> new OPC_OC
#
# Note that OPC_OC does not refer to the same population in the two sets, so establish
# which set applies before porting code between scripts. The DISPLAY_LABEL of old code
# (e.g. 20260507_cluster_coeff_group_boxplot_ggplot.py:51) shows the malignant OPC_OC as
# "OPC/OC", which is exactly the problem the rename fixes.
# ---------------------------------------------------------------------------

# old name -> new name. Applied once during preprocessing; everything downstream sees
# only the new names.
old2new_celltype = {
    'OPC_OC': 'OPC_OC_like',
    'AC_MES': 'AC_MES_like',
    'Oligo_pro': 'OPC_OC',
}

# new name -> the canonical display name used in figures
celltype_display = {
    'Cycling': 'Cycling',
    'OPC_OC_like': 'OPC/OC-like',
    'AC_MES_like': 'AC/MES-like',
    'Myeloid': 'Myeloid',
    'Neutrophil': 'Neutrophil',
    'T_NK': 'T/NK',
    'B_cell': 'B cell',
    'Endothelial_cell': 'Endothelial cell',
    'Astrocyte': 'Astrocyte',
    'OPC_OC': 'OPC/OC',
}

# Groups: Group1/Group2 internally, PL/TL in every figure
group_labels = {'Group1': 'PL', 'Group2': 'TL'}
group_colors = {'PL': '#990000', 'TL': '#2b8cbe'}
group_order = ['PL', 'TL']

malig_celltype = ['Cycling', 'OPC_OC_like', 'AC_MES_like']
non_malig_celltype = [c for c in celltype_order if c not in malig_celltype]


# ---------------------------------------------------------------------------
# Sample table: the single entry point to the sample manifest
# ---------------------------------------------------------------------------
def load_samples(in_analysis_only=True, group_column='group_as_run'):
    """Read sample_metadata.tsv and return the sample table indexed by CODEX_ID.

    With in_analysis_only=True, only the 31 samples that enter the analysis are
    returned (2 low-quality and 2 recurrent samples are dropped; the reason for each is
    in the exclusion_reason column).

    group_column selects which notion the returned 'group' column carries:
      'group_as_run' (default)  Group1/Group2, matching g1_samples/g2_samples;
                                **required to reproduce the outputs in the repository**
      'group'                   PL/TL from Table S1 of the paper, the authoritative
                                grouping. The two columns currently agree on all 31
                                samples, so this is equivalent today; it is kept
                                separate so that a future change to Table S1 shows up
                                as a difference instead of passing silently.

    In the returned table the group column has already been set from group_column; both
    original columns are still present.
    """
    import pandas as pd
    from pathlib import Path

    path = Path(__file__).resolve().parent / 'sample_metadata.tsv'
    df = pd.read_csv(path, sep='\t', dtype=str)
    df['in_analysis'] = df['in_analysis'] == 'True'
    if in_analysis_only:
        df = df[df['in_analysis']]
    df = df.set_index('CODEX_ID')

    if group_column not in ('group', 'group_as_run'):
        raise ValueError(f"group_column must be 'group' or 'group_as_run', got {group_column!r}")
    df['group'] = df[group_column].map(lambda g: {'PL': 'Group1', 'TL': 'Group2'}.get(g, g))
    return df


def sample_lists(group_column='group_as_run'):
    """Return (group1_ids, group2_ids), ordered as in g1_samples/g2_samples."""
    df = load_samples(group_column=group_column)
    order = {s: i for i, s in enumerate(g1_samples + g2_samples)}
    g1 = sorted(df.index[df['group'] == 'Group1'], key=lambda s: order.get(s, len(order)))
    g2 = sorted(df.index[df['group'] == 'Group2'], key=lambda s: order.get(s, len(order)))
    return g1, g2
# ---------------------------------------------------------------------------
# cell_type_2: the 17-class fine-grained annotation used for the nearest-neighbour
# distance analysis (pipeline/08_min_distance).
#
# It nests inside celltype_order above (the 10 ct_agg classes); verified to have no
# crossings:
#   Myeloid -> Dendritic_cell + Myeloid_CD68/CD163/M2_like/other
#   T_NK    -> Tcell_CD4 + Tcell_CD8 + Treg + NK
#   the other 7 classes map one to one
#
# old2new_celltype has already been applied here too, so OPC_OC below means normal
# oligodendrocyte precursors and the malignant states are OPC_OC_like / AC_MES_like.
# ---------------------------------------------------------------------------
celltype2_order = [
    'Cycling', 'OPC_OC_like', 'AC_MES_like',
    'Dendritic_cell', 'Myeloid_CD68', 'Myeloid_CD163', 'Myeloid_M2_like', 'Myeloid_other',
    'Neutrophil',
    'Tcell_CD4', 'Tcell_CD8', 'Treg', 'NK', 'B_cell',
    'Endothelial_cell', 'Astrocyte', 'OPC_OC',
]

celltype2_display = {
    'Cycling': 'Cycling',
    'OPC_OC_like': 'OPC/OC-like',
    'AC_MES_like': 'AC/MES-like',
    'Tumor': 'Tumor',                       # the merged class under map_tum
    'Dendritic_cell': 'DC',
    'Myeloid': 'Myeloid',                   # the merged class under map_tum_map_mye
    'Myeloid_CD68': 'CD68+ Myeloid',
    'Myeloid_CD163': 'CD163+ Myeloid',
    'Myeloid_M2_like': 'CD68+CD163+ Myeloid',
    'Myeloid_other': 'CD68-CD163- Myeloid',
    'Neutrophil': 'Neutrophil',
    'Tcell_CD4': 'CD4 T cell',
    'Tcell_CD8': 'CD8 T cell',
    'Treg': 'Treg',
    'NK': 'NK',
    'B_cell': 'B cell',
    'Endothelial_cell': 'Vascular',
    'Astrocyte': 'Astrocyte',
    'OPC_OC': 'OPC/OC',
}

# Class order after merging, one entry per variant produced by pipeline/08_min_distance
celltype2_order_map_tum = (
    ['Tumor'] + [c for c in celltype2_order if c not in malig_celltype]
)
celltype2_order_map_tum_map_mye = [
    c for c in celltype2_order_map_tum
    if c not in ('Myeloid_CD68', 'Myeloid_CD163', 'Myeloid_M2_like', 'Myeloid_other')
]
celltype2_order_map_tum_map_mye.insert(
    celltype2_order_map_tum_map_mye.index('Dendritic_cell') + 1, 'Myeloid')
