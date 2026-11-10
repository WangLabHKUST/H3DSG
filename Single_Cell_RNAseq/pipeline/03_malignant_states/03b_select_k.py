#!/usr/bin/env python
"""
Choose K for every sample by the plateau-max rule

Input : ../../results/malignant_states/03_cnmf/<sid>/<sid>.k_selection_stats.df.npz
        ../../sample_metadata.tsv   (renamed_id, written as a comment next to each sample)
Output: ../../results/malignant_states/selected_K_plateau.txt
            sid<TAB>K; the renamed_id is only a "#" comment, because 06_metaprograms.R reads
            the file with a two-column read.table(comment.char = "#")
        ../../results/malignant_states/03_cnmf/k_selection_summary.csv
            silhouette and prediction error of every K of every sample

Rule: the LARGEST K whose silhouette (cNMF stability) is at least 0.95 x the maximum
silhouette of that sample.
Why: a K that is too large only adds small fragment programs, which step 06 drops
automatically (dominant cells < 10 or < 1% of the sample); a K that is too small merges two
real states into one program and cannot be repaired afterwards. The S-phase program also
needs a larger K to separate. The prediction error falls monotonically with K and is only
listed for reference.

With the 19 samples the rule gives 121 programs in total.
Paths are relative to this directory: run the script from here.
"""
import glob
import os
import re
import warnings

import pandas as pd

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    pd.set_option("future.infer_string", False)

from cnmf.cnmf import load_df_from_npz  # noqa: E402

RES_DIR = "../../results/malignant_states"
CNMF_DIR = os.path.join(RES_DIR, "03_cnmf")
KFILE = os.path.join(RES_DIR, "selected_K_plateau.txt")
SAMPLE_META_FILE = "../../sample_metadata.tsv"
TOL = 0.05


def natkey(sid):
    """PL_1 < PL_2 < PL_10"""
    return [int(p) if p.isdigit() else p for p in re.split(r"(\d+)", sid)]


def main():
    files = sorted(glob.glob(os.path.join(CNMF_DIR, "*", "*.k_selection_stats.df.npz")),
                   key=lambda f: natkey(os.path.basename(os.path.dirname(f))))
    if not files:
        raise SystemExit("No k_selection_stats found, step 03 has not finished")
    m = pd.read_csv(SAMPLE_META_FILE, sep="\t", dtype=str)
    renamed = dict(zip(m["sid"], m["renamed_id"]))

    rows, picks = [], []
    for f in files:
        sid = os.path.basename(os.path.dirname(f))
        d = load_df_from_npz(f).sort_values("k")
        ks = [int(x) for x in d["k"]]
        sil = list(d["silhouette"])
        cutoff = (1 - TOL) * max(sil)
        k = max(ks[i] for i in range(len(ks)) if sil[i] >= cutoff)
        picks.append((sid, k))
        print("%-8s K=%-3d silhouette=%s" % (sid, k, " ".join("%d:%.3f" % (a, b) for a, b in zip(ks, sil))))
        for a, b, e in zip(ks, sil, d["prediction_error"]):
            rows.append(dict(sample=sid, renamed_id=renamed.get(sid, sid), k=a, silhouette=b,
                             prediction_error=e, selected=(a == k)))

    pd.DataFrame(rows).to_csv(os.path.join(CNMF_DIR, "k_selection_summary.csv"), index=False)
    with open(KFILE, "w", encoding="utf-8") as fh:
        fh.write("# written by 03b_select_k.py, rule = plateau-max (largest K with silhouette >= 0.95 x maximum)\n")
        fh.write("# columns: sid  K   # renamed_id (comment only: step 06 reads this file with a two-column read.table)\n")
        for sid, k in picks:
            rid = renamed.get(sid, sid)
            fh.write("%s\t%d\t# %s\n" % (sid, k, rid) if rid != sid else "%s\t%d\n" % (sid, k))
    print("K of %d samples written to %s (%d programs in total)" % (len(picks), KFILE, sum(k for _, k in picks)))


if __name__ == "__main__":
    main()
