#!/usr/bin/env python
"""
cNMF consensus at the selected K, and the top 100 genes of every program

Input : ../../results/malignant_states/03_cnmf/<sid>/   cNMF runs of step 03
        ../../results/malignant_states/selected_K_plateau.txt   K per sample (03b_select_k.py)
Output: ../../results/malignant_states/03_cnmf/<sid>/<sid>.usages.k_<K>.dt_0_02.consensus.txt
            usage of every program in every cell (cells x programs)
        ../../results/malignant_states/03_cnmf/<sid>.programs.genes100.k_<K>.csv
            the top 100 genes of every program
        A sample whose products already exist is skipped (--overwrite recomputes).

The consensus is run at two local density thresholds per sample:
  2.00  no filtering, to see how many outlier replicates there are
  0.02  strict filtering of outlier replicates   <-- the only one step 06 uses

Usage:
    python 04_run_consensus.py [sample ...] [--overwrite]

Runtime about 5 min for the 19 samples. cnmf 1.7.1 (base anaconda), pandas 3.
Paths are relative to this directory: run the script from here.
"""
import os
import sys
import time
import traceback
import warnings

import pandas as pd

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    pd.set_option("future.infer_string", False)

from cnmf import cNMF  # noqa: E402

RES_DIR = "../../results/malignant_states"
OUT_DIR = os.path.join(RES_DIR, "03_cnmf")
KFILE = os.path.join(RES_DIR, "selected_K_plateau.txt")
DENSITY_THRESHOLDS = [2.00, 0.02]
DT_FINAL = 0.02


def read_kfile(path):
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                parts = line.split()
                out.append((parts[0], int(parts[1])))
    if not out:
        sys.exit("no sample in %s" % path)
    return out


def consensus_done(obj, k, dt):
    dt_s = str(dt).replace(".", "_")
    keys = ("consensus_usages__txt", "consensus_spectra__txt", "gene_spectra_score__txt",
            "gene_spectra_tpm__txt", "clustering_plot")
    return all(os.path.exists(obj.paths[key] % (k, dt_s)) for key in keys)


def run_one(sid, k, overwrite):
    obj = cNMF(output_dir=OUT_DIR, name=sid)
    for dt in DENSITY_THRESHOLDS:
        if overwrite or not consensus_done(obj, k, dt):
            obj.consensus(k=k, density_threshold=dt, show_clustering=True, close_clustergram_fig=True)
    out_csv = os.path.join(OUT_DIR, "%s.programs.genes100.k_%d.csv" % (sid, k))
    if overwrite or not os.path.exists(out_csv):
        _, _, _, topgenes = obj.load_results(K=k, density_threshold=DT_FINAL)
        topgenes.head(100).to_csv(out_csv)
    return out_csv


def main():
    argv = sys.argv[1:]
    overwrite = "--overwrite" in argv
    only = [a for a in argv if a != "--overwrite"]
    ktab = read_kfile(KFILE)
    if only:
        ktab = [(s, k) for s, k in ktab if s in only]
    print("K from: %s\noutput to: %s\n%d samples: %s" % (KFILE, OUT_DIR, len(ktab),
          " ".join("%s(K=%d)" % x for x in ktab)), flush=True)
    failed = []
    for i, (sid, k) in enumerate(ktab, 1):
        t0 = time.time()
        try:
            csv = run_one(sid, k, overwrite)
            print("[%d/%d] %s K=%d done (%.1f min) -> %s" % (i, len(ktab), sid, k, (time.time() - t0) / 60, csv), flush=True)
        except Exception:
            traceback.print_exc()
            failed.append(sid)
    print("\n04 done.")
    if failed:
        print("!! these samples did not finish: %s" % " ".join(failed))
        sys.exit(1)


if __name__ == "__main__":
    main()
