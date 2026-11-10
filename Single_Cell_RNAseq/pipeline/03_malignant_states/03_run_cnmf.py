#!/usr/bin/env python
"""
cNMF of the malignant cells, one sample at a time, K = 3..10

Input : ../../results/malignant_states/03_cnmf_input/<sid>_filtered_malignant.tsv
        raw counts of the malignant cells (cells x genes, tab separated), written by
        02_filter_malignant.R
Output: ../../results/malignant_states/03_cnmf/<sid>/<sid>.k_selection.png,
        <sid>.k_selection_stats.df.npz and cnmf_tmp/
        A sample that already has its k_selection.png is skipped.

Every sample is decomposed on its own: never run the samples together (cNMF on a merged
matrix returns programs that track patient and batch, not cell state). Programs are
compared across samples afterwards, in 06_metaprograms.R.

cNMF parameters: K = 3..10, 100 iterations per K, 2000 over-dispersed genes, Frobenius
loss, seed 42. The per-iteration spectra files are deleted after combine().

Usage:
    python 03_run_cnmf.py                  # all samples found in 03_cnmf_input
    python 03_run_cnmf.py PL_1 TL_2        # only the listed samples
    NWORKERS=18 python 03_run_cnmf.py      # number of worker processes (default 16)
The 19 samples were run as two such processes at once (the samples sorted by matrix size and
dealt alternately to the two), NWORKERS=18 each, on a 40-core machine.

Runtime about 25 min for the 19 samples in that setup. cnmf 1.7.1 (base anaconda), pandas 3.
Paths are relative to this directory: run the script from here.
"""
import glob
import os
import sys
import time
import traceback
import warnings

import pandas as pd

# pandas 3 + anndata 0.10.9: the h5ad writer does not accept ArrowStringArray. This must be
# set at module top level (Windows spawns the worker processes, which import this module again).
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    pd.set_option("future.infer_string", False)

from cnmf import cNMF  # noqa: E402

RES_DIR = "../../results/malignant_states"
IN_DIR = os.path.join(RES_DIR, "03_cnmf_input")
OUT_DIR = os.path.join(RES_DIR, "03_cnmf")
SUFFIX = "_filtered_malignant.tsv"

K_RANGE = list(range(3, 11))   # K scanned: 3-10
N_ITER = 100                   # replicates per K
NUMGENES = 2000                # over-dispersed genes
SEED = 42
BETA_LOSS = "frobenius"
NWORKERS = int(os.environ.get("NWORKERS", 16))


def list_samples(argv):
    if argv:
        return argv
    files = sorted(glob.glob(os.path.join(IN_DIR, "*" + SUFFIX)))
    if not files:
        sys.exit("No *%s found in %s, run 02_filter_malignant.R first" % (SUFFIX, IN_DIR))
    return [os.path.basename(f)[: -len(SUFFIX)] for f in files]


def run_one(sid):
    obj = cNMF(output_dir=OUT_DIR, name=sid)
    obj.prepare(counts_fn=os.path.join(IN_DIR, sid + SUFFIX), components=K_RANGE, n_iter=N_ITER,
                seed=SEED, num_highvar_genes=NUMGENES, beta_loss=BETA_LOSS)
    obj.factorize_multi_process(total_workers=NWORKERS)
    obj.combine()
    # one spectra file per K x iteration; useless after combine and large
    for f in glob.glob(os.path.join(OUT_DIR, sid, "cnmf_tmp", sid + ".spectra.k_*.iter_*.df.npz")):
        os.remove(f)
    obj.k_selection_plot(close_fig=True)


def main():
    samples = list_samples(sys.argv[1:])
    print("input: %s\noutput: %s" % (IN_DIR, OUT_DIR), flush=True)
    print("samples (%d): %s" % (len(samples), " ".join(samples)), flush=True)
    print("NWORKERS=%d  K=%s  n_iter=%d\n" % (NWORKERS, K_RANGE, N_ITER), flush=True)
    failed = []
    for i, sid in enumerate(samples, 1):
        if os.path.exists(os.path.join(OUT_DIR, sid, sid + ".k_selection.png")):
            print("[%d/%d] %s already has k_selection.png, skipped" % (i, len(samples), sid), flush=True)
            continue
        print("\n========== [%d/%d] %s ==========" % (i, len(samples), sid), flush=True)
        t0 = time.time()
        try:
            run_one(sid)
            print(">> %s done in %.1f min" % (sid, (time.time() - t0) / 60), flush=True)
        except Exception:
            traceback.print_exc()
            failed.append(sid)
            print("!! %s failed, continuing with the next one" % sid, flush=True)
    print("\n03 done.")
    if failed:
        print("!! these samples did not finish: %s" % " ".join(failed))
        sys.exit(1)


if __name__ == "__main__":
    # Windows spawns the worker processes, which import this module again: keep this guard
    main()
