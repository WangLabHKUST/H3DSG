"""Shared by all steps: read the configuration and the sample table, and agree on the output folder of each step."""
import argparse
import importlib.util
import os
import sys
from pathlib import Path

import pandas as pd

PIPE_DIR = Path(__file__).resolve().parent.parent

# output subfolder of each step (relative to OUT_DIR)
STEP_DIRS = {1: 's1_loom', 2: 's2_inputs', 3: 's3_h5ad', 4: 's4_velocity', 5: 's5_embedding', 6: 's6_figures'}


def load_config(path=None):
    path = Path(path or os.environ.get('VELO_PIPE_CONFIG') or PIPE_DIR / 'config.py').resolve()
    spec = importlib.util.spec_from_file_location('velo_pipe_config', path)
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)
    cfg.CONFIG_PATH = path
    return cfg


def parse_args(description, extra=None):
    p = argparse.ArgumentParser(description=description)
    p.add_argument('--config', help='configuration file; default config.py (or the environment variable VELO_PIPE_CONFIG)')
    p.add_argument('--force', action='store_true', help='recompute even if the output exists')
    p.add_argument('--stems', nargs='*', help='overrides config.STEMS')
    if extra:
        extra(p)
    args = p.parse_args()
    cfg = load_config(args.config)
    if args.stems:
        cfg.STEMS = args.stems
    return cfg, args


def step_dir(cfg, n):
    d = Path(cfg.OUT_DIR) / STEP_DIRS[n]
    d.mkdir(parents=True, exist_ok=True)
    return d


def read_samples(cfg):
    s = pd.read_csv(cfg.SAMPLE_SHEET, sep='\t', dtype={'margin': int})
    need = {'sample', 'seurat_id', 'bam', 'barcodes', 'margin'}
    assert need <= set(s.columns), f'samples.tsv lacks columns: {need - set(s.columns)}'
    assert s['sample'].is_unique and s['seurat_id'].is_unique, 'sample and seurat_id must be unique in samples.tsv'
    # bam and barcodes are relative to the folder of the sample table
    base = Path(cfg.SAMPLE_SHEET).resolve().parent
    for c in ('bam', 'barcodes'):
        s[c] = [str(base / p) for p in s[c]]
    return s


def stem_samples(cfg, stem, samples):
    """Samples of an analysis set, in the order of the sample table."""
    if stem not in cfg.STEM_DEFS:
        raise SystemExit(f'unknown analysis set {stem}; config.STEM_DEFS has {list(cfg.STEM_DEFS)}')
    if cfg.STEM_DEFS[stem]['exclude_margin']:
        samples = samples[samples['margin'] == 0]
    return samples['sample'].tolist()


def skip_if_done(out, force, what):
    """True if the output exists and --force was not given."""
    if Path(out).exists() and not force:
        print(f'[skip] {what}: {out} exists (--force recomputes it)', flush=True)
        return True
    return False


def atomic_write(df, out):
    """Writes a temporary file, then renames it, so that a killed run leaves no half-written CSV."""
    out = Path(out)
    tmp = out.with_name(out.name + '.tmp')
    df.to_csv(tmp, index=False)
    tmp.replace(out)


def pin_threads():
    """UMAP and figures run on one thread, so that the result is reproducible. Call before importing numpy."""
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
        os.environ.setdefault(v, '1')
