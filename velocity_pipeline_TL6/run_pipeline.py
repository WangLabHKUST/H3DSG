#!/usr/bin/env python
"""Runs the whole pipeline, or part of it: BAM -> loom -> h5ad -> velocity -> UMAP -> figures.

    python run_pipeline.py                     # steps 1-6; a step whose output exists is skipped
    python run_pipeline.py --from 2            # from step 2 (the looms exist)
    python run_pipeline.py --from 5 --to 6     # only UMAP and figures
    python run_pipeline.py --to 1 --dry-run    # only print the velocyto commands
    python run_pipeline.py --force --from 6    # recompute step 6
    python run_pipeline.py --stems ZDK-inc-edge ZDK-exc-edge

Run it in the velo environment (environment.yml): conda activate velo
The output of each step also goes to <OUT_DIR>/logs/s<N>.log; package versions and parameters to <OUT_DIR>/run_info.json.
"""
import argparse
import json
import platform
import subprocess
import sys
import time
from importlib import metadata
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'steps'))
from _common import load_config  # noqa: E402

STEPS = {1: 's1_velocyto.py', 2: 's2_prepare_inputs.py', 3: 's3_build_h5ad.py',
         4: 's4_velocity.py', 5: 's5_umap.py', 6: 's6_figures.py'}
PKGS = ['velocyto', 'scvelo', 'scanpy', 'anndata', 'umap-learn', 'scikit-learn', 'numpy', 'pandas', 'matplotlib']


def snapshot(cfg, argv):
    def ver(p):
        try:
            return metadata.version(p)
        except metadata.PackageNotFoundError:
            return None
    keep = ('STEMS', 'MOMENTS', 'VELOCITY_MODE', 'VELOCITY_GRAPH_NEIGHBORS', 'UMAP', 'STREAM_KW', 'FIG_COLORS',
            'FIG_PREFIX', 'FIG_DPI', 'FIG_ALPHA', 'NN_K', 'GTF', 'META_CSV', 'PCA_CSV', 'UMAP_CSV')
    out = Path(cfg.OUT_DIR)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'run_info.json').write_text(json.dumps(dict(
        time=time.strftime('%Y-%m-%d %H:%M:%S'), argv=argv, python=platform.python_version(),
        packages={p: ver(p) for p in PKGS}, config={k: getattr(cfg, k) for k in keep if hasattr(cfg, k)}),
        indent=2, default=str))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--from', dest='first', type=int, default=1, choices=STEPS)
    p.add_argument('--to', dest='last', type=int, default=6, choices=STEPS)
    p.add_argument('--config')
    p.add_argument('--force', action='store_true')
    p.add_argument('--dry-run', action='store_true', help='step 1 only: print the velocyto commands without running them')
    p.add_argument('--stems', nargs='*')
    args = p.parse_args()

    cfg = load_config(args.config)
    log_dir = Path(cfg.OUT_DIR) / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    snapshot(cfg, sys.argv)

    for n in range(args.first, args.last + 1):
        cmd = [sys.executable, '-u', str(HERE / 'steps' / STEPS[n])]
        if args.config:
            cmd += ['--config', args.config]
        if args.force:
            cmd.append('--force')
        if args.stems:
            cmd += ['--stems', *args.stems]
        if args.dry_run and n == 1:
            cmd.append('--dry-run')
        print(f'\n===== step {n}: {STEPS[n]} =====', flush=True)
        t0 = time.time()
        with open(log_dir / f's{n}.log', 'w') as log:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in proc.stdout:
                sys.stdout.write(line)
                log.write(line)
            rc = proc.wait()
        if rc != 0:
            sys.exit(f'step {n} failed (exit code {rc}), log: {log_dir / f"s{n}.log"}')
        print(f'===== step {n} done in {time.time() - t0:.0f}s =====', flush=True)


if __name__ == '__main__':
    main()
