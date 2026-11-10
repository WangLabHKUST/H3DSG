"""Step 1: BAM -> loom (spliced / unspliced / ambiguous), one `velocyto run` per sample.

Input   bam (possorted_genome_bam.bam of Cell Ranger) and barcodes (barcodes.tsv.gz of the filtered matrix, used only
        to select the cells; the count matrix itself is not used, spliced and unspliced counts all come from the BAM)
        of samples.tsv, and config.GTF
Output  results/s1_loom/<sample>/<sample>.loom  +  velocyto.log

Default logic, no mask. -e <sample> fixes the name of the loom file (without it: possorted_genome_bam_<5 random
characters>.loom; step 3 reads *.loom and accepts both), and -@ / --samtools-memory speed up the sorting of the BAM.

velocyto needs a BAM sorted by cell barcode. If cellsorted_<bam name> is not in the folder of the BAM it writes one
there (about 25 GB per sample here, several hours). A sample that already has a loom is skipped.
"""
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from _common import parse_args, read_samples, skip_if_done, step_dir


def velocyto_bin():
    # prefer the velocyto of the conda environment of this python over another version on PATH
    cand = Path(sys.executable).parent / 'velocyto'
    return str(cand) if cand.exists() else (shutil.which('velocyto') or 'velocyto')


def build_cmd(cfg, row, out):
    return [velocyto_bin(), 'run',
            '-b', str(row.barcodes),
            '-o', str(out),
            '-e', row.sample,
            '-@', str(cfg.VELOCYTO_SAMTOOLS_THREADS),
            '--samtools-memory', str(cfg.VELOCYTO_SAMTOOLS_MEM_MB),
            str(row.bam), str(cfg.GTF)]


def run_one(cfg, row, root, force, dry):
    out = root / row.sample
    existing = sorted(out.glob('*.loom'))
    if existing and not force:
        print(f'[skip] {row.sample}: {existing[0]} exists', flush=True)
        return row.sample, 0
    for label, p in (('BAM', row.bam), ('barcodes', row.barcodes), ('GTF', cfg.GTF)):
        if not Path(p).exists():
            raise FileNotFoundError(f'[{row.sample}] {label} not found: {p}')
    out.mkdir(parents=True, exist_ok=True)
    cmd = build_cmd(cfg, row, out)
    print(f"[run] {' '.join(cmd)}\n", end='', flush=True)
    if dry:
        return row.sample, 0
    with open(out / 'velocyto.log', 'w') as log:
        rc = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT)
    if rc != 0:
        raise RuntimeError(f'[{row.sample}] velocyto exit code {rc}, see {out / "velocyto.log"}')
    return row.sample, rc


def main():
    cfg, args = parse_args(__doc__.split('\n')[0],
                           extra=lambda p: p.add_argument('--dry-run', action='store_true', help='only print the commands'))
    samples = read_samples(cfg)
    root = step_dir(cfg, 1)
    with ThreadPoolExecutor(max_workers=cfg.VELOCYTO_PARALLEL) as ex:
        futs = [ex.submit(run_one, cfg, r, root, args.force, args.dry_run) for r in samples.itertuples()]
        for f in futs:
            name, _ = f.result()
            print(f'[done] {name}', flush=True)


if __name__ == '__main__':
    main()
