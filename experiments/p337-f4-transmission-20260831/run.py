#!/usr/bin/env python3
"""Compile and run exactly one frozen fresh F4 size; standard library only."""
import argparse
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--freeze-commit', required=True)
    ap.add_argument('--threads', type=int, default=14)
    ap.add_argument('--background', action='store_true')
    args = ap.parse_args()
    package = Path(__file__).resolve().parent
    output = package / 'raw'
    output.mkdir(exist_ok=True)
    prefix = output / ('n' + str(args.n))
    if args.background:
        with prefix.with_suffix('.log').open('x') as log:
            command = [sys.executable, str(Path(__file__).resolve()), '--n', str(args.n),
                       '--freeze-commit', args.freeze_commit, '--threads', str(args.threads)]
            process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
        prefix.with_suffix('.pid').write_text(str(process.pid) + '\n')
        print(json.dumps({'pid': process.pid, 'N': args.n, 'log': str(prefix.with_suffix('.log'))}))
        return
    contract = json.loads((package / 'CONTRACT.json').read_text())
    if args.n not in contract['Ns'] or not 1 <= args.threads <= 14:
        raise ValueError('unplanned N or thread count')
    receipt_path = prefix.with_suffix('.run.json')
    if receipt_path.exists() or prefix.with_suffix('.hist.csv').exists():
        raise FileExistsError('refusing to overwrite an existing production domain')
    backend = package.parents[1] / 'src/threshold_rank_orientation_mc.cpp'
    sources = [package / name for name in ('CONTRACT.json', 'producer.cpp', 'score.py', 'run.py')]
    sources.append(backend)
    receipt = {'status': 'running', 'N': args.n, 'started_utc': utc(),
               'freeze_commit': args.freeze_commit, 'hostname': platform.node(),
               'architecture': platform.machine(), 'pid': os.getpid(), 'threads': args.threads,
               'source_sha256': {str(p.relative_to(package.parents[1])): sha(p) for p in sources},
               'contract': contract, 'old_data_pooled': False}
    receipt['cgroup_v1'] = {name: Path(path).read_text().strip() for name, path in {
        'cpu_quota_us': '/sys/fs/cgroup/cpu/cpu.cfs_quota_us',
        'cpu_period_us': '/sys/fs/cgroup/cpu/cpu.cfs_period_us',
        'memory_limit_bytes': '/sys/fs/cgroup/memory/memory.limit_in_bytes'}.items()}
    receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
    begin = time.monotonic()
    try:
        binary = package / 'producer'
        compile_cmd = ['g++', '-O3', '-DNDEBUG', '-std=c++17', '-fopenmp',
                       str(package / 'producer.cpp'), '-o', str(binary)]
        receipt['compiler'] = subprocess.check_output(['g++', '--version'], text=True).splitlines()[0]
        subprocess.run(compile_cmd, check=True)
        receipt['compile_command'] = compile_cmd
        receipt['binary_sha256'] = sha(binary)
        cmd = [str(binary), '--n', str(args.n), '--samples-per-batch', str(contract['samples_per_batch']),
               '--batch-begin', '0', '--batch-end', str(contract['batches']),
               '--seed', str(contract['seeds'][str(args.n)]), '--threads', str(args.threads),
               '--output-prefix', str(prefix), '--freeze-commit', args.freeze_commit]
        receipt['command'] = cmd
        receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
        subprocess.run(['/usr/bin/time', '-v', '-o', str(prefix.with_suffix('.resource.txt'))] + cmd, check=True)
        csv_path = prefix.with_suffix('.hist.csv')
        gzip_path = csv_path.with_suffix(csv_path.suffix + '.gz')
        with csv_path.open('rb') as source, gzip_path.open('xb') as destination:
            with gzip.GzipFile(filename='', mode='wb', fileobj=destination, mtime=0) as compressed:
                shutil.copyfileobj(source, compressed)
        receipt.update(status='completed', exit_code=0,
                       samples=contract['samples_per_N'][str(args.n)],
                       batch_begin=0, batch_end=contract['batches'],
                       csv_sha256=sha(csv_path), gzip_sha256=sha(gzip_path),
                       producer_metadata_sha256=sha(prefix.with_suffix('.metadata.json')))
    except Exception as error:
        receipt.update(status='failed', error=repr(error))
        raise
    finally:
        receipt.update(finished_utc=utc(), elapsed_seconds=time.monotonic()-begin)
        receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
        print(json.dumps({k: receipt.get(k) for k in ('status', 'N', 'samples', 'elapsed_seconds', 'error')}), flush=True)


if __name__ == '__main__':
    main()
