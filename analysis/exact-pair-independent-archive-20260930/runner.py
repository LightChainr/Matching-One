#!/usr/bin/env python3
"""One bounded replay of existing permutations; standard library only."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
FIXED = dict(L=512, N=262144, a=154646, b=155385, batches=14, samples_per_batch=10000)
OLD_HEADER = ['J1', 'J2', 'dx_b', 'dy_b', 'nu_b', 'tau', 'nu_tau']
NEW_HEADER = ['J1', 'rank_b', 'dx_b', 'dy_b', 'nu_b', 'synergy_edges', 'nonlocal_candidates']


def require(ok, message):
    if not ok:
        raise ValueError(message)


def utc():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')


def read(path):
    return json.loads(path.read_text())


def mapped(row):
    j1, j2, dx, dy, c, _, _ = map(int, row)
    rank = 0 if j1 > FIXED['b'] else (2 if j2 <= FIXED['b'] else 1)
    return [0 if rank == 0 else j1, rank, dx, dy, c]


def prepare():
    repo = ROOT.parent.parent
    old = ROOT.parent / 'completion-hazard-production-20260929'
    reference = ROOT.parent / 'exact-completion-pairs-20260929'
    manifest = read(old / 'data/run.json')
    require(manifest['status'] == 'completed', 'source incomplete')
    require(all(manifest['contract'][k] == v for k, v in FIXED.items()), 'fixed constants differ')
    batches = sorted(manifest['batches'], key=lambda b: b['batch'])
    require([b['batch'] for b in batches] == list(range(14)), 'source batch IDs')
    other = read(reference / 'data/run.json')
    require(not ({b['seed'] for b in batches} & {b['seed'] for b in other['batches']}), 'seed overlap')
    paths = [old / 'engine.cpp', reference / 'engine.cpp', old / 'data/run.json',
             reference / 'data/run.json', reference / 'results/result.json']
    paths += [old / 'data' / b['file'] for b in batches]
    records = []
    for src in paths:
        if src.suffix == '.cpp':
            dst = ROOT / 'source' / src.parent.name / src.name
        elif src.parent.name == 'data' and src.parent.parent == old:
            dst = ROOT / 'inputs/archive' / src.name
        elif src.parent.name == 'data':
            dst = ROOT / 'inputs/reference-run.json'
        else:
            dst = ROOT / 'inputs/reference-result.json'
        require(not dst.exists(), 'refuse overwrite: ' + str(dst))
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        relative = str(src.relative_to(repo))
        records.append(dict(original=relative, snapshot=str(dst.relative_to(ROOT)), sha256=digest(dst),
                            last_file_commit=subprocess.check_output(['git', '-C', str(repo), 'log', '-1', '--format=%H', '--', relative], text=True).strip()))
    for batch in batches:
        require(digest(ROOT / 'inputs/archive' / batch['file']) == batch['sha256'], 'source hash mismatch')
    save(ROOT / 'source-record.json', dict(prepared_utc=utc(), fixed=FIXED,
         head_commit=subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip(),
         branch=subprocess.check_output(['git', '-C', str(repo), 'branch', '--show-current'], text=True).strip(),
         old_acquisition_commit=manifest['source_commit'], files=records,
         reference_seeds_disjoint=True, prediction='deltaY < 0', new_random_prefixes=0,
         archive_status='Other outcomes in this existing independent block were previously read; e had not been measured. Not prospective sampling or untouched holdout.'))
    print('Prepared unchanged sources and all 14 archived input batches.', flush=True)


def environment():
    cg = Path('/sys/fs/cgroup')
    def value(*names):
        for name in names:
            path = cg / name
            if path.exists():
                return path.read_text().strip()
        return None
    quota = value('cpu/cpu.cfs_quota_us')
    period = value('cpu/cpu.cfs_period_us')
    maximum = value('cpu.max')
    cpu = int(quota) / int(period) if quota and int(quota) > 0 else None
    if cpu is None and maximum and not maximum.startswith('max'):
        q, p = maximum.split()
        cpu = int(q) / int(p)
    mem = value('memory/memory.limit_in_bytes', 'memory.max')
    require(cpu is not None and mem not in (None, 'max'), 'unmeasured resource limits')
    return dict(utc=utc(), platform=platform.platform(), python=sys.version,
                machine=platform.machine(), hostname=platform.node(), cpu_quota_cores=cpu,
                memory_limit_bytes=int(mem), affinity_cpus=len(os.sched_getaffinity(0)),
                compiler=subprocess.check_output(['chroot', str(ROOT / 'toolchain'), '/usr/bin/g++', '--version'], text=True),
                ps=subprocess.check_output(['ps', '-eo', 'pid,ppid,comm,etime,pcpu,pmem'], text=True),
                thread_environment={k: os.environ.get(k) for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')})


def benchmark():
    env = environment()
    build = ROOT / 'build'
    build.mkdir(exist_ok=False)
    # The container root loses packages at restart. Compile within this task's
    # persistent installroot, leaving every other remote workspace untouched.
    project = ROOT / 'toolchain/project'
    project.mkdir(exist_ok=False)
    shutil.copytree(ROOT / 'source', project / 'source')
    cmd = ['chroot', str(ROOT / 'toolchain'), '/usr/bin/g++', '-O3', '-std=c++17',
           '/project/source/exact-completion-pairs-20260929/engine.cpp', '-o', '/project/engine']
    subprocess.run(cmd, check=True)
    binary = build / 'engine'
    shutil.copyfile(project / 'engine', binary)
    binary.chmod(0o755)
    b = read(ROOT / 'inputs/archive/run.json')['batches'][0]
    out = build / 'same-seed-32.csv'
    run = [str(ROOT / 'toolchain/usr/bin/time'), '-v', str(binary), '512', '32', b['seed'], '155385', str(out)]
    start = time.monotonic()
    proc = subprocess.run(run, text=True, capture_output=True, check=True)
    seconds = time.monotonic() - start
    with gzip.open(ROOT / 'inputs/archive' / b['file'], 'rt') as f, out.open() as g:
        a, z = csv.reader(f), csv.reader(g)
        require(next(a) == OLD_HEADER and next(z) == NEW_HEADER, 'header mismatch')
        rows = list(z)
        require(len(rows) == 32, 'benchmark rows')
        for i, row in enumerate(rows):
            require(mapped(next(a)) == list(map(int, row[:5])), 'benchmark alignment row ' + str(i))
    rss = int(next(line.split(':', 1)[1] for line in proc.stderr.splitlines() if 'Maximum resident set size (kbytes)' in line))
    workers = min(14, math.floor(env['cpu_quota_cores']), env['affinity_cpus'], int(env['memory_limit_bytes'] * .8 // (rss * 1024)))
    require(workers >= 1, 'no feasible workers')
    save(ROOT / 'benchmark.json', dict(environment=env, build_command=cmd, binary_sha256=digest(binary),
         command=run, stdout=proc.stdout, time_log=proc.stderr, wall_seconds=seconds, max_rss_kib=rss,
         same_seed=b['seed'], batch=b['batch'], rows=32, all_mapped_first5_match=True,
         scientific_weight=0, recommended_workers=workers,
         approximate_replay_seconds=seconds / 32 * 10000 * math.ceil(14 / workers)))
    print(json.dumps(read(ROOT / 'benchmark.json')), flush=True)


def replay():
    benchmark_record = read(ROOT / 'benchmark.json')
    env = environment()
    workers = min(benchmark_record['recommended_workers'], math.floor(env['cpu_quota_cores']), env['affinity_cpus'], 14)
    source = read(ROOT / 'inputs/archive/run.json')
    data = ROOT / 'data'
    data.mkdir(exist_ok=False)
    record = dict(status='running', started_utc=utc(), pid=os.getpid(), workers=workers, environment=env,
                  fixed=FIXED, new_random_prefixes=0, source_manifest_sha256=digest(ROOT / 'inputs/archive/run.json'),
                  binary_sha256=digest(ROOT / 'build/engine'), runner_sha256=digest(Path(__file__)), batches=[])
    save(data / 'started.json', record)
    start = time.monotonic()
    def one(batch):
        file = 'exact-square-L512-b%02d.csv' % batch['batch']
        raw = data / file
        command = [str(ROOT / 'build/engine'), '512', '10000', batch['seed'], '155385', str(raw)]
        t0 = time.monotonic()
        proc = subprocess.run(command, capture_output=True, text=True)
        row = dict(batch=batch['batch'], seed=batch['seed'], samples=10000, original_file=batch['file'],
                   dependency_group=batch['dependency_group'], command=command,
                   returncode=proc.returncode, stdout=proc.stdout, stderr=proc.stderr, wall_seconds=time.monotonic()-t0)
        save(data / ('batch-%02d.json' % batch['batch']), row)
        require(proc.returncode == 0, 'engine failed: ' + file)
        compressed = data / (file + '.gz')
        with raw.open('rb') as f, compressed.open('xb') as output:
            with gzip.GzipFile(filename='', mode='wb', fileobj=output, mtime=0) as g:
                shutil.copyfileobj(f, g)
        row.update(file=compressed.name, sha256=digest(compressed))
        print(json.dumps(dict(batch=row['batch'], seconds=row['wall_seconds'], status='completed')), flush=True)
        return row
    try:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(one, b) for b in source['batches']]
            for future in as_completed(futures):
                record['batches'].append(future.result())
        record['batches'].sort(key=lambda b: b['batch'])
        record.update(status='completed', wall_seconds=time.monotonic()-start, finished_utc=utc())
        save(data / 'run.json', record)
        print(json.dumps(dict(status='completed', wall_seconds=record['wall_seconds'])), flush=True)
    except BaseException as exc:
        record.update(status='failed', error=repr(exc), finished_utc=utc())
        save(data / 'failed.json', record)
        raise


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode', choices=['prepare', 'benchmark', 'replay'])
    args = ap.parse_args()
    globals()[args.mode]()
