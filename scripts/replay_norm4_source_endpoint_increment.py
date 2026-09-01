#!/usr/bin/env python3
"""Add cluster marks to 900k unused OLD permutations at each norm-4 endpoint."""
from __future__ import annotations

import json
import platform
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import replay_norm4_source_thermal as old

ROOT = old.ROOT
CONTRACT = ROOT / 'analysis/norm4_source_endpoint_1m_contract.json'
OUTPUT = ROOT / 'results/norm4-source-endpoint-1m'


def main():
    contract = json.loads(CONTRACT.read_text())
    receipt = OUTPUT / 'run.json'
    paths = {n: OUTPUT / 'increment/raw' / f'n{n}.csv' for n in contract['Ns']}
    if receipt.exists() or any(path.exists() for path in paths.values()):
        raise ValueError('endpoint increment already exists; do not repeat or overwrite')
    original_runs = {n: old.source_run(n) for n in contract['Ns']}
    for n, run in original_runs.items():
        if run['old_counter_interval'] != contract['original_counter_interval'] or run['seed'] != contract['seeds'][str(n)]:
            raise ValueError('increment differs from immutable original production')
    backend = old.git_bytes(old.SOURCE_COMMIT, old.BACKENDS['integer_period'])
    record = {
        'schema': 'matching-one.norm4-source-endpoint-increment-run.v1',
        'status': 'running', 'started_utc': old.utc_now(), 'contract': contract,
        'execution_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'original_runs': original_runs,
        'source_backend': {'commit': old.SOURCE_COMMIT, 'path': old.BACKENDS['integer_period'], 'sha256': old.sha(backend)},
        'code': [{'path': old.display_path(path), 'sha256': old.sha(path.read_bytes())}
                 for path in (Path(__file__), old.CPP, CONTRACT)],
        'environment': {'python': platform.python_version(), 'machine': platform.machine(),
                        'compiler': subprocess.check_output(['c++', '--version'], text=True).splitlines()[0]},
        'new_samples': 0, 'old_permutation_reobservations': 1800000,
        'repeated_already_marked_permutations': 0, 'workers': 2,
        'server_actions': 0, 'test_suites': [], 'runs': []}
    OUTPUT.joinpath('increment/raw').mkdir(parents=True, exist_ok=True)
    with receipt.open('x') as handle:
        handle.write(json.dumps(record, indent=2) + '\n')

    def save():
        receipt.write_text(json.dumps(record, indent=2) + '\n')

    started = time.perf_counter()
    try:
        with tempfile.TemporaryDirectory(prefix='matching-norm4-source-endpoint-') as directory:
            build = Path(directory)
            archive = build / 'archived_integer.cpp'
            archive.write_bytes(backend)
            binary = build / 'endpoint-replay'
            command = ['c++', '-O3', '-std=c++17', f'-DMATCHING_NORM4_BACKEND="{archive}"',
                       '-DMATCHING_NORM4_INTEGER=1', str(old.CPP), '-o', str(binary)]
            begin = time.perf_counter()
            process = subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=True)
            record['compile'] = {'command': command, 'elapsed_seconds': time.perf_counter() - begin,
                                 'stdout': process.stdout, 'stderr': process.stderr,
                                 'binary_sha256': old.sha(binary.read_bytes())}
            save()

            def replay(n):
                command = [str(binary), str(n), str(paths[n]), '100000', '900000']
                begin = time.perf_counter()
                process = subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=True)
                return {'N': n, 'command': command, 'elapsed_seconds': time.perf_counter() - begin,
                        'stdout': process.stdout, 'stderr': process.stderr,
                        'output': old.display_path(paths[n]), 'sha256': old.sha(paths[n].read_bytes()),
                        'reobserved_counter_interval': contract['increment_counter_interval'],
                        'seed': contract['seeds'][str(n)], 'old_permutation_reobservations': 900000,
                        'new_samples': 0}

            with ThreadPoolExecutor(max_workers=2) as pool:
                for future in as_completed([pool.submit(replay, n) for n in contract['Ns']]):
                    row = future.result()
                    record['runs'].append(row)
                    save()
                    print(f"N{row['N']}: added 900k old-permutation source marks in {row['elapsed_seconds']:.2f}s", flush=True)
        record['status'] = 'completed'
    except Exception as error:
        record['status'] = 'failed'
        record['error'] = repr(error)
        raise
    finally:
        record['elapsed_seconds'] = time.perf_counter() - started
        record['finished_utc'] = old.utc_now()
        save()


if __name__ == '__main__':
    main()
