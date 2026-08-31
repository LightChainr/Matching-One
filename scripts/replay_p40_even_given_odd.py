#!/usr/bin/env python3
"""Compile once and reobserve only P40's frozen configurations, never new ones."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/'analysis/p40_even_given_odd_replay.json'

def sha(content):
    return hashlib.sha256(content).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=ROOT/'results/p40-even-given-odd')
    args=parser.parse_args()
    contract=json.loads(MANIFEST.read_text())
    outputs=[args.output_dir/'raw'/f'n{run["N"]}.csv' for run in contract['runs']]
    receipt=args.output_dir/'run.json'
    if receipt.exists() or any(path.exists() for path in outputs):
        raise ValueError('existing output: inspect it instead of replaying or overwriting')
    source=f'{contract["source_commit"]}:{contract["source_backend_path"]}'
    blob=subprocess.check_output(['git','rev-parse',source],cwd=ROOT,text=True).strip()
    if blob!=contract['source_backend_blob']:
        raise ValueError('archived backend differs from the specified immutable blob')
    archived=subprocess.check_output(['git','show',source],cwd=ROOT)
    compiler=subprocess.check_output(['c++','--version'],text=True).splitlines()[0]
    args.output_dir.joinpath('raw').mkdir(parents=True,exist_ok=True)
    record={'schema':'matching-one.p40-even-replay-run.v1',
            'started_utc':datetime.now(timezone.utc).isoformat(),
            'execution_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'source_backend':{'commit':contract['source_commit'],'path':contract['source_backend_path'],
                              'git_blob':blob,'sha256':sha(archived)},
            'code':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p.read_bytes())} for p in
                    (Path(__file__),MANIFEST,ROOT/'src/p40_even_given_odd_replay.cpp')],
            'environment':{'python':platform.python_version(),'machine':platform.machine(),
                           'platform':platform.platform(),'compiler':compiler},
            'new_samples':0,'server_actions':0,'test_suites':[], 'runs':[]}
    with tempfile.TemporaryDirectory(prefix='matching-p40-even-build-') as temporary:
        build=Path(temporary)
        backend=build/'archived_gaussian.cpp'
        backend.write_bytes(archived)
        binary=build/'p40-even-replay'
        command=['c++','-O3','-std=c++17',f'-DMATCHING_P40_BACKEND="{backend}"',
                 str(ROOT/'src/p40_even_given_odd_replay.cpp'),'-o',str(binary)]
        start=time.perf_counter()
        subprocess.run(command,check=True,cwd=ROOT)
        record['compile']={'command':command,'elapsed_seconds':time.perf_counter()-start,
                           'binary_sha256':sha(binary.read_bytes())}
        for run,output in zip(contract['runs'],outputs):
            command=[str(binary),str(run['N']),str(output)]
            start=time.perf_counter()
            completed=subprocess.run(command,check=True,cwd=ROOT,text=True,capture_output=True)
            record['runs'].append({'N':run['N'],'command':command,
                 'elapsed_seconds':time.perf_counter()-start,'stdout':completed.stdout,
                 'output':str(output.relative_to(ROOT)) if output.is_relative_to(ROOT) else str(output),
                 'sha256':sha(output.read_bytes()),'old_counters':contract['counter_interval']})
            receipt.write_text(json.dumps(record,indent=2)+'\n')
            print(completed.stdout.strip(),flush=True)
    record['completed_utc']=datetime.now(timezone.utc).isoformat()
    receipt.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'runs':record['runs'],'new_samples':0},indent=2))

if __name__=='__main__':
    main()
