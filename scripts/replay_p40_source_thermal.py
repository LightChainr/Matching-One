#!/usr/bin/env python3
"""Collect K-stratified source moments from the immutable old P40 counters."""
from __future__ import annotations
import argparse
import json
import platform
import subprocess
import tempfile
import time
from datetime import datetime,timezone
from pathlib import Path
from replay_p40_even_given_odd import sha

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/'analysis/p40_source_thermal_replay.json'

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=ROOT/'results/p40-source-thermal')
    args=parser.parse_args()
    contract=json.loads(MANIFEST.read_text())
    outputs=[args.output_dir/'raw'/f'n{r["N"]}.csv' for r in contract['runs']]
    receipt=args.output_dir/'run.json'
    if receipt.exists() or any(p.exists() for p in outputs):
        raise ValueError('existing output: inspect rather than repeat or overwrite')
    ref=f'{contract["source_commit"]}:{contract["source_backend_path"]}'
    blob=subprocess.check_output(['git','rev-parse',ref],cwd=ROOT,text=True).strip()
    if blob!=contract['source_backend_blob']:
        raise ValueError('immutable source backend does not match')
    backend=subprocess.check_output(['git','show',ref],cwd=ROOT)
    args.output_dir.joinpath('raw').mkdir(parents=True,exist_ok=True)
    record={'schema':'matching-one.p40-source-thermal-run.v1',
       'started_utc':datetime.now(timezone.utc).isoformat(),
       'execution_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
       'backend':{'commit':contract['source_commit'],'path':contract['source_backend_path'],
                  'git_blob':blob,'sha256':sha(backend)},
       'code':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p.read_bytes())} for p in
               (Path(__file__),MANIFEST,ROOT/'src/p40_source_thermal_replay.cpp',ROOT/'scripts/replay_p40_even_given_odd.py')],
       'environment':{'python':platform.python_version(),'machine':platform.machine(),
                      'platform':platform.platform(),'compiler':subprocess.check_output(['c++','--version'],text=True).splitlines()[0]},
       'new_samples':0,'server_actions':0,'test_suites':[],'runs':[]}
    with tempfile.TemporaryDirectory(prefix='matching-p40-source-thermal-') as directory:
        build=Path(directory)
        archived=build/'archived_gaussian.cpp'
        archived.write_bytes(backend)
        binary=build/'p40-source-thermal-replay'
        command=['c++','-O3','-std=c++17',f'-DMATCHING_P40_BACKEND="{archived}"',
                 str(ROOT/'src/p40_source_thermal_replay.cpp'),'-o',str(binary)]
        begin=time.perf_counter()
        subprocess.run(command,check=True,cwd=ROOT)
        record['compile']={'command':command,'elapsed_seconds':time.perf_counter()-begin,
                           'binary_sha256':sha(binary.read_bytes())}
        for run,output in zip(contract['runs'],outputs):
            command=[str(binary),str(run['N']),str(output)]
            begin=time.perf_counter()
            process=subprocess.run(command,check=True,cwd=ROOT,text=True,capture_output=True)
            record['runs'].append({'N':run['N'],'command':command,
                'elapsed_seconds':time.perf_counter()-begin,'stdout':process.stdout,
                'output':str(output.relative_to(ROOT)) if output.is_relative_to(ROOT) else str(output),
                'sha256':sha(output.read_bytes()),'old_counter_interval':contract['counter_interval']})
            receipt.write_text(json.dumps(record,indent=2)+'\n')
            print(process.stdout.strip(),flush=True)
    record['completed_utc']=datetime.now(timezone.utc).isoformat()
    receipt.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'runs':record['runs'],'new_samples':0},indent=2))

if __name__=='__main__':
    main()
