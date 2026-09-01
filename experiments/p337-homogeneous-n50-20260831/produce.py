#!/usr/bin/env python3
"""Run one frozen finite population with explicit CPU/RSS/state limits."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

PACKAGE=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--freeze-commit',required=True)
    ap.add_argument('--geometry',choices=('first','second'),required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args()
    repo=Path(subprocess.check_output(['git','-C',str(PACKAGE),'rev-parse','--show-toplevel'],text=True).strip())
    freeze=subprocess.check_output(['git','-C',str(repo),'rev-parse',args.freeze_commit+'^{commit}'],text=True).strip()
    if freeze!=args.freeze_commit:raise ValueError('full freeze SHA required')
    geometry=[5,5] if args.geometry=='first' else [1,7]
    name='producer/geometry_%s_%s.txt'%tuple(geometry)
    files=('CONTRACT.md','score.py','produce.py','producer/frontier.cpp',name)
    hashes={}
    for file in files:
        path=PACKAGE/file
        blob=subprocess.check_output(['git','-C',str(repo),'show',freeze+':'+path.relative_to(repo).as_posix()])
        if blob!=path.read_bytes():raise ValueError('producer differs from freeze: '+file)
        hashes[file]=hashlib.sha256(blob).hexdigest()
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    with tempfile.TemporaryDirectory(prefix='p337-compiled-') as tmp:
        executable=Path(tmp)/'frontier'
        compile_command=['clang++','-O3','-std=c++17',str(PACKAGE/'producer/frontier.cpp'),'-o',str(executable)]
        subprocess.run(compile_command,check=True)
        binary_sha=hashlib.sha256(executable.read_bytes()).hexdigest()
        command=[str(executable),'--graph',str(PACKAGE/name),'--output-prefix',str(out/'population'),
                 '--authorization-commit',freeze,'--cpu-seconds','55','--rss-mib','2048','--state-cap','5000000']
        subprocess.run(command,check=True)
    meta=json.loads((out/'population.json').read_text())
    meta.update(producer_freeze=freeze,producer_hashes=hashes,binary_sha256=binary_sha,
                compile_command=compile_command,command=command)
    if meta['complete']:
        rows=list(csv.DictReader((out/'population.csv').open()))
        meta['histogram']=[{key:int(value) for key,value in row.items()} for row in rows]
    (out/'table.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps({key:meta[key] for key in ('geometry','complete','stop','completed_layers','peak_rss_mib')},indent=2))


if __name__=='__main__':main()
