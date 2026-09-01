#!/usr/bin/env python3
"""Gated prospective shard driver. Compile-only never draws a permutation."""
import argparse,gzip,hashlib,json,os,platform,re,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PINNED=('CONTRACT.json','producer.cpp','run_production.py','score_production.py','archive_channel_split.py','vendor/primitive.cpp','vendor/integer_period.cpp')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def compile_engine(n):
    kind='integer_period' if n==340 else 'primitive';(ROOT/'build').mkdir(exist_ok=True)
    cmd=['g++','-O3','-std=c++17','-fopenmp',f'-DMATCHING_NORM4_BACKEND="{ROOT}/vendor/{kind}.cpp"']
    if n==340:cmd+=['-DMATCHING_NORM4_INTEGER=1']
    dest=ROOT/'build'/kind;cmd += [str(ROOT/'producer.cpp'),'-o',str(dest)]
    p=subprocess.run(cmd,text=True,capture_output=True,check=True)
    return dest,{'command':cmd,'stdout':p.stdout,'stderr':p.stderr,'binary_sha256':sha(dest)}
def authorize(path):
    if not path or not path.is_file():raise RuntimeError('NO FREEZE AUTHORIZATION: no prospective sample may be generated')
    a=json.loads(path.read_text())
    if a.get('authorization')!='root-explicitly-authorized-after-freeze' or not re.fullmatch('[0-9a-f]{40}',a.get('freeze_commit','')):
        raise RuntimeError('Actual root freeze commit and explicit production authorization required')
    for name in PINNED:
        if a.get('sha256',{}).get(name)!=sha(ROOT/name):raise RuntimeError('Frozen input changed: '+name)
    return a
def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,choices=[85,340],required=True)
    p.add_argument('--batch-begin',type=int,default=0);p.add_argument('--batch-end',type=int,default=200)
    p.add_argument('--workers',type=int,default=14);p.add_argument('--authorization',type=Path)
    p.add_argument('--compile-only',action='store_true');args=p.parse_args()
    if args.compile_only:
        binary,record=compile_engine(args.n);print(json.dumps({'status':'compiled_only_no_samples','binary':str(binary),**record}));return
    a=authorize(args.authorization)
    if not 0<=args.batch_begin<args.batch_end<=200 or not 1<=args.workers<=14:raise ValueError('fixed batch/worker limits')
    if args.n==85 and (args.batch_begin,args.batch_end)!=(0,200):raise ValueError('N85 fixed one5M shard')
    if args.n==340 and (args.batch_begin%25 or args.batch_end-args.batch_begin!=25):raise ValueError('N340 fixed25-batch20M shards')
    out=ROOT/'production';out.mkdir(exist_ok=True);tag=f'n{args.n}-b{args.batch_begin:03d}-{args.batch_end:03d}'
    receipt=out/(tag+'.run.json');raw=out/(tag+'.csv')
    if receipt.exists() or raw.exists() or raw.with_suffix('.csv.gz').exists():raise RuntimeError('Existing shard: no duplicate sampling/overwrite')
    record={'status':'starting','freeze_commit':a['freeze_commit'],'authorization_sha256':sha(args.authorization),
            'hostname':platform.node(),'python':platform.python_version(),'started_unix':time.time(),
            'N':args.n,'batch_begin':args.batch_begin,'batch_end':args.batch_end,'workers':args.workers,
            'frozen_sha256':a['sha256']};save(receipt,record);t=time.perf_counter()
    try:
        binary,build=compile_engine(args.n);record['compile']=build;record['status']='running';save(receipt,record)
        cmd=[str(binary),str(args.n),str(args.batch_begin),str(args.batch_end),str(raw),str(args.workers),a['freeze_commit']]
        process=subprocess.run(cmd,text=True,capture_output=True,check=True);record.update(command=cmd,stdout=process.stdout,stderr=process.stderr)
        record['uncompressed_sha256']=sha(raw);compressed=raw.with_suffix('.csv.gz')
        with raw.open('rb') as source,gzip.open(compressed,'wb',compresslevel=6) as target:
            while True:
                block=source.read(1024*1024)
                if not block:break
                target.write(block)
        record['gzip_sha256']=sha(compressed);record['samples']=(25000 if args.n==85 else 800000)*(args.batch_end-args.batch_begin)
        record['status']='completed';raw.unlink();print(process.stdout,flush=True)
    except Exception as e:
        record['status']='failed';record['error']=repr(e)
        if isinstance(e,subprocess.CalledProcessError):record['failure_stdout']=e.stdout;record['failure_stderr']=e.stderr
        raise
    finally:record['elapsed_seconds']=time.perf_counter()-t;record['finished_unix']=time.time();save(receipt,record)
if __name__=='__main__':main()
