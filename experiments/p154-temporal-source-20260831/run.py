#!/usr/bin/env python3
"""Build and acquire only the frozen old-counter joint marks, then score them."""
import csv, gzip, hashlib, json, os, platform, subprocess, time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def main():
    out=ROOT/'results';out.mkdir(exist_ok=True)
    if (out/'run.json').exists():raise RuntimeError('Acquisition already exists; refusing replay')
    contract=json.loads((ROOT/'CONTRACT.json').read_text())
    record={'status':'running','started_unix':time.time(),'contract':contract,
            'hostname':platform.node(),'platform':platform.platform(),'python':platform.python_version(),
            'cpu_count':os.cpu_count(),'threads':16,'new_counters':0,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and 'results' not in p.parts},
            'runs':[],'compile':[]}
    save(out/'run.json',record);started=time.perf_counter()
    (out/'raw').mkdir();(ROOT/'build').mkdir(exist_ok=True)
    try:
        for kind in ('primitive','integer_period'):
            cmd=['g++','-O3','-std=c++17','-fopenmp',f'-DMATCHING_NORM4_BACKEND="{ROOT}/vendor/{kind}.cpp"']
            if kind=='integer_period':cmd+=['-DMATCHING_NORM4_INTEGER=1']
            cmd += [str(ROOT/'replay.cpp'),'-o',str(ROOT/'build'/kind)]
            t=time.perf_counter();p=subprocess.run(cmd,text=True,capture_output=True,check=True)
            record['compile'].append({'command':cmd,'seconds':time.perf_counter()-t,'stdout':p.stdout,'stderr':p.stderr,'sha256':sha(ROOT/'build'/kind)})
        old=np.load(ROOT/'inputs/old_profiles.npz')
        for n in contract['Ns']:
            kind='integer_period' if n in (260,340) else 'primitive'
            path=out/'raw'/f'n{n}.csv';cmd=[str(ROOT/'build'/kind),str(n),str(path),'16']
            t=time.perf_counter();p=subprocess.run(cmd,text=True,capture_output=True,check=True)
            values=np.zeros((100,2,n+1,5),dtype=np.int64);seen=set()
            with path.open() as f:
                for row in csv.DictReader(f):
                    b,g,k,r=int(row['batch']),('first','second').index(row['orientation']),int(row['k']),int(row['early_rank'])
                    key=(b,g,k,r)
                    if key in seen:raise RuntimeError('duplicate joint mark')
                    seen.add(key)
                    values[b,g,k]+=[int(row[x]) for x in ('sum_q','sum_e','sum_s_now','sum_qs_now','sum_es_now')]
            if len(seen)!=100*2*(n+1)*3 or not np.array_equal(values,old[f'n{n}']):
                raise RuntimeError(f'N{n}: old batch/K source profile reconstruction failed')
            compressed=path.with_suffix('.csv.gz')
            with path.open('rb') as src,gzip.open(compressed,'wb',compresslevel=6) as dst:dst.write(src.read())
            row={'N':n,'command':cmd,'seconds':time.perf_counter()-t,'stdout':p.stdout,'stderr':p.stderr,
                 'raw_sha256':sha(path),'gzip_sha256':sha(compressed),'old_profile_integer_exact':True,
                 'old_permutations_reobserved':1000000 if n in (260,340) else 100000}
            path.unlink();record['runs'].append(row);save(out/'run.json',record);print(json.dumps(row),flush=True)
        record['status']='acquisition_complete';record['elapsed_seconds']=time.perf_counter()-started
        record['total_old_permutations_reobserved']=sum(r['old_permutations_reobserved'] for r in record['runs'])
        save(out/'run.json',record)
        subprocess.run([os.sys.executable,str(ROOT/'analyze.py')],check=True)
        record['status']='completed';record['total_elapsed_seconds']=time.perf_counter()-started
    except Exception as e:
        record['status']='failed';record['error']=repr(e)
        if isinstance(e,subprocess.CalledProcessError):record['failure_output']={'stdout':e.stdout,'stderr':e.stderr}
        raise
    finally:
        record['finished_unix']=time.time();save(out/'run.json',record)
if __name__=='__main__':main()
