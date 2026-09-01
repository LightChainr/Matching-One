#!/usr/bin/env python3
"""Generate one fixed cell00 extension budget; first batch is retained as preflight."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=16)
    args=parser.parse_args()
    root=Path(__file__).parent
    target=root/'extension'
    target.mkdir(exist_ok=False)
    started=time.time()
    for f in json.loads((root/'inputs/prefix_archive/manifest.json').read_text())['files']:
        assert hashlib.sha256((root/'inputs/prefix_archive'/f['local_path']).read_bytes()).hexdigest()==f['sha256']
    def run(task):
        n,b=task
        command=[str(root/'cell00_extension'),str(n),str(b),str(target)]
        p=subprocess.run(command,cwd=root,capture_output=True,text=True)
        (target/f'N{n}.batch{b:02}.log').write_text(p.stdout+p.stderr)
        if p.returncode: raise RuntimeError(f'Extension failed {n}/{b}: {p.stderr}')
        return json.loads((target/f'N{n}.batch{b:02}.metadata.json').read_text())
    first=run((325,0))
    seen={}
    with gzip.open(target/'N325.batch00.csv.gz','rt') as f:
        for r in csv.DictReader(f):
            key=(int(r['counter']),int(r['quartet']),int(r['group']))
            assert 8<=key[1]<=71 and int(r['first_rank'])==int(r['second_rank'])==0
            val=(int(r['next_label']),int(r['first_next_rank']),int(r['second_next_rank']),int(r['first_e']),int(r['first_c']),int(r['second_e']),int(r['second_c']))
            if key in seen: assert seen[key]==val
            seen[key]=val
    assert len(seen)==first['prefixes']*64*2
    assert first['new_tail_paths']==first['prefixes']*64*4
    # Domain proof uses actual bounds; splitmix64 and xor are bijections on64 bits.
    old_max=((19999*8+7)*8+1*4+2)
    new_max=((19999*64+63)*8+1*4+2)
    assert old_max < 2**31 and new_max < 2**31
    (target/'preflight.json').write_text(json.dumps({'first_batch':first,'label_group_keys':len(seen),'old_max_local_id':old_max,'new_low_max':new_max,
        'new_bit':31,'stream_domains_disjoint':True,'replica_label_and_immediate_contact_match':True,'preflight_is_part_of_fixed_budget':True},indent=2)+'\n')
    tasks=[(n,b) for n in (325,425) for b in range(20) if (n,b)!=(325,0)]
    with ThreadPoolExecutor(max_workers=args.workers) as pool: meta=[first]+list(pool.map(run,tasks))
    assert sum(x['prefixes'] for x in meta)==3053
    assert sum(x['new_tail_paths'] for x in meta)==781568
    receipt={'host':platform.node(),'workers':args.workers,'started_unix':started,'finished_unix':time.time(),'elapsed_seconds':time.time()-started,
       'prefixes':3053,'new_prefixes':0,'new_quartets_per_prefix':64,'new_tail_paths':781568,'same_old_population':True,
       'producer_sha256':hashlib.sha256((root/'src/cell00_extension.cpp').read_bytes()).hexdigest(),'executable_sha256':hashlib.sha256((root/'cell00_extension').read_bytes()).hexdigest(),
       'files':[{ 'name':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(target.glob('*.csv.gz'))]}
    (target/'run_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='files'}),flush=True)

if __name__=='__main__':main()
