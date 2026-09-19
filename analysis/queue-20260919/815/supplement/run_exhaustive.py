"""Compile and run the bounded #815 census; never substitutes for #775."""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import platform
import random
import subprocess
import tempfile
from time import perf_counter

ROOT = Path(__file__).resolve().parent
REFERENCE_SHA = '9fc042535aa0707651a5e38a2e2121669d99d59a'


def main():
    start = perf_counter()
    data = (ROOT/'reference773.py').read_bytes()
    assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest() == REFERENCE_SHA
    spec = importlib.util.spec_from_file_location('reference773', ROOT/'reference773.py')
    ref = importlib.util.module_from_spec(spec); spec.loader.exec_module(ref)
    with tempfile.TemporaryDirectory(prefix='issue815-') as temp:
        executable = Path(temp)/'check'
        command = ['g++','-std=c++17','-O3','-Wall','-Wextra','-pedantic',
                   str(ROOT/'exhaustive_check.cpp'),'-o',str(executable)]
        subprocess.run(command,check=True,capture_output=True,text=True)
        compiler = subprocess.check_output(['g++','--version'],text=True).splitlines()[0]
        oracle = []
        for L in range(1,6):
            if L <= 4:
                masks = list(range(1 << (L*L))); sampling = 'all masks'
            else:
                rng = random.Random(81520260919)
                masks = [0,(1 << 25)-1]+rng.sample(range(1 << 25),2000)
                sampling = '2000 fixed-seed random masks plus empty/full'
            p = subprocess.run([str(executable),'--ranks',str(L)],
                               input='\n'.join(map(str,masks))+'\n',text=True,capture_output=True,check=True)
            actual = list(map(int,p.stdout.split()))
            expected = [ref.torus_rank(m,L) for m in masks]
            assert actual == expected, L
            oracle.append({'L':L,'masks':len(masks),'selection':sampling,'mismatches':0})
        systems = []
        for L in range(1,6):
            p = subprocess.run([str(executable),'--census',str(L)],text=True,capture_output=True,check=True)
            systems.append(json.loads(p.stdout))
        classification_controls = []
        for L in (3,4):
            result = ref.classify(L); candidate = systems[L-1]
            assert result['jump2_total'] == candidate['jump2']
            for key in ('T3','T4plus','Rsplit'):
                assert result['type_counts'].get(key,0) == candidate[key]
            classification_controls.append({'L':L,'jump2':result['jump2_total'],
                                             'type_counts':result['type_counts'],'mismatches':0})
        guard = subprocess.run([str(executable),'--census','6'],capture_output=True,text=True)
        assert guard.returncode == 2
    output = {'issue':815,'reference_git_blob':REFERENCE_SHA,
              'rank775_kernel_retrieved':False,'integer_algorithms':['integer lift DFS','integer unit-capacity max flow'],
              'oracle_checks':oracle,'reference_classification_rerun':classification_controls,
              'systems':systems,'L3_to_L5_all_jump2_checked':sum(s['jump2'] for s in systems if s['L']>=3),
              'L3_L4_arm_contract':'thin annuli; disjoint single-vertex paths allowed, nonempty ordered endpoints required',
              'L1_scope':'rank census only; periodic self-loop quotient has no outside-site arm event',
              'L2_scope':'cover-germ check only; the radius-one box is not injective in the quotient',
              'R1_limit_claim':'no square-site exponent or continuum b(lambda)=0 upgrade',
              'L6_guard':'rejected as required','compiler':compiler,'platform':platform.platform(),
              'available_cpus':os.cpu_count(),'total_wall_seconds':round(perf_counter()-start,6)}
    (ROOT/'exhaustive-result.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(output,ensure_ascii=False,indent=2))

if __name__ == '__main__':
    main()
