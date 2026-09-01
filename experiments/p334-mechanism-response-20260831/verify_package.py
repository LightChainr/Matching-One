#!/usr/bin/env python3
"""Local integrity and retained covariance checks; no scientific resampling."""
import hashlib
import json
import math
from pathlib import Path

root=Path(__file__).parent
checked=0
for folder in ('inputs','census','inputs/prefix_archive'):
    manifest=json.loads((root/folder/'manifest.json').read_text())
    for f in manifest['files']:
        assert hashlib.sha256((root/folder/f['local_path']).read_bytes()).hexdigest()==f['sha256']
        checked+=1
extension=json.loads((root/'extension/run_receipt.json').read_text())
for f in extension['files']:
    assert hashlib.sha256((root/'extension'/f['name']).read_bytes()).hexdigest()==f['sha256']
    checked+=1
assert extension['new_tail_paths']==781568 and extension['prefixes']==3053
meta=[json.loads(p.read_text()) for p in (root/'extension').glob('*.metadata.json')]
assert len(meta)==40 and sum(x['new_tail_paths'] for x in meta)==781568
assert {n:sum(x['prefixes'] for x in meta if x['N']==n) for n in (325,425)}=={325:1502,425:1551}
score_summaries={}
for folder in ('results','results-exact-score','results-extension'):
    score=json.loads((root/folder/'score.json').read_text())
    max_factor_error=0.
    for n,s in score['sizes'].items():
        assert s['batch_ids']==list(range(20))
        assert len(s['labels'])==len(s['estimate'])==len(s['se'])
        assert len(s['LOO'])==len(s['factor'])==20
        for j in range(len(s['labels'])):
            avg=math.fsum(row[j] for row in s['LOO'])/20
            for row,factor in zip(s['LOO'],s['factor']):
                target=math.sqrt(19/20)*(row[j]-avg)
                max_factor_error=max(max_factor_error,abs(target-factor[j]))
            se=math.sqrt(math.fsum(row[j]**2 for row in s['factor']))
            assert math.isclose(se,s['se'][j],rel_tol=1e-12,abs_tol=1e-25)
    assert max_factor_error<1e-12
    score_summaries[folder]={'factor_check_max_absolute_difference':max_factor_error,'score_sha256':hashlib.sha256((root/folder/'score.json').read_bytes()).hexdigest()}
final=json.loads((root/'results-extension/score.json').read_text())
qa={'input_files_hash_verified':checked,'new_fork_files':40,'new_tail_paths':781568,'new_prefixes':0,'all_covariance_factor_checks':score_summaries,
    'fast_fourth_order_algebra_check':final['algebra_vs_direct_four_subset_error'],
    'old8_vs_enumeration_max_error':{n:x['fast_vs_old8_enumeration_max_error'] for n,x in final['sizes'].items()},
    'no_clipping_of_signed_det_square':True,'baseline_and_extension_same_original_population':True}
(root/'QA.json').write_text(json.dumps(qa,indent=2)+'\n')
print(json.dumps(qa,indent=2))
