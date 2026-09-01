#!/usr/bin/env python3
"""Verify and recover completed artifacts; leave the frozen analysis untouched."""
from datetime import datetime
import hashlib
import json
from pathlib import Path
import tarfile

ROOT=Path(__file__).resolve().parents[1]
SPECS={
    'HZsCM6':('/private/tmp/p334-prospective-HZ-delivery.tar','f624fc688b1d66de75dfc241505b01cbfc528298277d1be866332c3536681e15'),
    'TV2N0X':('/private/tmp/p334-prospective-TV-delivery.tar','3c0145790864f951c4fe19f8da1c0eafacd5d0b8f918b5dbad06bca3b99a252c'),
}
receipts={}
for host,(name,expected) in SPECS.items():
    data=Path(name);assert hashlib.sha256(data.read_bytes()).hexdigest()==expected
    local=ROOT/'execution'/host;local.mkdir(parents=True,exist_ok=True)
    with tarfile.open(data) as archive:
        for member in archive.getmembers():
            relative=Path(member.name)
            if not member.isfile() or relative.is_absolute() or '..' in relative.parts:
                raise ValueError('unexpected archive member')
            target=ROOT/relative if relative.parts[0]=='production' else local/relative
            if target.exists():raise ValueError('refuse local overwrite')
            target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes(archive.extractfile(member).read())
    checks=json.loads((local/'delivery_checksums.json').read_text())
    for path,digest in checks.items():
        target=ROOT/path if path.startswith('production/') else local/path
        if hashlib.sha256(target.read_bytes()).hexdigest()!=digest:
            raise ValueError('recovered file hash mismatch')
    receipts[host]=json.loads((local/'remote_run.json').read_text())
    receipts[host].update(delivery_tar_sha256=expected,verified_files=len(checks))
runfiles=sorted((ROOT/'production').glob('*.run.json'))
assert len(runfiles)==120
starts=[];seals=[];ends=[]
for p in runfiles:
    r=json.loads(p.read_text())
    assert r['freeze_commit']=='4b3c21b7c8c33a5df7eab7eaa2a9f04af18d1277'
    start,seal,end=[datetime.fromisoformat(r[k]) for k in ('features_started','predictions_sealed','tails_completed')]
    assert start<seal<end
    assert start>datetime.fromisoformat('2026-08-31T19:36:13+08:00')
    starts.append(start);seals.append(seal);ends.append(end)
    stem=p.name.removesuffix('.run.json')
    for suffix,key in [('.prefix.csv.gz','prefix_sha256'),('.census.csv.gz','census_sha256'),('.prediction.csv.gz','prediction_sha256'),('.tails.csv.gz','tail_sha256'),('.sufficient.json','sufficient_sha256')]:
        assert hashlib.sha256((ROOT/'production'/(stem+suffix)).read_bytes()).hexdigest()==r[key]
freeze=json.loads((ROOT/'FREEZE.json').read_text())
for path,digest in freeze['file_sha256'].items():
    assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
result={'status':'all120shards_recovered_and_verified_before_final_score','frozen_commit':'4b3c21b7c8c33a5df7eab7eaa2a9f04af18d1277',
        'commit_time':'2026-08-31T19:36:13+08:00','first_prefix_job_started':min(starts).isoformat(),
        'first_prediction_sealed':min(seals).isoformat(),'last_prediction_sealed':max(seals).isoformat(),
        'first_tail_completion':min(ends).isoformat(),'last_tail_completion':max(ends).isoformat(),
        'all_prediction_files_precede_their_tail_phase':True,'frozen_file_hashes_unchanged':True,'machines':receipts,
        'lifecycle':'scientific production ended; per parent instruction keep HZ/TV Running and tunnels available for P154 handoff'}
(ROOT/'execution'/'RECOVERY.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='machines'},indent=2))
