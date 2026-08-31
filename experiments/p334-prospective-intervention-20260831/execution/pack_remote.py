#!/usr/bin/env python3
"""Package completed raw shards and receipts without viewing effect estimates."""
import hashlib
import json
from pathlib import Path
import socket
import tarfile

root=Path.cwd()
r=json.loads((root/'remote_run.json').read_text())
assert r['exit_code']==0
assert len(list((root/'production').glob('*.sufficient.json')))==60
paths=sorted(p for p in (root/'production').rglob('*') if p.is_file())
paths += [root/n for n in ('remote_run.json','run.log','input_hash_check.txt')]
checks={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
(root/'delivery_checksums.json').write_text(json.dumps(checks,indent=2)+'\n')
target=root/'delivery.tar'
if target.exists():raise SystemExit('refuse overwrite delivery')
with tarfile.open(target,'w') as out:
    for p in paths+[root/'delivery_checksums.json']:
        out.add(p,arcname=str(p.relative_to(root)),recursive=False)
print(json.dumps({'hostname':socket.gethostname(),'file_count':len(paths)+1,'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}))
