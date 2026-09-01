#!/usr/bin/env python3
"""Record one execution of the already frozen shard runner; no analysis changes."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import socket
import subprocess
import time

p=argparse.ArgumentParser()
p.add_argument('--python',required=True)
p.add_argument('--shards',required=True)
args=p.parse_args()
root=Path.cwd()
if (root/'remote_run.json').exists():
    raise SystemExit('a previous execution receipt already exists')
def read(path):
    try:return Path(path).read_text().strip()
    except OSError:return None
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
cmd=[args.python,'run_shards.py','--freeze-commit','4b3c21b7c8c33a5df7eab7eaa2a9f04af18d1277','--shards',args.shards,'--workers','14']
r={'hostname':socket.gethostname(),'architecture':platform.machine(),'command':cmd,
   'started_utc':datetime.now(timezone.utc).isoformat(),'freeze_commit':'4b3c21b7c8c33a5df7eab7eaa2a9f04af18d1277',
   'quota':read('/sys/fs/cgroup/cpu/cpu.cfs_quota_us'),'period':read('/sys/fs/cgroup/cpu/cpu.cfs_period_us'),
   'memory_limit':read('/sys/fs/cgroup/memory/memory.limit_in_bytes'),
   'binary_sha256':hashlib.sha256((root/'prospective').read_bytes()).hexdigest(),
   'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0]}
start=time.monotonic();child=subprocess.Popen(cmd,env=env);r['pid']=child.pid
(root/'remote_run.json').write_text(json.dumps(r,indent=2)+'\n')
peak=0
while True:
    usage=read('/sys/fs/cgroup/memory/memory.usage_in_bytes')
    if usage:peak=max(peak,int(usage))
    try:code=child.wait(timeout=1);break
    except subprocess.TimeoutExpired:pass
r.update(exit_code=code,finished_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-start,
         sampled_container_memory_peak_bytes=peak,child_maxrss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
         memory_measurement='container current usage sampled each second; includes container services; child maxrss is not sum over concurrent processes')
(root/'remote_run.json').write_text(json.dumps(r,indent=2)+'\n')
raise SystemExit(code)
