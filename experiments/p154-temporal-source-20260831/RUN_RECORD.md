# Run and retrieval record

- Frozen reviewed source: `7da1eeb0e51cf430987dbf204d23713c2ab5a46c`.
- Original RNG/geometry/UF source: `bfab0330f5f56ca4d746b45d737f1607e3d229a0`; vendor bytes and source-result hashes are in SOURCES.json.
- Assigned machine: DevEnvC_ZyTrST, instance `f415a4bcbd9a438b85f5f29e4a507ea4`, 16vCPU/32GiB ARM. Only this machine was used by this experiment.
- Root established the authorized session and restored authentication before handoff. This subtask did not reset a key or start a second tunnel. Canonical port10023 and the existing managed tunnel PID77144 were reused.
- GCC10.3.1 was installed with `yum -y install gcc-c++`. Python used `/workspace/matching-one-p398-rate-20260831/.venv/bin/python`, Python3.9.9/NumPy1.26.4/SciPy1.13.1.
- Remote package: `/workspace/p154-temporal-source-20260831`.
- Actual command: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=16 /workspace/matching-one-p398-rate-20260831/.venv/bin/python -u run.py > run.log 2>&1`.
- The driver acquired each old permutation once. Totals: 2,400,000 old permutations, 4,800,000 paired geometry observations, zero new counters. Every all-K path contributes dependent rows; these are not additional independent samples.
- Acquisition/compilation/compression:41.643954s; score7.248716s; driver50.200440s. Per-N driver spans: .741750/.981858/1.474418/1.927060/12.871352/17.360297s, in increasing N order.

## Frozen execution hashes

| Object | SHA256 |
|---|---|
| CONTRACT.json | `e2c15461d91250faaad064e9862157d84b1b9c32069cb271fad7501cb3a466be` |
| replay.cpp | `cd6d9386db6e23ee04381644cb25b94222f91e2cbdb9b8fad26ccd8195a21a21` |
| run.py | `88bff4931e2eeb47a39ee512a47662d4f2d27b649c65eb0cd1f29ae1c7fe6861` |
| analyze.py | `d435c81268766ed656c1a38528a9432a93720c606b063802fa2d3bc5910c3722` |
| results/latest.json | `93f753a9a447773a1367fb728304e4a1028d83fc4d3c685bce506ded6280534c` |
| Retrieved results tar.gz | `c7531e1b06614d8ee5d9b8e9e368d27df3cbd2fb04f13d358f582ff4de5b4b65` |

The retrieved tar includes results and run.log. Its archive hash, result JSON hash and all six compressed raw hashes matched the remote receipt. The executed replay/scorer/contract were not changed afterward. README, REPORT and this record were written locally around/after the completed execution.

The separate local `verify_moving_root.py` check used the managed Python3.11 research environment on the retrieved data. It added no configurations; its48 root solves are explicitly an auxiliary derivative check, recorded in results/moving_root_verification.json. All12 source/size comparisons passed, maximum absolute derivative discrepancy1.8223e−7.

## Machine closure

After driver completion and hash-verified retrieval, authenticated `ps -eo pid,comm,etime,pcpu,pmem` showed no experiment process. The official command `hdspace devenv stop --instance-id=f415a4bcbd9a438b85f5f29e4a507ea4` returned Stopping. A subsequent official `devenv view` confirmed **Ready**. The instance-specific helper stop reported that its managed tunnel was no longer running. No unrelated VM or tunnel was stopped.

No source-team worktree was modified. No commit, push, GitHub mutation or other-team status message was made by this subtask.
