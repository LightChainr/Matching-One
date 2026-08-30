# P250 fresh projective-leg production result

Status: the fresh 10k block confirms mesoscopic charged two-point propagation,
but does not confirm a nonzero cubic vector at the second separation.  The
phase branch remained locked and was not computed.

## Frozen execution

The sample-size grid selected 10k before the run.  XPk2PZ executed seed
`25033433720260930`, counters `[0,10000)`, 50 batches and 16 workers at commit
`da65598`.  Raw-file SHA256 values are

```text
response_10k.json        1e42be903d15b8bd791dc3bcb4d329816ed657720bcdc814c47d38d7edc313de
response_10k.batches.csv 917d78a1a4b1f5d5ade8b6c38d46cf0f435dd6018b82abeffad82533ee9c74dc
```

Two failed launch logs are preserved.  Both failures occurred before replica
generation: the image lacked `/usr/bin/time`, then the old 2k cap correctly
rejected an unauthorized 10k call.  The final runner retained the 2k default
cap and accepted 10k only after exact manifest equality checks.

## Support-first reveal

| d | weakest pair z | pair gate | cubic support | support gate |
|---:|---:|---|---:|---|
| 1 | 15.397 | pass | `24.918/8`, p `0.00160` | pass |
| 2 | 8.608 | pass | `9.633/8`, p `0.29172` | fail |

Both charged two-point tails are substantially more resolved than the frozen
`5 sigma` threshold.  Both cubic covariance matrices are nondegenerate.  The
d1 eight-real cubic vector rejects zero at the frozen one-percent level; d2
does not.

Because the support gate required both rows, the scorer emitted

```text
phase_closure.status = locked_support_gate_failed
phase_closure.computed = false.
```

No phase point, covariance or p-value exists in the result.

## What changed scientifically

The projective-leg operator solved the earlier geometric problem: its charged
two-point mode genuinely propagates beyond contact.  What failed is more
specific.  The fresh d2 three-body charged cumulant did not reproduce the
large 2k smoke score, despite an even stronger pair denominator.  The smoke
alternative was therefore an upward fluctuation or unstable cubic direction,
not a reason to buy more samples.

This creates a cleaner mechanism question than phase fitting.  The fresh block
is compatible with a pair-bearing but rapidly Gaussianizing, charge-selected
mesoscopic mode: two-point structure survives while the neutral 113/122 cubic
vector is unresolved at d2.  That is a hypothesis, not an exact Gaussianity
claim.  Testing it would require a newly frozen higher-cumulant or separation-
decay score; the present phase lock must not be bypassed by a d1-only reveal.

## Decision

Do not extend this block merely to unlock phase closure.  Record the operator
as a successful charged-pair observable and a failed two-separation cubic/OPE
observable.  If P250 continues, promote the pair-nonzero/cubic-null mechanism
as a new question with a fresh observable, rather than reusing the locked
phase target.
