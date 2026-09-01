# P537 held-out N65 contact-stage gate

Status: `CONTACT_FUSION_COMPLETION_TRANSMITS`

The frozen 20,000,000-sample N65 production reproduces the N25 sign rotation
in a new population and in the canonical selected-carrier allocation of the
pooled-root original-U score.  With rows
`0->1, 1->2` and columns `single contact, double contact`, the primary tensor is

```text
-1.57857290796e-7  -9.22765780050e-8
-3.09129846134e-7  +3.69680208602e-7
```

Thus the four preregistered signs are `[-,-,-,+]`.  The nonfactorization score
is

```text
Delta_cs = -8.68821605512e-14
SE       =  1.41397768630e-14
95% CI   = [-1.14596123203e-13, -5.91681978997e-14]
theta_cs = -1.0
```

The strict negative interval and the frozen sign pattern pass the prospective
gate.  The N25 local phenomenon therefore survives in the same frozen N65
allocation:
single-arm contact remains negative at both births, while double-arm contact
changes sign only at topological completion.  A scalar contact loading or a
separable `contact x birth-stage` law cannot produce this rank-two tensor.

This is a rejection of scalar/separable structure for that allocation.  It is
not yet a coordinate-invariant operator statement: selected cells generally
move under a common thermal gauge change even though their complete sum does
not.  Also, `theta=-1` is forced by the open sign cone `[-,-,-,+]`; its zero
jackknife SE is not an infinite-precision effect size.

The run used the independent 100-batch P45 N65 baseline, four counter-keyed
shards, 100 new batches, seed `20260901537`, and no top-up.  The raw TSV shards
are intentionally represented by `SHA256SUMS` rather than committed.

[`AUDIT.json`](AUDIT.json) completes the preregistered positive-exposure and
full two-group covariance retention without changing the frozen primary
result.  It also reports the selected-carrier total and labels conditional
density and marginal cell diagnostics as exploratory.
