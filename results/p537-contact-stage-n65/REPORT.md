# P537 held-out N65 contact-stage gate

Status: `CONTACT_FUSION_COMPLETION_TRANSMITS`

The frozen 20,000,000-sample N65 production reproduces the N25 sign rotation
in a new population and in the full pooled-root original-U score.  With rows
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
gate.  The N25 local phenomenon therefore survives as a global N65 response:
single-arm contact remains negative at both births, while double-arm contact
changes sign only at topological completion.  A scalar contact loading or a
separable `contact x birth-stage` law cannot produce this rank-two tensor.

The run used the independent 100-batch P45 N65 baseline, four counter-keyed
shards, 100 new batches, seed `20260901537`, and no top-up.  The raw TSV shards
are intentionally represented by `SHA256SUMS` rather than committed.

## Two-scale operator fingerprint

The preregistered positive exposures show that the exceptional entry-double
cell is not suppressed because the event becomes unusually rare.  Comparing
the exact N25 tensor with N65 gives the zero-parameter power matrix

```text
N^-3      N^-29/8
N^-3      N^-3
```

Its N25-to-N65 prediction has joint `Q=0.6364` on 4 df (`p=0.95893`).  A
common `N^-3` law is rejected (`Q=11.706`, `p=0.01967`), as is an arbitrary
common rescaling of the N25 matrix (`p=0.01591`).  The entry-double exposure
decays as `N^-5.0095`, while its conditional signed density grows only as
`N^1.3879`, giving the net `N^-3.62157 ~= N^-29/8`.  Relative to the three
`N^-3` cells, the missing `N^5/8=L^5/4` conditional amplification is the
thermal-dimension fingerprint.

The projective cross-ratio, invariant under every row and column rescaling,
obeys

```text
N^(-5/8) * (-chi_N) = 0.1528319  at N25
                    = 0.1505870 +/- 0.0485956 at N65.
```

The finite tensor is therefore naturally a noncommuting contact-fusion and
birth-completion transfer, with an asymptotically triangular candidate block,
not a scalar contact coefficient.  Complete numbers are in
`scale-fingerprint.json`; the operator interpretation is developed in
`notes/p537-contact-completion-commutator.md`.
