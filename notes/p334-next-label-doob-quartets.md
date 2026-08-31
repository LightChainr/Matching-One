# P334: a frozen next-label covariance readout for new suffix production

This scorer is prepared before the final fresh-suffix archive is read. It
does not launch production or consume partial files. The population is the
original N325/N425 set of 20,000 paired prefixes per size; the randomness being
added is their continuation, not a new prefix sample.

For one prefix, a quartet samples independent next labels U,V (possibly equal)
and two independent suffixes a,b under each label. Each complete tail is
shared by both orientations. X is the paired H4-normalized vector
`[F1,F2,A,E]` at p_ref and after integration, followed by raw-rank `K1,K2,W`.
Write `a=X_Ua-X_Va`, `b=X_Ub-X_Vb`. Then

```
Vtot   = (a a^T+b b^T)/4
Dnext  = (a b^T+b a^T)/4
Vafter = (a-b)(a-b)^T/4
Vtot   = Dnext+Vafter.
```

The target averages are total conditional suffix covariance, next-label Doob
covariance and after-next conditional covariance. Dnext is an unbiased signed
cross-product estimate: finite-sample negative entries/eigenvalues are not
clipped or converted into state counts. Exact/degenerate zero denominators are
reported as not scoreable, never silently replaced.

Eight quartets mean sixteen independent next-label draws, each with two
suffix replicas. The conditional covariance of the 32-tail mean is therefore
`Dnext/16+Vafter/32`, **not** `Vtot/32`. This is the precise resource question:
whether next-label randomness or the remaining suffix limits conditional
averaging.

The removed covariance is `(15/16)Dnext+(31/32)Vafter`. In contrast, an old
baseline minus a fresh finite-fork mean contains both sources of suffix noise;
its outer product is not removed noise. The single signed first/completion
cross readout is `Gamma=(Dnext_AA-Dnext_EE)/4`, with the same complete covariance.

Every matrix and mean is retained for all nine original joint-rank cells and
their six transpose groups, with the full-population denominator. The main
allocation readouts are the 01+10 shares in K1, canonical E and integrated E.
The fresh 32-tail mean is compared to old complete baseline and safe estimates
from `bb79fd47` using their common twenty original batches. The low-rank
covariance factor retains all cross-group and matrix-entry errors without a
large covariance inversion or an independence approximation.

Frozen interface (run only after the root supplies the completed source SHA):

```
/Users/lc/python-envs/research-py311/bin/python scripts/p334_next_label_doob_quartets.py \
  --source-commit FINAL_SHA --source-directory FINAL_DIRECTORY
```

The committed directory must contain exactly the forty named complete batch
files `N325.batch00.csv.gz` through `N425.batch19.csv.gz`; a `batches/`
subdirectory is allowed. The 16-column header and fixed quartet semantics are
recorded in `analysis/p334_next_label_doob_quartets.json`. No scientific result
is asserted by this interface-only commit.
