# P226 norm-5 chiral Hecke production result

Status: scored under the pre-reveal manifest and joint 4D GLS contract.

Run: `N=325`, `p=0.592746050790`, radius 1, 200,000 samples, 100
batches, 16 Huawei aarch64 workers. The primary vector is one correlated
four-vector `(plus_re,plus_im,minus_re,minus_im)`; its components are not
counted as independent votes.

## Result

```text
primary = (2.7243, -0.9509, 1.1561, -1.3039)
observed handed-ratio phase = 29.20 +/- 49.88 degrees
```

The frozen model ranking is

```text
model  chi2/df   p(df=2)   relative likelihood   normalized weight
H8     1.074/2   0.585     1.000                 0.712
H12    3.539/2   0.170     0.292                 0.208
H4     5.430/2   0.066     0.113                 0.081
```

H8 is the clear winner of this acquisition, but 200,000 samples do not make
the comparison decisive: neither H12 nor H4 is rejected at a strict 5% gate.
The correct conclusion is a 71/21/8 likelihood ordering, not a discovered
spin-8 field.

The scientific tension is nevertheless useful. The ordinary norm-5 global
transfer favored H4, whereas this charged same-parent handed coordinate favors
H8. The two observables therefore should not be assumed to probe one scalar
angular amplitude. A charged/defect sector mixture or a different operator
projection is now a live positive explanation, rather than merely extra fit
noise.

The true reflected-pair conjugacy relation is an exact configurationwise null
and passed identically. The same-parent plus/minus pair is not a reflection
null; its nonzero phase is the intended signal.
