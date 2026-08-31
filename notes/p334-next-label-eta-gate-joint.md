# P334: the common next-label mode cancels 23–32% of the marginal E innovation

The positive first/completion cross response has a natural dimensionless size.
For the same next-label covariance B in the paired `(F1,F2)` coordinates,

```
eta_next = 2 B12/(B11+B22)
B_AA = (B11+B22)(1+eta_next)
B_EE = (B11+B22)(1-eta_next).
```

Thus eta measures cancellation in E **relative to the sum of the two marginal
next-label variances**. It is not a fraction of all suffix variance or of the
signed E mean. No independence of F1/F2 is presumed: their covariance is the
quantity being measured.

| Same-source magnitude | N325 | N425 |
|---|---:|---:|
| canonical eta_next | 0.231366 +/- 0.007062 | 0.247563 +/- 0.008706 |
| integrated eta_next | 0.305290 +/- 0.018451 | 0.322626 +/- 0.021075 |

The common label therefore cancels about 23–25% of that canonical reference
variance and 31–32% after integration, while increasing A by the same
fraction. This quantifies the newly observed positive Gamma; it is a
post-reveal magnitude of that one direction, not a second independent test.
The data still concern a finite original-prefix population under this shared
label/paired-CRN coupling, not an asymptotic field identity.

## Exact reuse and the gate comparison boundary

Eta was computed solely from the already saved Dnext batch matrices in
`24872eef`, without reading the long CSV again. Each LOO ratio uses the same
nineteen retained original batches. No clipping or eigenvalue selection was
applied.

The fixed 26-column gate readout from `c6ee37a8` is now appended to the complete
Doob covariance factor. Its covariance and cross covariance with Gamma and eta
are retained on the same batches. Those gate columns describe 01/10, whereas
the all-prefix Gamma includes all five R0-containing cells. Binary immediate
gate covariance and the continuous future first/completion conditional
response are distinct estimands; this join imposes no proxy equivalence or
population rate-product factorization.

Both providers read the same fresh continuation source
`e32a85939279b8574278024d647b56d2d1485247`: 640,000 new tails per size on the
original 20,000 prefixes. The new suffix randomness does not make the reused
prefixes or gate/Doob reports independent evidence blocks.

Artifact: `results/p334-next-label-mechanism-joint/score.json`, with the
focused Gamma/eta/gate covariance and two complete low-rank factor files.
For every saved coordinate, the common covariance is exactly `factor.T @
factor`, rank at most 19; no inverse is computed. The primary nine-cell full
matrices and old/new mean comparisons remain in
`results/p334-next-label-doob-quartets/score.json`.

Reproduce only this compact matrix/gate join with
`/Users/lc/python-envs/research-py311/bin/python scripts/p334_next_label_mechanism_joint.py`.
No new MC, DP, raw CSV read, model family or test suite is introduced here.
