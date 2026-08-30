# Held-out N145 reveal of the natural A-current coordinate

Status: fresh third scale, produced only after commit `3d238af` froze the
observable, primitive Gaussian pair, sample design, seed and four numeric
targets.  No N145 quantity changed a target or normalization.

## Frozen experiment

The coordinate remains

```text
K_A = p(1-p) Jminus_A/W_A = d_eta log W_A,
Delta_K_A = K_A(second)-K_A(first).
```

The primitive equal-area N145 pair is `12+i` versus `9+8i`.  Both quotient
representations have norm 145, gcd one and Smith invariants `(1,145)`; they are
not Gaussian associates or conjugates, and reduction modulo 3 is
nondegenerate.  Production used 2.4M samples/shape, 40 aligned batches, seed
`202608337145` and replica offset zero.

Before production the N145 targets were frozen as:

| target | value | role |
|---|---:|---|
| zero | 0 | primary |
| N85-fitted scale-neutral | 0.0255926129 | primary |
| N85-fitted project H4 | 0.0107447763 | primary |
| N65-to-N85 descriptive-ratio continuation | 0.00350404814 | secondary |

The secondary target was explicitly registered as a diagnostic, not a newly
promoted exponent model.

## Held-out result

At the unchanged `p_ref=0.592746050790`:

```text
N145 first  (12+i): W_A=0.337792502452,
                     Jminus_A=-0.015005967019,
                     K_A=-0.010723781470;

N145 second (9+8i): W_A=0.337209109436,
                     Jminus_A=+0.013393047720,
                     K_A=+0.009587692386.

Delta_K_A(N145)=0.020311473856 +/- 0.003136452035,
z=6.476.
```

The third scale remains clearly nonzero.  Its ratio to the already revealed
N85 point is `0.79365`: attenuation from N85 to N145 is weak, not a continuation
of the fast N65-to-N85 ratio.

## Frozen target score

| target | residual | measurement quadratic | predictive quadratic |
|---|---:|---:|---:|
| zero | +0.0203115 | 41.938 | 41.938 |
| N85-fitted scale-neutral | -0.00528114 | 2.835 | 0.384 |
| N85-fitted project H4 | +0.00956670 | 9.304 | 4.376 |
| secondary fast continuation | +0.0168074 | 28.716 | 11.504 |

The predictive column adds only the target-fit variance declared before N145.
Scale-neutral propagation is the closest primary frozen target.  The secondary
fast continuation is disfavored, and the pure N85-anchored H4 transfer also
underpredicts N145.

The result answers the preregistered question asymmetrically: the fast N85
attenuation is not stable at the third scale.  The live reading is correction
curvature or state/geometry mixing.  It is not evidence for a new fitted
exponent, and one point does not establish scale-neutral asymptotics.

## Dependency and exact gates

- the N145 seed/counter block is independent of the N65 and N85 archives;
- every production field matches the `3d238af` freeze;
- `dW_A/dp=J_A,birth-J_A,exit` closes to `1.36e-14`;
- the Huawei engine self-test passed and stderr is empty;
- raw files are retained with SHA256 checksums;
- observed SE `0.00314` is 1.371 times the N85 variance projection.

## Reproduction

```bash
python3 scripts/score_natural_current_third_scale.py \
  --preregistration analysis/p337_natural_current_third_scale_preregistration.json \
  --births results/server-20260830/P337-natural-current-third-scale-N145/raw/n145_2p4m.births.csv \
  --metadata results/server-20260830/P337-natural-current-third-scale-N145/raw/n145_2p4m.metadata.json \
  --json results/server-20260830/P337-natural-current-third-scale-N145/score.json \
  --markdown results/server-20260830/P337-natural-current-third-scale-N145/REPORT.md

python3 -m unittest discover -s tests \
  -p 'test_score_natural_current_third_scale.py'
```
