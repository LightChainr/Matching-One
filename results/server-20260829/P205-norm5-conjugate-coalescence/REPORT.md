# P205 same-parent norm-5 conjugate coalescence — 10M prospective score

## Result

The frozen H4 affine law **survives** the first same-parent conjugation and
cyclic-to-noncyclic quotient control:

| frozen model | z at N325 | z at N425 | joint chi-square / 2 | survival |
|---|---:|---:|---:|---:|
| H4 (primary) | +1.1794 | -1.1541 | 2.72284 | 0.25630 |
| H12 | -0.0321 | +0.8431 | 0.71189 | 0.70051 |
| H8 | +0.4997 | -0.4693 | 0.46995 | 0.79059 |

This removes the grossest quotient-class failure mode: changing the new C node
from a cyclic quotient to Smith `(5,65)` / `(5,85)` does not break the frozen
one-harmonic affine response at 10M samples per pair.  It does **not** identify
H4 uniquely.  H8 is the smallest-chi-square fixed model, H12 is nearly as good,
and the H4 displacement from H8 is only `Delta chi-square = 2.25290`.  The bold
mechanism reading is therefore: a low-dimensional angular transfer survives
the noncyclic coalescence experiment, while the harmonic label remains aliased
at this precision.

## Frozen score

The fixed observable is

```text
M_G(p_ref) = R_G,cross(p_ref) - R_matching,cross(1-p_ref)
p_ref = 0.59274605079
```

The raw point estimates were:

| N | M_C | M_A | M_B |
|---|---:|---:|---:|
| 325 | 1.8448614e-4 | -3.7359881e-5 | 2.0440622e-4 |
| 425 | -1.3085175e-4 | -1.9667171e-4 | 5.0227889e-5 |

The primary residuals used the preregistered integer nulls without fitting:

```text
N325:  5 M_C - 11 M_A +  6 M_B = 0
N425: 20 M_C + 13 M_A - 33 M_B = 0
```

H12 and H8 were then scored in the frozen order using the exact YAML weights.
At each size the scorer reconstructs the full 3x3 C/A/B delete-one covariance
under the common priority field and propagates it through the affine null.  The
two sizes have independent seeds, so their residual covariance is block
diagonal.  No radial exponent, fitted center, source amplitude, target-selected
harmonic or quotient offset enters the score.

## Common-field and provenance result

The duplicated C-first histogram rows and C-first moment rows are byte-identical
inside both sizes at both preflight and production scale.  All remote and local
SHA-256 values agree after intake.  Production used:

| size | Huawei DevEnv | jobs | seed | counter interval | wall seconds |
|---|---|---|---:|---|---:|
| N325 | DevEnvC_ZyTrST | C-A / C-B, 8 threads each | 2026105501 | `[9300000000,9310000000)` | 84.93 / 83.55 |
| N425 | DevEnvC_HZsCM6 | C-A / C-B, 8 threads each | 2026105502 | `[9300000000,9310000000)` | 110.70 / 109.06 |

Both machines used commit
`1668bddd9e8b10cbe5c57add79eb61f74861fb4b`, GCC 10.3.1, and the identical
ARM64 binary SHA-256
`ee9010f524935099ba22f1820fddc05a79dd309d98784dfd5ba7da28129b6856`.

## Information-value consequence

The 10M block answers the quotient-control question but not the harmonic
selection question.  A post-reveal variance projection, conditional on the
current H4 residual means persisting, places the two-degree-of-freedom 5%
crossing at about **22.0M total samples per pair**.  Thus a clean next increment
would be 20M new replicas per pair (30M total), using disjoint counters and the
same frozen residuals.  This is a planning projection, not part of the
prospective 10M evidence and not a reason to relabel the present H4 pass.

## Reproduction

```bash
python3 scripts/score_p205_norm5_conjugate_coalescence.py \
  --pair 325:A:results/server-20260829/P205-norm5-conjugate-coalescence/raw/n325_C_A_10m.hist.csv:results/server-20260829/P205-norm5-conjugate-coalescence/raw/n325_C_A_10m.moments.csv:results/server-20260829/P205-norm5-conjugate-coalescence/raw/n325_C_A_10m.metadata.json \
  --pair 325:B:results/server-20260829/P205-norm5-conjugate-coalescence/raw/n325_C_B_10m.hist.csv:results/server-20260829/P205-norm5-conjugate-coalescence/raw/n325_C_B_10m.moments.csv:results/server-20260829/P205-norm5-conjugate-coalescence/raw/n325_C_B_10m.metadata.json \
  --pair 425:A:results/server-20260829/P205-norm5-conjugate-coalescence/raw/n425_C_A_10m.hist.csv:results/server-20260829/P205-norm5-conjugate-coalescence/raw/n425_C_A_10m.moments.csv:results/server-20260829/P205-norm5-conjugate-coalescence/raw/n425_C_A_10m.metadata.json \
  --pair 425:B:results/server-20260829/P205-norm5-conjugate-coalescence/raw/n425_C_B_10m.hist.csv:results/server-20260829/P205-norm5-conjugate-coalescence/raw/n425_C_B_10m.moments.csv:results/server-20260829/P205-norm5-conjugate-coalescence/raw/n425_C_B_10m.metadata.json \
  --dps 70 \
  --output results/server-20260829/P205-norm5-conjugate-coalescence/analysis/score.json
```
