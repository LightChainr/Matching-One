# Held-out low-rank state transfer on N=145 -> 290

Status: frozen prospective sensitivity analysis before the N145/N290 threshold-rank full-curve reveal. The already-frozen P50 center/slope/root scorer remains the primary chronology. This analysis reuses the same raw block and is therefore nonadditive in the evidence ledger.

## Why freeze this now

After the completed norm-5 reveal, the old H4/H12 angular alias is no longer the main uncertainty. H4 is compatible with the N325/N425 transfer while H12/H8 are rejected, but the child block alone is still compatible with zero. The derivative sector is also narrower: S, D and D-prime obey their leading laws on new geometries, whereas S-prime requires a correction and norm-5 does not distinguish analytic q=2 from rank-2/Jordan.

Two simple scalar explanations have now failed on the same expensive norm-5 curves: a constant rank-gap boundary correction and a single scalar width rescaling of the higher thermal jet. The natural next model is therefore the already-discovered low-rank transfer state, but that model was developed after earlier reveals. N=145 -> 290 is the first clean opportunity to turn it into a held-out prediction.

## State

At the intrinsic matching center define

```text
I_S  = N P4[S]
I_Du = N P4[D'] / Mbar'
T_D  = N^(13/8) P4[D]
T_Su = N^(13/8) P4[S'] / Mbar'
```

The first two coordinates are the matching-even/identity-like block. The last two are the matching-odd/thermal-like block. The discovery model says that the first block is rank one and the thermal block is rank two, with the second thermal state almost invisible in central D but visible in S-prime.

For a doubling N -> 2N, the compact closure predicts zero leading increment for

```text
I_S, I_Du, T_D.
```

The live correction mechanisms differ only in the fourth coordinate:

```text
analytic: T_Su(N) = A + C/N
Jordan:   T_Su(N) = A + B log N.
```

Using the already-committed source fit, without N145/N290 target data,

```text
N145 -> N290
analytic Delta T_Su = +0.1264093880 +/- 0.0080286542
Jordan   Delta T_Su = +0.2555067617 +/- 0.0156669062
```

The separation is large compared with the source-fit uncertainty; actual target power is determined by the full-curve jackknife covariance and is not assumed in advance.

## Score

For each size, recompute the intrinsic center, slope and all P4 projectors inside every delete-one replicate, then transform to the four-component state. The planned N145 and N290 streams are independent, so

```text
Cov(delta_state) = Cov(state_145) + Cov(state_290).
```

For each frozen model, subtract its expected increment vector and add the source-fit increment variance to the T_Su diagonal entry. Report the four-dimensional chi-square and the marginal T_Su signed residual. Do not drop a failed component to improve the model.

This score answers two questions at once:

1. does a compact three-state transfer closure survive a genuinely held-out full curve?;
2. conditional on that closure, does the thermal companion behave more like an ordinary inverse-N eigenmode or a rank-2 logarithmic generalized eigenvector?

## Decision logic

- If both models pass, retain the compact state but call q2/Jordan unresolved; norm-4 dyadic closure is then the cleanest next discriminator.
- If only analytic passes, test its exact dyadic second-difference law before assigning an operator interpretation.
- If only Jordan passes, require an independent multiplier/dyadic closure before calling the state logarithmic.
- If both fail, reject the three-state closure as predictive. Use the held-out residual direction to choose between quotient-arithmetic controls, the self-matching RG tangent, and modulus-shape spectroscopy. Do not add a free scalar exponent first.

## Evidence boundary

The P50 raw block has one primary evidential role. Center/slope/root, this state score, Krawtchouk/Hermite modes, rank-gap statistics, metric-free ratios and multi-u profiles are correlated coordinate systems on the same histograms. They may be compared jointly or used as sensitivity diagnostics, but they are not independent votes.

Executable artifacts:

```text
predictions/p50_matrix_state_transfer_20260829.yaml
scripts/score_p50_matrix_state_transfer.py
tests/test_p50_matrix_state_transfer.py
```
