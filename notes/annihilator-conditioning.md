# Conditioning of matching annihilators

The matching-annihilator idea is useful only if correction cancellation is not overwhelmed by noise amplification or exponent misspecification. This note records the first arbitrary-precision conditioning checks produced by `scripts/matching_annihilator.py`.

Environment used for the check:

```text
Python 3.13.5
mpmath 1.3.0
100 decimal digits
```

## Two-size leading-term cancellation

Command:

```bash
python scripts/matching_annihilator.py \
  --sizes 15 16 \
  --cancel 13/4
```

Result:

```text
w_15 = -4.285042429789350863560506750489816953660
w_16 =  5.285042429789350863560506750489816953660
L1 noise factor = 9.5700848595787017271
L2 independent-noise factor = 6.8039152044810750696
constraint condition = 7.0206136474206056176e4
```

Up to an overall normalization, this is the known equation

\[
16^{13/4}M_{16}(p)-15^{13/4}M_{15}(p)=0.
\]

Even the two-size accelerator substantially amplifies independent noise. The practical variance can nevertheless be better because the matching functions and neighboring sizes may be strongly correlated; this must be measured, not assumed.

## Naive adjacent three-size cancellation

Canceling both `L^-13/4` and `L^-25/4` on sizes 14, 15, 16 gives

```text
w_14 =   7.683689116355361306960299177243173585683
w_15 = -22.17590527813922058983841071594256889423
w_16 =  15.49221616178385928287811153869939530854
L1 noise factor = 45.351810556278441180
L2 independent-noise factor = 28.121497380922620216
constraint condition = 1.1316518205595397342e9
```

This construction is numerically dangerous. It can turn a formally higher-order estimator into a worse statistical estimator.

## Wider and overdetermined windows

For the same two canceled powers:

| sizes | L1 norm | L2 norm | max abs weight | constraint condition |
|---|---:|---:|---:|---:|
| `14,15,16` | 45.35 | 28.12 | 22.18 | `1.13e9` |
| `12,14,16` | 10.69 | 6.90 | 4.85 | `1.55e8` |
| `8,12,16` | 2.82 | 2.07 | 1.86 | `8.81e6` |
| `12,13,14,15,16` | 11.09 | 5.52 | 3.74 | `1.59e8` |
| `8,10,12,14,16` | 2.42 | 1.35 | 0.97 | `6.26e6` |
| `8,9,...,16` | 2.72 | 1.02 | 0.56 | `6.50e6` |

The apparent lesson is not simply “use small sizes.” Wider windows reduce noise amplification but expose the estimator to larger unmodeled corrections. The correct optimization problem is therefore covariance- and bias-aware:

\[
\min_w\; w^T\Sigma w + \lambda\,B(w)^2
\]

subject to the declared annihilation constraints. Here `Sigma` is the measured covariance of the matching observables and `B(w)` is a conservative model for uncanceled corrections.

## Research consequence

Higher matching annihilators should be treated as constrained experimental design, not mechanical Richardson extrapolation. The next implementation must support:

1. generalized least-variance weights using an empirical covariance matrix;
2. explicit penalties for omitted powers and exponent uncertainty;
3. training-window selection without viewing the largest sizes;
4. noise-amplification and condition-number rejection thresholds;
5. out-of-sample comparison against the ordinary and two-size roots.

A useful result may be a robust estimator of lower formal order rather than the highest algebraic cancellation available.
