# Finite-volume Russo bridge for matching slopes

Status: exact finite-system identity and regression layer for Issue #100. This
note does not identify a continuum operator.

## Exact identity

Let `R_G(p)` be any monotone wrapping event for occupied sites on a finite
periodic quotient, and let `R_hat(q)` be the same declared topology channel on
the matching graph. Define

```text
M(p) = R_G(p) - R_hat(1-p).
```

For a vertex `v`, let `Piv_G(v)` be the event that changing only `v` from
vacant to occupied changes the primal wrapping indicator from zero to one.
Russo's formula on the finite product space gives

```text
d R_G(p)/dp = sum_v P_p(Piv_G(v)).
```

Writing `q=1-p` and differentiating the second term gives a second plus sign:

```text
M'(p)
  = sum_v P_p(Piv_G(v))
  + sum_v P_(1-p)(Piv_hat(v)).
```

The plus sign is an exact chain-rule fact, not a scaling assumption. It makes
the finite matching slope a sum of two nonnegative total pivotal masses.

## Independent regression paths

`scripts/exact_pivotal_russo.py` compares three calculations:

1. differentiate the exact channel-specific matching polynomial;
2. enumerate every environment of every candidate pivotal site on the primal
   and matching graphs;
3. for `cross`, reconstruct the derivative from every threshold-rank
   permutation using the frozen `K_minus/K_plus` convention.

Tests cover tiny axis, diamond, and primitive Gaussian tori. They also check
the direct subset identity for `either` and monotonicity for every channel
currently exposed by the homology classifier.

## Interpretation boundary

The exact result justifies calling `M'(p)` a total pivotal mass. It does not by
itself prove any of the following:

- that finite-size corrections are governed by one particular LCFT module;
- that a matching-even or matching-odd angular derivative has a scalar parity
  assignment in the continuum;
- that a local spin-4 moment of the pivotal measure has already been defined.

The observed leading law `M' ~ L^(3/4)=N^(3/8)` is consistent with the known
near-critical pivotal exponent, while the repository's resolved finite-size
drift remains part of the empirical score. A continuum upgrade still requires
an explicitly defined angular pivotal observable and new-size or exact-control
evidence.

## Reproduction

```bash
python scripts/exact_pivotal_russo.py \
  --geometry axis --L 2 --p 0.317 --channel cross

python -m unittest tests.test_exact_pivotal_russo -v
```
