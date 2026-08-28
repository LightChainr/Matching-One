# P44 exact central-parity gate

## Result

The smallest nondegenerate compatible Gaussian quotient, `(a,b)=(3,1)`, has
`N=10`, 30 edges, and 1024 site configurations. Exhaustive enumeration gives:

- the lifted graph agrees edge-for-edge with the cyclic checkerboard
  triangulation;
- primal and site-matching edge sets are identical;
- all 1024 masks satisfy complement antisymmetry for each recorded wrapping
  channel;
- every channel has anti-palindromic integer Bernstein coefficients;
- `M(1/2)=0` exactly for `cross`, `both`, `either`, `direction_0`, and
  `direction_1`;
- the C4-invariant channels `cross`, `both`, and `either` have zero failures
  under the exact quarter-turn automorphism.

For all five channels the black-minus-complement integer coefficients by
occupation count `k=0,...,10` are

```text
[-1, -10, -45, -100, -100, 0, 100, 100, 45, 10, 1]
```

The individual wrapping coefficients are retained in `exact.json`; unlike the
difference coefficients, they depend on the wrapping channel.

## Evidence boundary

This proves a finite-quotient implementation invariant at the critical
center. It does **not** prove `M(p)` is identically zero. On this same quotient,
every recorded channel gives the exact off-center value

```text
M(1/3) = -47/81.
```

It also does not measure the proposed `N^-1 cos(4 theta)` amplitude; that
requires the fresh N=130/170 stochastic orientation experiment in Issue #44.

## Reproduction

```bash
python3 scripts/c4_self_matching_exact.py \
  --a 3 --b 1 \
  --json results/local-20260828/P44-c4-self-matching-exact/exact.json

python3 -m unittest discover -s tests \
  -p 'test_c4_self_matching_exact.py' -v
```

The three focused tests cover lifted/cyclic graph equivalence, exhaustive
central parity, and the explicit rejection of the false `M(p) == 0` claim.
