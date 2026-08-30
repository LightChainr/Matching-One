# Formal Gaussian Jordan conjugation certificate

The existing Issue #145 semigroup oracle proves multiplicativity for norm-2 and
norm-5 Gaussian multipliers.  This independent continuation adds the compatible
conjugation structure.

For each declared multiplier `z` in `{1+i, 2+i, 1+3i}`, Fraction arithmetic
checks

```text
chi(conj z) = conj chi(z),
T(conj z) = conj T(z),
T(z) T(conj z) = T(conj z) T(z) = T(N(z)).
```

Here conjugation on `T` acts entrywise on the exact complex coefficients and
leaves the independent formal symbols `log(2)` and `log(5)` fixed.  The
norm-10 multiplier is also checked both directly and as `(1+i)(2+i)`, including
the conjugated composite path.

Run:

```text
python3 scripts/formal_jordan_conjugate_semigroup.py
python3 -m unittest tests/test_formal_jordan_conjugate_semigroup.py -v
```

## Boundary

This is a formal semigroup/*-compatibility certificate.  It does not adjoin or
construct group inverses, identify a Hecke operator, connect the matrices to a
physical coarse-graining map, or establish that the observed sector is a
unique Jordan pair.  Issue #145 remains open.
