# Cyclic deck invariant-monomial census

Issue #244 gives the exact charge-neutrality rule for cyclic deck characters.
This certificate turns that rule into a bounded invariant-ring census.

For `C_Q`, introduce one variable `z_r` for each nontrivial charge
`r=1,...,Q-1`.  A monomial with exponent vector `e` is invariant exactly when

```text
sum_r r e_r = 0 mod Q.
```

Every exponent vector through total degree five is enumerated exactly.

## Frozen counts

Including the constant monomial, the degree-zero through degree-five counts
are

```text
C2: 1, 0, 1, 0, 1, 0
C5: 1, 0, 2, 4, 7, 12.
```

The two C5 quadratic invariants are precisely the conjugate pairs `z1*z4`
and `z2*z3`.  There are four primitive cubic generators, four primitive
quartic generators, and at degree five the four pure powers
`z1^5,z2^5,z3^5,z4^5`.  Together with the two quadratics these give 14
primitive neutral generators visible through the bound.  Composite neutral
monomials, such as products of the quadratic pairs, are deliberately not
misreported as new generators.

## Reproduction

```text
python3 scripts/cyclic_deck_invariant_monomials.py
python3 -m unittest tests/test_cyclic_deck_invariant_monomials.py -v
```

## Boundary

The census states which monomials symmetry allows.  It does not prove nonzero
lattice overlap, compute a transfer amplitude, read measured responses,
identify an operator, or recommend production.  Issue #244 remains open.
