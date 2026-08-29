# Galois groups of committed axis matching polynomials

Status: C5 finite certificates for axis L=2,3,4. Not a statement about infinite-volume `p_c`. Axis L=5 / PR #84 is excluded.

## Why this slice exists

Issue #113 closed the local complex-zero scaling route on the same polynomials. The remaining exact algebraic question in #104 is whether the finite physical roots live in a small solvable extension or in a full symmetric group.

Coefficients are the committed power-basis polynomials in `results/exact_small_matching_polynomials.md`. Pairwise gcds over `Q` are 1.

## L=2 is C4, not S4

```text
M(p) = -2 p^4 + 4 p^2 - 1
```

The polynomial is even in `p`. Writing `u=p^2` produces `g(u)=-2u^2+4u-1` with discriminant 8, not a square, so `Q(u)=Q(√2)`. The four roots `±√u_±` live in a degree-4 splitting field `Q(√2, √u_+)`. Combined with irreducibility (reduction mod 3 is a 4-cycle) the transitive subgroup of `S4` of order 4 that contains 4-cycles is `C4`.

This is the `p → -p` artifact of a tiny biquadratic, not a closed-form mechanism for `p_c`.

## L=3 is S9 and L=4 is S16

Dedekind–Frobenius, using only square-free reductions that do not drop degree:

| L | irreducible | (n−1)-cycle | transposition-bearing | group |
|---|---|---|---|---|
| 3 | p=5, type 9 | p=23, type 8+1 | p=11, type 5+2+1+1 | S9 |
| 4 | p=5, type 16 | p=19, type 15+1 | p=31, type 11+3+2 | S16 |

An (n−1)-cycle in a transitive group is primitive: the unique fixed point cannot sit in a nontrivial block. An element of type 5+2+1+1 (resp. 11+3+2) raised to the 5th (resp. 33rd) power is a transposition. A primitive permutation group containing a transposition is `S_n`.

## Boundary

These finite roots at L=3,4 are not solvable by radicals. That does **not** imply the infinite threshold is transcendental, nor that it lacks a special-function expression. Diamond polynomials and axis L=5 are not certified here.
