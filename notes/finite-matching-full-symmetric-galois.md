# Full symmetric Galois groups for the current exact finite matching polynomials

## Status

This note records exact finite-size algebra. It does **not** make a statement about the algebraic/transcendental nature of the infinite-volume threshold.

For the square-site matching polynomials currently available exactly, finite-field certificates show:

```text
axis L=3:    degree 9   Galois group S_9
axis L=4:    degree 16  Galois group S_16
axis L=5:    degree 25  Galois group S_25

diamond L=2: degree 8   Galois group S_8
diamond L=3: degree 18  Galois group S_18

gaussian (3,1), N=10: degree 10, Galois group S_10
```

The Gaussian `(3,1)` coefficient vector was derived independently and is accompanied by a repository regression that re-enumerates all 1024 configurations with the canonical integer-period topology engine. Treat that row as branch evidence until the regression passes in the integration review.

Axis `L=2` is irreducible of degree 4 but has the special even form `-2p^4+4p^2-1`; it is deliberately not included in the full-symmetric claim.

The machine verifier is `scripts/certify_matching_galois_groups.py`.

## 1. Irreducibility / transitivity

Finite-field Rabin certificates prove the corresponding integer polynomials irreducible over `Q`:

```text
axis L=3 mod 5
axis L=4 mod 5
axis L=5 mod 19
diamond L=2 mod 3
diamond L=3 mod 79
gaussian (3,1) mod 31
```

Therefore the Galois group acts transitively on the `N` roots.

The polynomial degree equals the number of sites in every case listed above.

## 2. A second unramified prime supplies a transposition and a large prime cycle

At a second prime the squarefree factor degrees are:

```text
axis L=3,  mod 13:   2 + 7
axis L=4,  mod 331:  1 + 2 + 13
axis L=5,  mod 863:  2 + 23

diamond L=2, mod 127: 1 + 2 + 5
diamond L=3, mod 241: 2 + 3 + 13
gaussian (3,1), mod 13: 1 + 2 + 7
```

The verifier stores the actual monic finite-field factors, multiplies them back to the source polynomial modulo the stated prime, and proves each factor irreducible with Rabin's criterion. Hence the factorization is squarefree and Dedekind/Frobenius gives a Galois-group element with exactly the corresponding cycle type.

All cycle lengths other than the unique `2` are odd. Raising the Frobenius element to the least common multiple of those odd lengths kills all of them while leaving the 2-cycle. Thus the Galois group contains a transposition.

A different power kills the other cycles while leaving the large prime cycle:

```text
7-cycle, 13-cycle, 23-cycle, 5-cycle, 13-cycle, 7-cycle
```

respectively.

## 3. The large prime cycle forces primitivity

Suppose a transitive degree-`n` group preserved a nontrivial block system with block size `d` and `k=n/d` blocks.

For all cases above, the certified prime cycle length `r` is larger than both `d` and `k` for every possible nontrivial divisor `d` of `n`:

```text
n=9:   r=7  > 3
n=16:  r=13 > 8
n=25:  r=23 > 5
n=8:   r=5  > 4
n=18:  r=13 > 9
n=10:  r=7  > 5
```

Because `r` is prime and `r>k`, an `r`-cycle cannot act nontrivially on the `k` blocks; it would have to fix each block setwise. But then all `r` moved points would have to lie inside one block, contradicting `r>d`.

So no nontrivial block system exists: the Galois group is primitive.

## 4. Primitive + transposition gives the full symmetric group

Let `tau` be the certified transposition. Consider the graph on the roots whose edges are the conjugates `g tau g^-1` as `g` ranges over the Galois group.

The connected components of this graph form a Galois-invariant block system. The graph has an edge, so components are not all singletons. Primitivity therefore forces the graph to be connected.

Transpositions along the edges of any connected graph generate the full symmetric group. All those edge transpositions already belong to the Galois group. Therefore

\[
\boxed{G=S_n}.
\]

This argument is implemented as certificate metadata; no numerical root approximation is used.

## 5. Consequences for the physical finite matching root

Because each polynomial is irreducible, its physical root has algebraic degree exactly `N`.

Because its normal closure has Galois group `S_N`, and `S_N` is non-solvable for `N>=5`, the physical root is not expressible by radicals over `Q` in all full-symmetric cases above.

Thus the finite sequence already displays two forms of algebraic complexity growth:

1. minimal-polynomial degree grows with the number of sites;
2. the Galois group is the full symmetric group across axis, diamond, and a genuinely different primitive Gaussian orientation.

This is stronger than merely observing complex roots or failing to find a low-degree factor.

## 6. Contrast with the self-matching control

The C4 self-matching N=10 triangulation has

\[
M_{10}(p)=2I_p(3,3)-1
=(2p-1)(6p^4-12p^3+4p^2+2p+1).
\]

The physical linear factor `2p-1` is guaranteed by exact complement antisymmetry on a self-matching finite quotient. This gives the factor/GCD program a positive control: an exact self-matching mechanism produces a persistent simple physical factor, while the square-site target examples above exhibit full-degree/full-symmetric algebraic complexity instead.

## 7. What this does and does not say about the infinite threshold

It supports the negative finite-cell conclusion:

> the observed finite matching roots are not converging through a sequence of bounded-degree algebraic formulas; their exact algebraic complexity is instead growing rapidly.

It does **not** prove that the limit `p_c` is transcendental, non-algebraic, or lacks another kind of exact representation. Limits of degree-growing algebraic numbers can be algebraic or transcendental.

The appropriate next exact tests are:

- integrate/reproduce the primitive Gaussian N=10 vector from the canonical topology engine;
- add one larger primitive Gaussian quotient if exhaustive cost remains reasonable;
- keep the self-matching control as the persistent-factor positive case;
- continue modular cycle-type certificates at every newly reachable exact size.
