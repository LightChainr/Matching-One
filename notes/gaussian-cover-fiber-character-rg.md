# Exact Gaussian-cover fiber character oracle

Issue #226 proposes using the deck group of a Gaussian cover as a canonical
microscopic coarse/detail basis.  The exact-first layer is finite integer
algebra: for `m=a+ib`, multiplication on the Gaussian lattice is

\[
 M_m=\begin{pmatrix}a&-b\\b&a\end{pmatrix},\qquad
 K_m=\mathbb Z^2/M_m\mathbb Z^2.
\]

For a nonsingular `2x2` matrix the Smith invariants are

\[
 d_1=\gcd(M_{ij}),\qquad d_2=|\det M|/d_1.
\]

The oracle derives rather than assumes:

| multiplier | degree | Smith invariants | additive deck group |
|---|---:|---:|---|
| `1+i` | 2 | `(1,2)` | `Z/2` |
| `2+i` | 5 | `(1,5)` | `Z/5` |
| `2-i` | 5 | `(1,5)` | `Z/5` |
| `2i` | 4 | `(2,2)` | `Z/2 x Z/2` |
| `3+i` | 10 | `(1,10)` | `Z/10` |

`3-i` is included only to close the conjugation action of the required
`3+i` quotient.

Each element is represented by a deterministic shortest integer vector.  A
dual representative `q` defines the exact character

\[
 \chi_q(x)=\exp(2\pi i\,q^T M_m^{-1}x).
\]

The machine artifact stores the rational phase modulo one and its integer
exponent relative to the Smith group exponent.  Thus no floating root-of-unity
comparison is needed.  Exact table checks verify the group homomorphism law
and character orthogonality.

## D4 and conjugation

Multiplication by the four Gaussian units acts within each quotient.
Reflection/conjugation maps `K_m` to `K_conj(m)`.  Consequently

```text
2+i <-> 2-i
3+i <-> 3-i
```

while the ideals of `1+i` and `2i` are self-conjugate up to a unit.  The JSON
stores both the permutation of element representatives and the contragredient
pushforward permutation of characters for all eight D4 operations.  At `2i`,
conjugation is the identity on the four classes, whereas a 90-degree rotation
swaps the two coordinate generators and fixes their sum.

## The exact `(1+i)^2=2i` composition

Reduction gives

\[
 0\longrightarrow \mathbb Z/2_{\rm detail}
 \longrightarrow K_{2i}\longrightarrow K_{1+i}\longrightarrow0.
\]

In the declared Gaussian basis every class has the exact split

\[
 z= r+(1+i)s\pmod {2i},\qquad r,s\in\{0,1\}.
\]

The four characters are therefore

```text
1, (-1)^s, (-1)^r, (-1)^(r+s),
```

and their exponent table is the `4x4` Hadamard table, exactly the tensor
product of the two degree-2 tables.  This identifies the coarse pullback, the
new detail, and their product without a fitted basis.

The split is basis-specific, not an assertion that every abstract extension
has a canonical splitting.  It is exactly the basis used by Gaussian cover
composition in this repository.

## Accurate norm-4 conclusion

The norm-4 group is not merely *chosen* to be `Z2 x Z2`; its Smith form forces
that structure and exponent two.  Moreover every Gaussian integer of norm 4
is an associate of `2` or `2i`, so every scalar Gaussian degree-4 cover has
the same `(2,2)` Smith structure.  A cyclic `Z4` degree-4 comparator is only
hypothetical outside this scalar Gaussian-multiplier family.

This exact oracle does not show that any nontrivial fiber character survives
as a continuum RG direction.  That question belongs to a later response
experiment, not to this commit.

Reproduce with:

```bash
python3 scripts/gaussian_cover_character_modes.py \
  --output results/exact-cover-character-oracles/gaussian_cover_characters.json
python3 -m unittest discover -s tests -p 'test_gaussian_cover_character_modes.py'
```
