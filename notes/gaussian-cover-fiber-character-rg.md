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

## Finite-volume character selection rule

Issue #244 supplies the exact gate that must precede such a response
experiment.  Let a deck translation \(T_t\) act on a regular cover with deck
group \(K\).  At uniform Bernoulli \(p\), the product measure is invariant.  A
nontrivial character score obeys

\[
 S_\chi(T_t\omega)=\chi(-t)S_\chi(\omega).
\]

If \(O\) is an unmarked deck-invariant observable, changing variables in its
linear response gives

\[
 R_\chi(O)=\frac{E_p[O S_\chi]}{p(1-p)}
           =\chi(-t)R_\chi(O).
\]

For some \(t\), a nontrivial \(\chi\) has \(\chi(-t)\ne1\), hence
\(R_\chi(O)=0\) exactly at every finite size.  Global wrapping, homology and
site-summed pivotal rows therefore cannot see an unmarked linear detail mode.
This supersedes any proposal to estimate that symmetry-forced column by a
larger Monte Carlo run.

Two channels remain legal.  First, an opposite-character marked row
\(O_{\bar\chi}\) makes \(E[O_{\bar\chi}S_\chi]\) invariant and potentially
nonzero.  Second, \(\chi\otimes\bar\chi\) permits an invariant quadratic
susceptibility.  For a real perturbation
\(p_i(\epsilon)=p+\epsilon h_i\), the required Bernoulli second score is

\[
 H_h=\left[\sum_i\frac{h_i(X_i-p)}{p(1-p)}\right]^2
 -\sum_i h_i^2\left[\frac{X_i}{p^2}
                    +\frac{1-X_i}{(1-p)^2}\right].
\]

The diagonal likelihood correction is essential: raw \(S_h^2\) is not the
second derivative.  This quadratic object is a composite susceptibility, not
a new linear RG tangent.

The norm-two executable oracle exhausts all \(2^{10}\) configurations of the
\((2+i)\to(2+i)(1+i)=1+3i\) cover at \(p=2/5\).  It verifies:

- matching-odd cross wrap has exactly zero linear detail response;
- a balanced but non-equivariant sign registry produces a nonzero response on
  a deck-even anchored control, so the test detects labeling errors;
- the opposite-character marked pivotal row has response
  \(-10944/390625\), equal to direct symbolic differentiation;
- the invariant Hessian response is \(109056/78125\), also equal to direct
  symbolic differentiation, while the raw score-product term alone differs.

These are exact finite-volume statements and require no CFT assumption.

Reproduce with:

```bash
python3 scripts/gaussian_cover_character_modes.py \
  --output results/exact-cover-character-oracles/gaussian_cover_characters.json
python3 scripts/gaussian_cover_selection_rule.py \
  --output results/exact-cover-character-oracles/norm2_selection_rule.json
python3 -m unittest discover -s tests -p 'test_gaussian_cover_character_modes.py'
```
