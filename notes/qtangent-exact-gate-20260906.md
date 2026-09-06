# The two typed pieces of the Q tangent are separated at Boolean degree one, exactly

**Date:** 2026-09-06
**Ticket:** #581 (exact gate, items 1–6)
**Claim level:** C5 — finite exact identities on tiny tori. No lattice amplitude claim.
**Artifact:** `results/qtangent-scale-decomposition/latest.json`
**Script:** `scripts/qtangent_scale_decomposition.py`

#581 asks a deterministic gate before any Monte Carlo. This is that gate, run to
exhaustion on the L=2 and L=3 square-bond tori — 256 and 262,144 configurations,
every expectation an exact dyadic rational. It adds no evidence. Its job is to
establish that the objects are what they are claimed to be, and that the
decomposition has the power to tell them apart.

## The gate

On an `L×L` square-bond torus at `p = 1/2`, with `A` the open edges, `A*` the
geometric dual transport, `T = k + |A|/2` and `X = r(A) − 1`:

| check | L=2 | L=3 |
|---|---:|---:|
| `T − T* = X` | 0 failures / 256 | 0 failures / 262,144 |
| `T = V/2 + B/2` | 0 | 0 |
| `B − B* = 2X` | 0 | 0 |
| image rank = largest component rank | 0 | 0 |
| `X` values realised | −1, 0, +1 | −1, 0, +1 |

Two of those deserve a word. Only `T − T* = X` is independent; the other two
follow algebraically from it and are checked rather than presented as separate
evidence. And the fourth is the one I would otherwise have assumed: `X` needs the
rank of the *image* in `H₁(T²)`, the span of every component's winding lattice,
and the standard argument that two disjoint non-contractible cycles on a torus
must be parallel makes that equal the largest single component's rank. That is an
argument; the table is the measurement.

The gate has power. Substituting naive bit-complement for the dual transport —
the wrong dual this repository has already had to rule out once — breaks
`T − T* = X` on most configurations. A test pins that, so passing the gate is not
something the construction gives away for free.

## The split, exactly

`T = V/2 + B_even/2 + X/2` with `B_even = [B(A) + B(A*)]/2`, so for every
observable `Cov(O,T)` splits into a duality-even Betti piece and an ambient
homology piece. Exact for all six observables at both sizes:

| observable | Betti-even (L=3) | ambient X (L=3) |
|---|---:|---:|
| `wrap_either` | −0.01353 | +0.14393 |
| `wrap_cross` | +0.01353 | +0.14393 |
| `wrap_direction_0` | −0.00480 | +0.14393 |
| `open_edges` | **0 exactly** | +0.59882 |
| `components` | +0.31289 | −0.15548 |
| `cycle_rank` | +0.31289 | +0.44334 |

`open_edges` is duality-odd on a self-dual torus, so its even piece is exactly
zero at both sizes — a free check that the even projection is even. And
`components` is the reminder that the "topological fraction" is not a share of
variance: its two pieces have opposite signs and largely cancel, so the ratio
sits outside `[0,1]`.

## The result: degree one separates the two pieces, by parity

This is the sharpest thing the gate produced, and it was not what I went looking
for.

The dual transport sends a single-bond Walsh character to minus the character of
the crossing bond, `χ_i → −χ_σ(i)`. So at degree one the duality-**even** part of
any function keeps the *antisymmetric* combination `(f̂({i}) − f̂({σ(i)}))/2` and
the duality-**odd** part keeps the symmetric one. On the square torus every bond
is equivalent to its crossing partner, so the even piece's degree-one weight
vanishes identically and the odd piece's doubles.

Both ingredients measured, not argued:

```text
B_even : every degree-one coefficient is exactly 0            (L=2 and L=3)
X      : every degree-one coefficient equals one value        -27/128   (L=2)
                                                          -8721/65536   (L=3)
f^({i}) = f^({sigma(i)}) for every bond, both pieces, both sizes
```

The consequence, stated carefully: **the degree-one part of any observable's
response to the Q score is entirely ambient homology.** Exactly, at both sizes
enumerated. The parity mechanism explains it and makes it plausible for general
`L`; it is not proved here for general `L`.

The full cross-spectra bear it out, and the separation strengthens with size
(L=3, normalised by the total covariance):

```text
wrap_either x X       degree 1: +0.554   2: +0.135   3: +0.153   4: +0.074  ...  monotone decay
wrap_either x B_even  degree 1:  0       2: +1.987   3: -0.810   4: +0.083  ...  alternating, long tail to 18
```

The ambient channel is low-degree and single-signed. The Betti-even channel has
*no* degree-one weight at all, is dominated by a degree-2 term almost twice the
total, and pays for it with an alternating tail out to the maximum degree.

## The scale axis agrees

Nested spatial filtration by torus radius, `Γ_j = E[D_j D_jᵀ]` on the joint vector
`(wrap_either, 2·B_even, X)`. Positive semi-definite at every level and
telescoping exactly to the total covariance, checked in rational arithmetic:

```text
L=3   level  bonds   Gamma[wrap,B_even]   Gamma[wrap,X]
        0      4        +0.00955            +0.03932
        1     12        -0.08891            +0.10340
        2     14        +0.02330            +0.04036
        3     18        +0.00193            +0.10477
```

Same qualitative picture from the other coordinate system: the ambient coupling
is positive at every scale, the Betti-even coupling alternates. #581 warns not to
count the two axes as independent evidence — they reuse one enumeration — and
they are not counted that way here. The point is that a mechanism with a real
signature should show it in both, and this one does.

## What this does and does not license

**Does:** #581's stop rule says to abandon the programme if the decomposition has
poor power on square-bond controls. It does not. The two typed pieces are
separated at the first Boolean degree by a parity argument, which is about as
clean a separation as a decomposition can have, and the scale axis agrees.

**Does not:** every number here is a toy. L=2 and L=3 have two and four
filtration levels, and no `O_H4` appears — the square-*site* matching-odd
observable has no canonical lift to a torus of side 2 or 3, and #581 says
explicitly not to invent one. So this gate says the question is *askable*. It
does not answer it.

Also not established: that `B_even` is a local energy field (it is a duality-even
measure tangent); that a covariance with `X` demonstrates a topological defect
theory; that the degree-one vanishing holds for general `L`.

## What comes next in #581's order

The first empirical control: one exact-critical square-bond block large enough to
carry three to five predeclared Euclidean scales, with a predeclared `O_H4` that
the repository already supports. The gate says the instrument works. The
measurement is the next thing to build, and the degree-one result gives it a
sharp first question — **how much of the declared response's Q-score coupling
sits at degree one**, since whatever does is ambient by parity and needs no fit at
all.
