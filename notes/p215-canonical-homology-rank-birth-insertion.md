# Issue #215: canonical homology-rank-birth insertion

Status: exact finite insertion algebra plus exhaustive tiny/rooted oracles. This
is the minimal insertion-level continuation of the filtration result merged at
`3bb4054`; it does not duplicate that path oracle.

## 1. From the persistence clock to a local insertion

The digital-Alexander theorem and the essential-birth filtration give

\[
 M(p)=\mathbb E_p[r]-1,
 \qquad
 r=\operatorname{rank}\operatorname{im}
 [H_1(\text{black})\to H_1(T^2)]\in\{0,1,2\}.
\]

For a site `v`, freeze all other sites and write `r0=r(v=0)` and
`r1=r(v=1)`. Define the two elementary **gate crossings**

\[
 I_{01}=\mathbf 1\{r_0=0,\ r_1\ge1\},\qquad
 I_{12}=\mathbf 1\{r_0\le1,\ r_1=2\}.
\]

Because `r` is monotone and takes only the values zero, one and two,

\[
 \boxed{\Delta_vr=r_1-r_0=I_{01}+I_{12}}
\]

configuration by configuration. This includes a direct `0->2` insertion: it
crosses both gates at the same site and contributes one unit to each.

Define

\[
 f_{01}(p)=\sum_v\mathbb E_p[I_{01}],\qquad
 f_{12}(p)=\sum_v\mathbb E_p[I_{12}].
\]

The product-measure differentiation formula then gives the exact decomposition

\[
 \boxed{M'(p)=f_{01}(p)+f_{12}(p)}.
\]

Equivalently, `r=1[r>=1]+1[r=2]`: the two terms are precisely the local Russo
influences of the first and second essential-birth CDFs from #269/#276.

### Why “gate” rather than “strict transition” matters

The exhaustive oracle finds direct `0->2` site insertions on both the degenerate
axis `L=2` control and the honest axis `L=4` control. On axis `L=4`, the
translation-weighted rooted table contains

```text
0->1 : 80128
1->2 : 45312
0->2 :  4624
```

among `524288` weighted root environments. A definition using only strict
`0->1` and strict `1->2` transitions would miss twice the `0->2` mass and fail
Russo exactly. The two-gate definition is therefore not cosmetic.

## 2. Canonical line and Smith-index typing

For a nonsimultaneous `0->1` birth, record

```text
ell  = primitive saturated line of the rank-one image after insertion,
iota = index of the actual integral image subgroup inside Z ell.
```

For a nonsimultaneous `1->2` birth, there is no canonical second line in
`Z^2`: a lift of the new quotient generator is only defined modulo the existing
line. The canonical mark is instead

```text
ell  = rank-one plateau line immediately before insertion,
iota = its endpoint saturation index.
```

This is exactly the line born at `K1` and destroyed at `K2`. The filtration
theorem from `3bb4054` says it is constant across every nonempty plateau, so
both endpoint insertions carry the same projective label.

A simultaneous `0->2` insertion has no intermediate rank-one state and hence
no canonical `ell`; the schema records `ell=null` instead of choosing a
basis-dependent artificial line.

All nonsimultaneous births in the declared controls have `iota=1`. Thus the
oracle closes the data structure and the primitive line, but finds no
nontrivial Smith-index memory at these sizes.

## 3. Exact H4 marks on the same insertion

Each non-null primitive line is mapped back to a physical lifted period vector
`(x,y)` and receives the exact rational spin-four mark

\[
 \cos4\theta=\frac{x^4-6x^2y^2+y^4}{(x^2+y^2)^2},\qquad
 \sin4\theta=\frac{4xy(x^2-y^2)}{(x^2+y^2)^2}.
\]

The axis `L=4` fixed-root control additionally attaches the repository's
existing radius-one landing-sector mark

```text
axis, diagonal, both, landed, h4=axis-diagonal.
```

This produces a typed insertion containing, without changing the Russo mass,

```text
birth gate: 0->1 or 1->2,
topological direction: ell,
integral memory: iota,
homology-direction H4,
local landing H4.
```

The local landing mark is not equivalent to rank birth. On axis `L=4`, only
`36608` of `84752` first-gate environments and `20448` of `49936` second-gate
environments are landed at radius one. The old marked four-arm observable is
therefore a genuine geometric sub-selection/coordinate, not the canonical
rank-birth mass itself.

## 4. Exact polynomial certificates

The oracle expands every configuration weight in the ordinary power basis.
For each geometry it independently constructs

```text
M(p) from every full configuration,
M'(p) by exact polynomial differentiation,
f01(p), f12(p) from every fixed-root environment.
```

All coefficients are `Fraction` objects and the arrays agree identically, not
only at sampled values. At `p=1/2`:

| geometry | `M'` | `f01` | `f12` |
|---|---:|---:|---:|
| axis `L=2` degenerate | `3` | `3/2` | `3/2` |
| Gaussian `(2,1)` | `25/8` | `15/8` | `5/4` |
| axis `L=4` | `4209/1024` | `5297/2048` | `3121/2048` |

For axis `L=4`, translation symmetry reduces the calculation to `2^15`
environments of the registered origin root and multiplies by `N=16`; agreement
with the independently enumerated full `M(p)` derivative is an exact check of
that reduction.

## 5. Increment relative to `3bb4054`

`3bb4054` establishes permutation/path semantics:

```text
K1, K2, the complete R_k path, plateau ell, endpoint reflection.
```

This follow-up adds exactly the missing local layer:

```text
two-state site insertion on a frozen environment,
pointwise 0->1/1->2 gate decomposition,
direct 0->2 handling,
Russo influence polynomials,
endpoint ell/iota typing,
homology and landing-sector H4 marks.
```

It makes the unmarked sum an algebraically exact local-to-global bridge, while
retaining line and geometry labels for later charged/singlet comparisons.

## Boundary

- The gate identity and Russo decomposition are general for any monotone
  rank-valued observable in `{0,1,2}`.
- Exhaustive geometric results cover only the declared controls.
- `iota=1` here does not exclude nontrivial integral-index evolution on larger
  or covered tori.
- No canonical line is claimed for a simultaneous `0->2` jump.
- H4-resolved insertions are typed sufficient statistics, not yet a continuum
  operator identification or an independent evidence row.

Reproduce with:

```bash
python scripts/homology_rank_birth_insertion.py \
  --output results/homology-rank-birth-insertion/latest.json
python scripts/homology_rank_birth_insertion.py --format markdown \
  --output results/homology-rank-birth-insertion/latest.md
python -m unittest discover -s tests -p 'test_homology_rank_birth_insertion.py'
```

