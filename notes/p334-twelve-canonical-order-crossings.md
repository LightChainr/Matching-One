# Complete clock reversals survive into the canonical thermal response

The 12 fixed physical prefix clocks at
[`bd95f2a0`](https://github.com/LightChainr/Matching-One/commit/bd95f2a048d5780568b689bd42e0a684daf74315)
have a new directly usable consequence: **all 11 rank-survival crossing
pairs retain their crossings after the binomial canonical readout**. There
are 20 simple thermal crossings and 21 intervals of fixed canonical order.
The other 55 pairs remain ordered. This is a finite result for these 12
old-source prefixes, not an estimated ensemble crossing rate.

The most informative named pair, 83/1006, has a unique crossing close to
the project's p reference:

\[
\boxed{p_*=0.594353897611717.}
\]

Here and below a short label means `counter-43042500000`. No network was
reconstructed or solved, and no new prefix, continuation or p-grid search
over candidate populations was sampled.

## The exact rank-to-canonical map

Each prefix has N=425, k0=252 and d=173 remaining sites. Let f_i(k) be its
already-solved number of safe k-subsets and
`S_i(k)=f_i(k)/choose(173,k)`. The complete waiting law has
`Pr(T=t)=S_i(t-1)-S_i(t)` and K2=252+T. Consequently

\[
F_{2,i}(p)=E_T\Pr[\operatorname{Bin}(425,p)\ge252+T]
=\sum_{k=0}^{173}[1-S_i(k)]B_{425,252+k}(p),
\]

where earlier occupancy ranks have zero completion coefficient. Thus

\[
\boxed{\Delta_{ij}(p)=F_{2,i}(p)-F_{2,j}(p)
=\sum_{k=0}^{173}[S_j(k)-S_i(k)]B_{425,252+k}(p).}
\tag{1}
\]

The Bernstein coefficients are exactly the *negative* survival difference.
Canonicalization averages the whole rank sequence; it is not the substitution
`p=(252+k_cross)/425`.

## Why the single rank reversal forces exactly one thermal zero

Every Bernstein basis term is positive for 0<p<1. A coefficient sequence
with no sign changes therefore has no interior zero. Under the logit
change `x=p/(1-p)`, (1), after a positive endpoint factor is removed, is a
polynomial in positive x with the same coefficient signs. The Bernstein/
Descartes sign-variation bound limits its interior roots, counting
multiplicity, to the number of sign changes.

If there is one coefficient sign change, the leading and trailing signs
are opposite. Continuity forces at least one zero, and the variation bound
allows at most one. **A single rank-order reversal must survive as exactly
one simple thermal crossing.** It cannot disappear under this binomial
average. This says nothing about the crossing's amplitude or experimental
visibility.

For two sign changes, either zero or two simple crossings can survive
(a double contact is also an algebraic possibility). Survival is not a
general theorem in that case. We therefore subdivide the exact reduced
Bernstein polynomial until each dyadic interval has variation zero or one.
The former excludes roots and the latter certifies exactly one. Both roots
survive for each of the nine two-change pairs in this fixed set.

Specifically, if a and b are the first/last nonzero coefficient indices,
factor out `p^(252+a)*(1-p)^(173-b)`. The remaining degree b-a Bernstein
coefficients are

\[
\gamma_r=
\frac{\binom{425}{252+a+r}}{\binom{b-a}{r}}
\frac{f_j(a+r)-f_i(a+r)}{\binom{173}{a+r}}.
\]

Their denominators are cleared exactly and a positive gcd removed. The
machine certificate records the primitive coefficient hash and complete
dyadic zero/one-variation partition. Root counts are exact; their decimal
locations and slopes are numerical evaluations inside certified intervals.

## 83 versus 1006: early completion advantage, later collective disadvantage

Define `Delta=F2_83-F2_1006`. The source has
`S83-S1006<0` at steps1..10 and `>0` at11..153. Therefore

- **0<p<p_star:** Delta>0; prefix83 has the larger canonical probability of
  having completed its second birth.
- **p_star<p<1:** Delta<0; prefix1006 has the larger completion probability.

The immediate gate and the later collective clock favor different prefixes.
Their canonical ordering accordingly depends on p, despite fixed N, k0,
age and ambient line.

| Readout | Value |
|---|---:|
| p_ref | 0.59274605079 |
| F2_83(p_ref) | 0.1332277813666 |
| F2_1006(p_ref) | 0.1329298729224 |
| **Delta(p_ref)** | **+0.000297908444124** |
| F2_83'(p_ref) | 6.83721846961 |
| F2_1006'(p_ref) | 7.00731592047 |
| **Delta'(p_ref)** | **-0.170097450862** |
| F2 at crossing | 0.1444973884997 |
| Delta'(p_star) | -0.200799013623 |
| integral Delta | -0.00275241089307 |

The positive value at p_ref coexists with a negative integrated difference:
the mean wait of83 is longer by1.17252704 steps, and
`integral Delta=-(E[T83]-E[T1006])/426`. This is a conditional thermal-order
reversal, not a crossing of a new global critical-point estimator.

## All retained crossings in the fixed twelve

| Pair | Rank sign changes | Canonical crossing p values |
|---|---:|---|
| 48/106 | 2 | 0.5459069193, 0.9264269389 |
| 48/622 | 1 | 0.9960045012 |
| 83/835 | 2 | 0.7164491607, 0.9916221251 |
| 83/881 | 2 | 0.7174239222, 0.8091408713 |
| 83/1006 | 1 | 0.5943538976 |
| 83/1010 | 2 | 0.6181331818, 0.9624629094 |
| 83/1410 | 2 | 0.6834000395, 0.8336877384 |
| 835/881 | 2 | 0.7199435795, 0.7412243720 |
| 835/1010 | 2 | 0.5603724358, 0.9540825139 |
| 835/1410 | 2 | 0.6612518943, 0.7936163634 |
| 1006/1010 | 2 | 0.6365169859, 0.9219982596 |

None of these 11 crossing pairs disappears, and the 55 other pairs gain
no crossings. This does **not** make all 20 crossings useful production
targets. For 48/622 at p approximately0.9960, both completion probabilities
are one to ordinary displayed precision: their uncompleted tails are only
about `6.71e-46`, and the difference slope is `1.63e-43`. We evaluate those
tails directly from the exact survival coefficients, never by subtracting
a rounded completion probability from one.

The output also records the full order on all 21 open p intervals. Each
simple crossing swaps adjacent entries in that order, providing a compact
description of the thermal ordering of the same twelve fixed states.

## Artifacts and scientific card

Run `python3 scripts/p334_twelve_canonical_crossings.py`. Output
`results/p334-twelve-canonical-crossings/score.json` includes all66 pair
certificates, root positions/slopes, rank sign runs, reference-p differences,
integrated differences, stable rare tails, and the21 order cells. It reads
only the committed full-clock coefficients. Focused algebraic examples
include a one-change polynomial and two-change polynomials whose roots
either survive or disappear.

- **Mechanism changed:** the complete physical continuation law induces
  genuine, sometimes near-reference, canonical second-birth order reversals.
  A prefix's immediate gate or mean wait alone cannot order its whole
  thermal response.
- **Not established:** no population crossing frequency, global MC gain,
  critical-field identity, or independent replication. The83/1006 example
  was named after the earlier fixed-set clock analysis.
- **Observer / source:** canonical K2 response, 12 frozen old-source real
  N425 second-orientation prefixes, k0=252/age10/ell=(12,-19), bd95f2a0.
- **Dependency:** the same twelve exact conditional laws; all p readouts
  reuse them. No covariance rescore of the earlier two-prefix examples.
- **Next discriminant:** conditional noise weight and cost across the
  complete available prefix stratum, rather than treating one crossing
  example or eleven pair crossings as a population rate.
