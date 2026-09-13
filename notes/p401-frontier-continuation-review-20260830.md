# Matching One frontier: continuation geometry, typed arms, and specialized extensions

Date: 2026-08-30. Main snapshot: `6f711e524d9b33693703109521e592ca893e9cb0`
(commit time 12:11:25 UTC; main was checked again before this contribution).
Research discussion was followed through `625d602`, `cfb3ead`, `e7e6c80`, and
the concurrently created #398/#400. This is a bounded research review, not a
claim to have read or retested every repository file.

New entrypoints: [#401](https://github.com/LightChainr/Matching-One/issues/401)
(continuation process), [#408](https://github.com/LightChainr/Matching-One/issues/408)
(typed-arm probability), and the [specialization control on #398](https://github.com/LightChainr/Matching-One/issues/398#issuecomment-5468824176).

## 1. Where the frontier actually moved

The prior `fee33287` proposal is no longer the frontier. The production age
analysis `2e99533` and saturation/site controls `742a8b0` already show a
resolved age dependence at N325/N425. The current-k0 pilot `625d602` already
computes the exact number of one-site rank-two triggers. Its `h2` means target
rank two, not a two-step time horizon.

Likewise, `1c93bed` already distinguishes theta and figure-eight direct births,
and `cfb3ead` supplies globally typed local surgeries. The remaining lower-bound
problem is the probability of the required exterior type, not the existence
of the local completion. The four-generation analysis `b887ef3` has already
rejected a correction-free common-amplitude N^(-5/6) line, while the largest
size ratio is close to that conditional scaling target.

The scalar-mark family exclusions `081a5ed`/`e7e6c80` are also stronger than
an unsuccessful fit. Concurrent #398 proposes Q-adic/finite-algebra structure;
#400 correctly separates additive-translation rank from dilation/field rank.
Neither is proposed again here.

The useful next change is from scalar summaries to **typed conditional
prediction**. Three different questions need different objects:

- occupation growth: a finite-horizon continuation kernel;
- rare topology creation: a global landing-type probability conditional on arms;
- Q-degeneration: an extension that survives specialization to the physical algebra.

Their ranks or nilpotents must not be identified just because they admit similar
matrix notation.

## 2. Exact new continuation theorem and the 1/44 witness

For a rank-one configuration S, put d=N-|S|. Let x be its number of vacant
one-site triggers. Let c2 count unordered pairs of vacancies that create rank
two together but not separately. For d>=2, uniform continuation gives exactly

```
q2(S) = P(K2<=|S|+2 | S)
      = 1-C(d-x,2)/C(d,2)+c2/C(d,2).
```

Proof: pairs containing at least one existing trigger are counted by the
binomial difference; all other successful pairs are precisely the cooperative
pairs. Monotonicity makes a hit at the first insertion persist at the second.

For independent Bernoulli continuation with vacancy-opening probability t,
`Q_S(t)=P(r(S union A_t)=2)` satisfies

```
Q_S'(0)=x,
Q_S''(0)=2[c2-C(x,2)].
```

This second response separates cooperation from redundant triggers exactly.

On the N10 Gaussian quotient 3+i, fix `k=5, ell=(0,1), x=1`:

| K1 | ordered prefixes | weighted c2 | E[c2] | one-step hit | hit within two steps |
|---|---:|---:|---:|---:|---:|
| 4 | 960 | 2400 | 5/2 | 1/5 | 13/20 |
| 5 | 2640 | 7200 | 30/11 | 1/5 | 37/55 |

The two-step difference is **1/44**. A separate enumeration of all 30,240
ordered length-five prefixes reproduces the DP counts and the successful-pair
weighted counts 6240/17760. There are 30 current subsets in this stratum.
Masks 155 and 157 in the inherited oracle's site order have `c2=2,3` and
`q2=3/5,7/10`, despite the same current `(k,rank,line,x)`.

On N13 (3+2i), `k=6, ell=(0,1), x=0`, the K1=5/6 cohorts give `q2=3/14,17/77`,
a difference **1/154**. These are exact finite-volume witnesses, not a statement
that the numerical gap survives scaling.

## 3. The hierarchy continues: a 2/315 witness after controlling c2

This refinement was discovered after #401 was opened. At N13, fix
`k=6, ell=(0,1), x=0, c2=5`. Then:

| K1 | ordered prefixes | successful-triple weighted count | hit within three steps |
|---|---:|---:|---:|
| 5 | 9360 | 199680 | 64/105 |
| 6 | 28080 | 605280 | 194/315 |

The difference is **2/315**. Independent enumeration of all 37,440 prefixes
in this declared current-state stratum verifies the recurrence. This is not
an assertion that all 13P6 prefixes were enumerated.

The next exact identity identifies the missing structure. Remove the x
one-site triggers and let t=d-x. On the remaining vacancies form the graph
of cooperative pairs. Let e=c2, w=sum_v C(deg(v),2), z=number of its triangles,
and c3=number of minimal successful triples (no successful singleton or pair).
For d>=3,

```
q3(S) = 1-C(t,3)/C(d,3)
        +[e(t-2)-w+z+c3]/C(d,3).
```

Inclusion-exclusion over edges proves the numerator: each edge occurs in t-2
triples; two overlapping edges are counted by w; a triangle is added back
once. Minimal successful triples supply the remaining disjoint class.

The N13 masks 655 and 693 both have x=0, e=5, d=7. Their `(w,z,c3)` are
`(4,0,0)` and `(6,0,3)`, so q3 is `3/5` versus `22/35`. The new information is
**incidence/overlap and higher cooperation**, not simply another count of the
same single-site pivotal observable.

## 4. A continuation state and a direct test of its predictive sufficiency

Define `b_m(S)=#{A subset vacant(S): |A|=m, r(S union A)=1}`. Then

```
P(K2>|S|+m | S)=b_m(S)/C(d,m).
```

The whole vector is a reliability signature of future rank-one survival from
that checkpoint. Its spatial refinement is the hypergraph of inclusion-minimal
successful additions. Singleton and pair hyperedges recover x and c2.
Sufficiency for the future rank path at one checkpoint does NOT imply that the
signature evolves as a closed Markov process after each update.

A stronger experiment than another age regression is to clone a current
configuration and draw two independent future orders. If Y1,Y2 indicate
absorption by horizon m, and Z is a frozen current descriptor, then

```
E[Y1 Y2 | Z]-E[Y1 | Z]E[Y2 | Z] = Var(q_m(S) | Z).
```

This tests missing predictive geometry even when age is an insensitive proxy.
Use checkpoint-level covariance and distinct-checkpoint U-statistics for the
product-of-means term; the two futures are not independent source states.

The N10 stratum above supplies an exact positive control: under its uniform
current-subset distribution, `E[q2]=2/3`, `E[q2^2]=67/150`, and unresolved
variance given `(k,line,x)` is **1/450**. Adding c2 reduces that variance to zero
at horizon two. It does not close horizon three in general.

For a candidate near-critical limit, keep horizons on the critical occupation
scale `m~N^(5/8)` when `Delta p~N^(-3/8)` is justified. Fixed two/three-site
horizons are infinitesimal controls, not direct tests of a macroscopic memory
limit. Also, the eventual rank-two hitting probability is identically one:
ordinary infinite-horizon metastable committors are trivial here.

## 5. Typed arms: the missing amplitude, not another free exponent

#408 asks for a finite occupied/vacant corridor complex in a punctured torus
that gives the required rank-zero exterior and independent deck addresses.
If its aspect ratios stay bounded, arm separation, near-critical RSW and the
existing 24-site surgery could supply a scale-uniform typed-extension bound.
This is a proposed proof strategy; simultaneous vacant barriers and exclusion
of unwanted cycles are part of the proof, not an assumed picture.

Conditional on the arm exponents, positive typed amplitudes and sufficient
near-critical tail control, the distinct path-mass predictions are

```
D_theta = N^(-5/6+o(1)),
D_figure8 = N^(-2+o(1)),
D_figure8/D_theta = N^(-7/6+o(1)).
```

The figure-eight intensity at criticality has candidate power L^(-13/4), and
its integrated mass has power L^-4. These coincide with familiar matching
powers but do not identify mechanisms. Scalar eight-arm weights have
`(h,hbar)=(21/8,21/8)`; thermal Q4 has `(37/8,5/8)` or its conjugate. Equal
dimension 21/4 is not equal spin, source, parity or observable.

The useful acquisition is the conditional frequency of the global carrier/
deck type among separated arm configurations. It separates local rarity from
global topology. Do not fit an extra power to the total mass before classifying
actual events into the channels already justified by the topology theorem.

## 6. A specialization control for #398

Every first-order family has a dual-number jet

```
J(A)=[[A0,0],[A1,A0]], N=[[0,0],[I,0]], N^2=0, [N,J(A)]=0.
```

If `A(t)^T G(t)=G(t)A(t)`, then the coefficient-of-t pairing
`B=[[G1,G0],[G0,0]]` obeys `J(A)^T B=B J(A)`. This is an exact bookkeeping
construction, not a solution of the repository's constrained endpoint/radical
problem on the same module.

Already `D(Q)=Q, G(Q)=1` gives `J(D)=[[1,0],[1,1]]`, B-self-adjoint, but

```
J(D)^2-J(D)=N !=0,
J(D)^2=(I+N)J(D).
```

It represents the parameter-thickened relation, not the Q=1 relation D^2=D.
The manufactured Jordan block is not a physical extension at fixed Q.
Distinguish base-parameter nilpotence, fixed-Q radical extensions, context
multiplication nilpotence and physical dilation nilpotence; identify them only
with explicit specialized relations and intertwiners. This refinement was
posted to #398 rather than creating another parallel Jantzen issue.

## 7. arXiv reading map and what transfers

The following are primary sources checked during this review. They are not
all new publications; old foundations remain necessary for the 2026 work.
Dates below are arXiv submission/revision dates, not HTML rendering dates.

| Source | Relevant input | Boundary for this project |
|---|---|---|
| [Smirnov--Werner, math/0109120](https://arxiv.org/abs/math/0109120) | triangular percolation exponents | no square-site exponent theorem |
| [GPS, 1008.1378](https://arxiv.org/abs/1008.1378) | pivotal measure and conformal rate scaling | not a rank-only Markov theorem |
| [GPS, 1305.5526](https://arxiv.org/html/1305.5526), Sections 6, 11 | six-arm stability; Markov property of the full near-critical configuration | projection to rank/line can lose predictive information |
| [Du--Gao--Li--Zhuang, 2205.15901](https://arxiv.org/abs/2205.15901) | sharp triangular-site arm asymptotics | no typed torus landing-frequency bound |
| [Sun--Xu--Zhuang, 2410.04767v2](https://arxiv.org/html/2410.04767v2), revised 2025-02-09 | exact annulus channels, including logarithmic expansions | one-/two-arm formulas are not six-/eight-arm torus gluing |
| [Camia--Feng, 2407.04246](https://arxiv.org/abs/2407.04246) | multipoint arm/pivotal probabilities and OPE | cooperation c2 has not been identified with a specific insertion |
| [Camia--Feng, 2508.16047v2](https://arxiv.org/html/2508.16047v2), revised 2026-06-01 | explicit energy/log-partner lattice observables and correlations | a precise control, not automatic overlap with global matching |
| [Roux--Ribault--Jacobsen, 2604.24491](https://arxiv.org/html/2604.24491), 2026-04-27 | torus one-point/sphere four-point correspondence and modular solutions | representation/modulus test only after the insertion dictionary |
| [Gao--Liu--Tse, 2311.07795](https://arxiv.org/abs/2311.07795) | backward equations and controlled transition paths | use finite horizons here; eventual absorption is certain |
| [Naguszewski--Quigley, 2607.08207](https://arxiv.org/abs/2607.08207), 2026-07-09 | rate accuracy versus pointwise committor accuracy in Ising nucleation | different model; motivates a precisely declared predictive objective |
| [Martin--Senecal--Spencer, 2601.17445](https://arxiv.org/html/2601.17445), 2026-01-24 | TL cell-module structure and Jantzen-like filtrations | #398 already owns this direction; compute the actual module's filtration |
| [Bernardi--Jelisiejew--Reig Fite, 2606.30600](https://arxiv.org/html/2606.30600), 2026-06-29 | Hankel and multiplication-tensor completion equivalence | finite-algebra/flatness hypotheses and marked generators matter |
| [Giaquinto--Mastnak, 2607.12247](https://arxiv.org/html/2607.12247), 2026-07-14 | ordinary TL center via cell filtration, radicals and deformation | not automatically periodic TL or the project's Q-lift |
| [Ikhlef--Morin-Duchesne, 2312.14837](https://arxiv.org/abs/2312.14837) | periodic-TL irreducibles/fusion and non-generic modules | a fixed-algebra control, not proof for a new rooted module |
| [Jacobsen--Ribault--Saleur, 2208.14298](https://arxiv.org/abs/2208.14298) | Potts/O(n) state spaces and symmetry labels | preserves the distinction between equal exponents and equal sectors |

The actionable common lesson is not "everything is a new Jordan mode".
Choose a generator, a source and a future/geometry question first, and require
the resulting conditional object to predict an orthogonal observation.

## 8. Reproduction, provenance and claim boundary

The unchanged source `scripts/p334_birth_age_collision_review_20260830.py`
from `fee33287` is included by its existing Git blob
`62e06795fdfa91a956aedd62b7344e84aa5efc5c`; no second topology implementation is
introduced. It was not present on the main snapshot, so the research branch
carries this explicit dependency without merging its unrelated ancestry.

Run from the repository root:

```sh
python scripts/p401_cooperative_continuation.py --output /tmp/p401-exact.json
python -m unittest discover -s tests -p 'test_p401_*.py' -v
```

Exact arithmetic: integers and fractions.Fraction. Four bounded quotients
supply 9,760 subset states, 2,822 rank-one states, 2,822 two-step gates and
2,812 three-step gates. Full survival signatures are checked on all 2,822
rank-one states. Thirteen focused tests pass, including the inherited 1/57
witness, the independent N10 prefix census, the conditioned N13 prefix census,
future-order verification, the full Bernstein expansion and corrupt-input
rejection. See `results/p401-cooperative-continuation/exact.json`.

These are finite exact statements. No asymptotic memory law, global typed-arm
lower bound, square-site universality theorem, physical Q4 overlap, production
simulation or repository-wide CI is claimed. Source hashes and local execution
metadata accompany the certificate. No old result or frozen prediction is
rewritten.
