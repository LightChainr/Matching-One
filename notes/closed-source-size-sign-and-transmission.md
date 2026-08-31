# The fixed source has a sign cycle, a winding-size law, and opposing defect channels

This turn changes the mechanism of the already specified source
`S*=C+F+Bvac`. It adds no source coefficient, finite-coupling grid point,
Monte Carlo stream or replay of the full old defect response.

## 1. The global response eventually reverses sign

The [exact N25 analytic germ](closed-source-angular-strong-coupling-law.md)
gives, with lambda=exp(-t), A=25^(13/8)/2 and DeltaCos4=1152/625,

```text
U/A=-(625/1152)lambda^11+(390625/1152)lambda^13+O(lambda^15).
```

This goes beyond the preceding U->0 bound: the eventual sign is negative,
and U_t is eventually positive. The previously computed U(log16)>0 now
forces at least one later zero crossing and a negative minimum. Along
with the earlier positive local maximum, the qualitative sequence includes

```text
weak-coupling increase -> positive peak -> decrease
-> at least one sign crossing -> negative valley -> recovery toward0.
```

This establishes existence, not uniqueness, of the extrema and crossings.
No extra coupling point was evaluated to find them. In particular the two
displayed series terms are not used to locate a zero by cancellation.

The negative tail is a specified geometric effect. Axis rank-one weight
starts at lambda^9, but reciprocal normalization cancels its thermal slope.
At lambda^11 the fixed-root slope5025/2 loses2500 to the actual root shift,
leaving25/2. The original E=1-P(rank1) projector then fixes the negative
sign. The stronger positive degree13 term explains why four positive
finite-coupling values did not reveal the eventual sign.

## 2. The next size is predicted by geometry, not by another scan

The [winding-barrier theorem](closed-source-winding-barrier.md) proves

```text
min_(rank1) g = 2 ell1(Lambda)-1,
ell1((a+ib)Z[i])=a+b,  a>=b>=0.
```

Thus under uniform period dilation k, the two equal-area quotients have
barriers10k-1 and14k-1. This is a cost linear in length, not in area.
Boundary smoothing and extra components are explicitly included in the
proof; a shortest winding cycle attains the bound.

The [square-family leading-law derivation](closed-source-square-family-leading-law.md)
also counts the next possible cost shell. For an axis LxL torus, L>=5,
paired with the same-area companion whose ell1>=L+2,

```text
U/A=-(L^2-6L+6)/DeltaCos4 * lambda^(2L+1)
     +O(lambda^(2L+3)).
```

The two added powers come from the normalized minimal-strip cancellation.
They cannot be read off from the barrier alone. For our fixed-angle family
L=5k this gives the following new size predictions:

| N | Axis / tilted rank-one barriers | Leading U/A |
|---:|---:|---|
| 25 | 9 / 13 | `-lambda^11/DeltaCos4` |
| 100 | 19 / 27 | `-46 lambda^21/DeltaCos4` |
| 225 | 29 / 41 | `-141 lambda^31/DeltaCos4` |

The first is also obtained from the saved complete N25 histogram. The last
two are combinatorial predictions, not sampled or enumerated new results.
For each fixed k the leading sign is negative. The fixed-N expansion is
not being exchanged with a large-N limit.

## 3. A local bulk model does not determine the topological finite root

The [local-colour representation](closed-source-local-colour-gas.md) gives
m² active colours and a vacant state, with nonnegative local weights and
the global projection m^(-r). That projection changes the pressure density
by at most2t/N at fixed t, yet exactly cancels the m²-fold full-colour
degeneracy. If it is dropped, the matching root changes by

```text
logit p_star-logit p_drop=2t/N+O(exp(-2t)).
```

This is an explicit size-dependent distinction between bulk-equivalent
laws and a finite topological observer. It does not assume a familiar
Potts critical line or a universality class for the constrained colour gas.

There is also a strict limit-order obstruction: literal empty/full
concentration along a growing-size sequence requires
`N exp(-2t)->0`. This follows already from the exact relative weight of
the N one-site configurations. It is a necessary condition, not a
complete double-scaling law or a condition for every kind of rank-one
suppression.

## 4. The prescribed jump-only defect model is now excluded

The overview's fixed follow-up from7132f0c2 is
[completed](checkerboard-defect-reweighting-decision.md). Keep the same
S*, the same parent pair and the physical chart pA=s+(1-s)p, pB=p.
The exact one-hole insertion has two terms:

```text
E[w DeltaO] + Cov(w,O_intact),  w=exp(t DeltaS).
```

The model omitting the second term predicts its original-U mixed
contribution to vanish. The frozen new cross-moment calculation gives

| Fixed operator contribution to U_st | Value |
|---|---:|
| Baseline reweighting | -4.550327123236791 |
| Weighted rank jump, previous total minus new primary | +15.306045530800864 |
| Earlier complete total, imported without rescoring | +10.755718407564073 |

The reweighting rational enclosure strictly excludes zero, so the
jump-only model fails. The two terms oppose each other. These are signed
operator contributions, not population shares: reweighting includes both
rank-preserving and rank-changing configurations. Ignoring covariance
centering would even reverse its sign. The original root, slope, dose
derivative and per-geometry normalization are included.

## Scientific card and continuation

- **Mechanism changes:** a negative strong-coupling tail and unavoidable
  sign crossing; an explicit size law from winding shells; rejection of
  the prescribed rank-jump-only defect model.
- **Lifecycle:** exact-series code09042093/result9b88a49b uses the previously
  locked N25 histograms; barrierac2d29cd and size-law79799dbc are pure theory.
  The local-colour/root-shift argument starts at935abc89. The separate
  reweighting contracte6a900d9 precedes codedb348346 and resultd7f385b6.
- **Dependency:** series and old four-point values share the same N25
  exhaustive populations. The two new cross moments are paired within
  the already defined N50 one-hole populations. Neither is independent
  random evidence; the larger-size formula is theory, not data.
- **Cost:** exact32-degree series0.118seconds; new paired counting and
  fixed reweighting score3.70seconds. Local CPU only, no cloud or test suite.
- **Boundary:** finite lattice and fixed named source/observer; no continuous
  H4-field identity, universal critical exponent, unique zero/valley or
  thermodynamic transition. Prior F4, lag-one and P334 stop decisions stand.
- **Next substantive target:** combine the fixed defect operator with the
  winding-shell size law, rather than add descriptors or locate another
  point on the completed N25 coupling grid. A bulk-identical model which
  omits its topological projection is already a distinct finite-observer
  mechanism, not a harmless source gauge.
