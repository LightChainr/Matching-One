# Direct colour channel zero, original global-U transmission nonzero

**Completed finite decision.** At the fixed N25(5,0)/(4,3) pair, in the
original Sstar law at m=2, the canonical S4 `[2,2]` seam insertion has
identically zero q and E numerators but a strictly positive response of
the original, separately normalized, pooled-root U:

```text
V_beta = d_epsilon U |epsilon=0
       = +5.440121494634842e-6.
```

The deciding reduced quantity has the exact rational enclosure

```text
582084561188942894896198083834109 / 10^40
  <= V_beta/A_N <=
582084561188942894896198083834110 / 10^40,
A_N=25^(13/8)/2.
```

Both endpoints are strictly positive. The null that this specified
normalization-only trace has zero original-U transmission is excluded.
This closes a genuine finite lattice interface: its channel is no longer
just a possible categorical trace. It is not a Q1 field-activation
measurement or an assignment of a continuum four-leg state.

## 1. What was fixed, and what was actually varied

The [contract](../analysis/p337_s4_trace_transmission_contract.json) was
committed as `55fdba789a576d8d4c507372b7834f92cf506c80` before the missing
seam counts were computed. All of the following were held fixed:

- the two N25 geometries and their period bases `(a,b),(-b,a)`;
- Q=4, m=2 and the original closed-source weights `h^K 2^-g`;
- a colour seam along the first period only;
- the S4 `[2,2]` central character, with no optimized coefficient;
- the saved original m2 pooled root, with separately normalized geometries.

The central trace is the exact combination

```text
Z_[22]=(1/6)Z_id+(1/2)Z_(12)(34)-(2/3)Z_(123).
```

For an occupation A, let b2 indicate any component with a nonzero first
deck-cycle coordinate modulo2, and let n3 count components nontrivial
modulo3. Its relative insertion is

```text
beta(A)=1/6+(1/2)1_(b2=0)-(2/3)4^-n3.
```

The varied coordinate is the **coefficient epsilon of this fixed trace**:
`w_epsilon(A)=w_star(A)[1+epsilon beta(A)]`, near epsilon=0. Beta is
bounded, so these weights are positive for sufficiently small epsilon.
This is not a derivative with respect to t, Q or a newly fitted local
source. The colour number and seam were not changed after seeing results.

## 2. The exact direct-numerator zero and the surviving channel

The [topological character proof](closed-source-s4-rank-one-filter.md)
shows that beta vanishes configuration by configuration on topology
rank0 and rank2. On rank1, q=E=0. Hence

```text
sum_A q(A) w(A) beta(A)=0,
sum_A E(A) w(A) beta(A)=0,
```

as complete thermal polynomials, not merely at the measured root.
The reason is the actual seam closure: rank0 gives a trivial colour
representation; rank2 gives the point-colour representation of its unique
full-rank component. Neither contains `[2,2]`. This is stronger than
testing an invariant linear endpoint.

The partition trace survives on rank1. Two separated essential clusters,
for example, give exactly beta=1/8 at Q4, as shown in the
[physical two-cluster construction](closed-source-two-winding-cluster-trace.md).
Its full contribution is signed in general and must not be called a
sector probability.

Let `f_g=Z_[22],g/Z_g`. Then
`delta q_g=-f_g q_g`, `delta E_g=-f_g E_g`. The new counts give, at the
unchanged pooled root and with h as common thermal coordinate,

| Geometry | f | d_h f |
|---|---:|---:|
| axis | +6.459778331458177e-7 | -6.573706791139839e-7 |
| tilted | +3.476525412128839e-8 | +1.778438505610881e-8 |

The trace fractions are positive at this point; the construction itself
does not guarantee positivity for every geometry or connection pattern.

## 3. The response is carried by the geometric thermal difference

The [exact transmission formula](closed-source-rank-one-trace-transmission.md)
retains the complete pooled root and thermal denominator:

```text
V_beta/A_N=C_c d_h f_c+C_dh d_h f_d+C_d f_d,
f_c=(f_axis+f_tilted)/2,  f_d=(f_axis-f_tilted)/2.
```

The coefficients are fixed by the old q/E curves, not fitted to the new
trace. The three signed contributions to the same total are

| Prescribed term | Contribution to V_beta |
|---|---:|
| common fractional trace's thermal derivative | -3.145535454562741e-8 |
| geometric fractional trace's thermal derivative | +5.531554355406087e-6 |
| geometric fractional trace value, including root motion | -5.997750622561810e-8 |
| **Total** | **+5.440121494634842e-6** |

This is an operator decomposition, not population shares or separate
evidence votes. The root derivative is also retained:
`d_epsilon h0=+1.688717003470896e-10`. Setting each geometric q to zero
or freezing the pooled root would change the estimand.

In particular a common, thermally constant fractional trace would cancel
from U exactly. The present result survives because the actual two
geometries have different thermal trace variations. The finite mechanism
has therefore been localized to a specified normalization transmission,
without adding a descriptor to explain a residual.

## 4. Relation to the weak-Q and endpoint questions

The regular unlabelled endpoint identity `ell P_[2](Q)=0` and all its
regular Q derivatives remain zero. The present insertion is a concrete
closed torus trace at Q4, not that endpoint. Nor does an S4 isotypic
projection isolate a unique CFT state: several lattice states and scaling
fields can occur in the same representation.

For Q1 activation the [removable twist-jet interface](closed-source-removable-twist-jet-interface.md)
now specifies a separate requirement. Since
`J=R/(sqrt(Q)-1)`, one needs
`J0=2R_Q|1` and `J_logQ|1=R_QQ|1+R_Q|1/2`. Its complete finite
landing must have a nonzero thermal-quotient derivative. A nonzero value
at Q4 does not provide those Q1 jets, justify a fitted sqrt(N) term, or
reverse the completed regular-endpoint exclusion.

The next question is accordingly concrete: obtain the specified generic-Q
trace landing in that same physical family and its orientation-resolved
thermal derivative. No further Q4 point or seam search is needed to
establish the finite transmission already decided here.

**Subsequent result:** that specific Q1 continuation and transmission are
now completed in the [two-kernel landing](closed-source-two-trace-kernels-q1.md)
and [fixed Q1 score](p337-q1-closed-trace-transmission-result.md). The
primary response is strictly negative, `-0.001904836180602413`. Its trace
already has a nonzero Q1 baseline; it is not a first derivative activation
of the regular endpoint. The old Q4 result alone was not used to infer it.

## Lifecycle and scientific card

- **Freeze:** `55fdba78`; scorer `0d58d3ab` with the metadata-only change
  `37e9114b`, before this score. Source-count code `0b2ad6ff` was committed
  before production; raw receipt `0eaf5985` was imported as `3dbdd91e`.
- **Computation:** one exhaustive pass per geometry for the missing
  per-component seam constraints,1.835s and1.618s locally. One rational
  score,0.141s. Old root, D and U intervals were imported unchanged;
  baseline jets needed by the new response came from the locked histogram.
  No old four-point experiment, root search, MC, cloud job or test suite
  was run.
- **Result:** `54352b2e`; full fractions, source hashes and exact trace
  polynomials in [score.json](../results/p337-s4-trace-transmission/score/score.json).
- **Dependency:** the same finite N25 occupation population, now with one
  predetermined character insertion. It is not independent stochastic
  evidence and does not rescue the stopped P154/P334/F4 production lines.
- **Changed mechanism:** direct q/E character coupling is zero, yet the
  fixed closed-seam normalization route has nonzero original-U response.
  The universal claim that normalization-only colour traces cannot affect
  this ratio is false.
