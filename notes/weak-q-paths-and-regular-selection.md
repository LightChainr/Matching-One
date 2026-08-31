# A fixed Q-path control reverses the tangent; a regular endpoint zero stays zero

## The two named Q paths are not the same finite response

Execution's [fixed Q-family identity, f85e6d6d](https://github.com/LightChainr/Matching-One/blob/f85e6d6d8c56afcad27ff73dbf0a3578c2157829/notes/closed-source-q-lift-and-thermal-quotient.md)
selects the positive finite measure

`P proportional to a^K exp(eta B) Q^(C_B-r/2)`.

The closed source follows eta=log Q. The rank-projected ordinary site-RC
path has eta=0. Both coincide at the same iid Q=1 baseline. For the
original root/slope-normalized U, constants and the common occupation K
source vanish exactly after following the moving pooled root. Consequently

```text
J_tied = (1/2) V_Sstar,
J_site = V_(C_B-r/2) = (1/2) V_Sstar - V_B,
V_B = V_Bvac,             Bvac=2N-4K+B.
```

These are derivatives with respect to log Q at Q=1 (also d/dQ there),
not derivatives with respect to t; Q=exp(2t) fixes the factor1/2.
The [completed exact reduction, e87d5de2](../results/weak-q-path-comparison/REPORT.md)
gives `J_tied=+.06308268170708457`, `J_site=-.2698280267134874` and
`V_B=+.33291070842057197`, with strict rational signs before the positive
area factor is evaluated numerically. The fixed local edge contribution
reverses the Q-path response. This excludes equality of the two named
finite tangents, without adding a fitted coefficient.

The saved [complete N25 response packet, ec01768f](https://github.com/LightChainr/Matching-One/blob/ec01768f520e85f1acfd9d3fde9bcf855477254e/results/p337-closed-source-n25/latest.json)
already supplies the two needed rational enclosures. Their prescribed
linear reduction is [reproducible here](../analysis/weak_q_path_comparison_contract.json).
No added source, independent evidence, enumeration or root search is needed.
The ordinary path here retains Q^(-r/2); it is not unprojected site-RC.

This is a direct Q-path comparison, not a field assignment from one size.
It explains why the exact local B control in the new cross-size proposal
must be retained instead of treating the tied source as changing only a
putatively universal Q. Strong-coupling Sstar/Sdrop tails concern a different
comparison and limit order; their still-open finite-coupling window is not
silently declared solved by this weak-Q calculation.

## Regular ordinary four-leg activation is excluded throughout this endpoint family

The [existing pair-carrier projectors, d006f9c1](https://github.com/LightChainr/Matching-One/blob/d006f9c1c2f19933e5510fe2289ae418f3424be0/notes/p262-confluent-potts-projector-tomography.md)
and [endpoint selection rule, 93206494](https://github.com/LightChainr/Matching-One/blob/932064943255b8d8506f200612ab0f298c9b7092/notes/global-matching-spin4-selection-rule.md)
give an explicit rational continuation, not just zero at the isolated
percolation value. In the unordered distinct-colour pair carrier put

```text
n=Q(Q-1)/2,
P0=2J/[Q(Q-1)],
P1=(X-4J/Q)/(Q-2),
P2=I-P0-P1.
```

For the regular invariant unlabelled endpoint ell, contraction of the
diagrams gives `ell I=ell`, `ell J=n ell`, `ell X=2(Q-1)ell`. Hence

```text
ell P0=ell,
ell P1=[2(Q-1)-4n/Q]/(Q-2) ell=0,
ell P2=0.
```

The last equality is a rational diagram identity in Q. Its contracted
singularity at Q=1 is removable; every regular Q derivative is zero.
This reasoning does not invoke the false claim that infinitely many
integer values alone determine an arbitrary analytic continuation.

For each occupation configuration, the added local-edge and ambient-rank
weights are colour-invariant scalars. The finite positive occupation sum
is analytic near Q=1 on the positive-real branch. Multiplying the zero
endpoint by these weights, summing, dividing by a nonzero partition,
taking thermal derivatives, taking the angular linear combination and
following a simple pooled root all preserve zero. Thus **if the proposed
four-leg mechanism is this regular one-insertion endpoint**, its overlap
satisfies `a_F(Q) identically0`, not merely `a_F(1)=0`. Its proposed
`beta=2 a_F'(1)/a_T(1)` is zero whenever the thermal denominator is nonzero.

This excludes the regular ordinary-endpoint version of the Q-activated
sqrt(N) mechanism in the [new weak-colour proposal](https://github.com/LightChainr/Matching-One/blob/3f37daeeff8526f1b7f90e9dea56b77b2bec6391/notes/closed-source-weak-colour-spectroscopy.md).
It does not reject every possible power activation in a torus observable.
The generic conditional sqrt(N) formula is retained as a shape diagnostic,
with its regular endpoint origin removed from the active model list.

## The actual missing interface is trace/confluence to original U

The selection theorem itself distinguishes a linear endpoint from a torus
trace: `tr P2=Q(Q-3)/2` is nonzero. For example the formal scalar
`(Q-1)tr P2` vanishes at Q=1 but has derivative -1 there. This is not
a computation of lattice U; it proves why being an unmarked scalar alone
does not establish the stronger endpoint assertion for every trace term.

Therefore a surviving activated explanation must explicitly supply its
trace, marked or finite confluent lattice-to-continuum interface. For
ordinary root-normalized U that interface has not been delivered here.
Calling the missing interface a fourth adjustable numerical model would
not make it a defined competitor.

Singular projectors are not a ready-made rescue: `(Q-1)P0 -> +2J` and
`(Q-1)P2 -> -2J`. A finite collision must combine both terms, endpoints,
bare operator and measure. For the regular zero, the full derivative is
`d_Q(ell P2)=ell' P2+ell P2'=0`; retaining P2' alone changes the question.
The four-leg dimension17/4 and thermal Q4 dimension21/4 are not equal
at Q=1, so these two states alone are not the equal-dimension collision
required by the existing confluent construction. No new spectrum is
inferred from the projector pole.

## One decision-focused next output

The [hypergraph/twist realization, 977fea92](https://github.com/LightChainr/Matching-One/blob/977fea9272c780aea19cc47f8d33324c28a1293e/notes/closed-source-hypergraph-rc-twist-projection.md)
already reconstructs original q/E from five positive partitions at m=2.
Its literal finite twists do not alone define a derivative at m=1; the
real-Q occupation-sector completion is a separately explicit prescription.
Use this existing representation to settle **which trace/confluent matrix
element, if any, contributes to the original thermal-differentiated U**.
If it reduces to the regular endpoint, stop that four-leg activation model.
If a different interface survives, state its finite normalization and
coefficient prediction before treating its sqrt(N) shape as evidence.

The regular thermal/Jordan alternatives remain typed hypotheses, not an
identified operator. A later matched-size weak-Q comparison must include
the fixed B control and declared correction allowance. No strong-m point,
new source scan, generic certificate catalogue or production job is attached
to this delivery. Completed P154/P334/F4 decisions are unchanged.
