# Weak colour response of the same closed source: logarithm or activated power?

The original root/slope normalization does **not** add a thermal-exponent
velocity to the proposed area-scale contrast: the thermal factors cancel,
under the precise single-field conditions below. Its regular thermal
target is

```text
[V_(cN)/U_(cN)-V_N/U_N]/log(c) -> 9 sqrt(3)/(16 pi)
                              = 0.31012250367475802785... .
```

There is a sharper live alternative than reopening the ordinary four-leg
endpoint overlap. If that overlap is zero at Q=1 but its Q derivative is
nonzero, the source response acquires a **sqrt(N)** term. Three matched
sizes distinguish its scaling from a regular thermal logarithm and from
a fixed-rank logarithmic/Jordan mixture without fitting a new exponent.

This derivation starts from `85fd4923`, imports the existing Q velocities,
and performs no Monte Carlo, finite-coupling score or velocity-oracle run.
The new numerical packet only evaluates those constants and their new
area-scale contrasts using the existing research Python environment.

## 1. The exact positive Q family, including its two extra factors

Write B for occupied NN edges, C_B for occupied NN components, r for
ambient rank, and y=p/(1-p). The same source satisfies

```text
S_star=2C_B+2B-5K-r+2N+1,
Q_c=exp(2t),  a=y Q_c^(-5/2),
mu_(p,t)(omega) proportional to a^K Q_c^(C_B+B-r/2).
```

This occupation measure is positive for every real Q_c>0, even when there
is no literal integer number of colours. For integer Q_c its local colour
representation assigns a vacant state or one of Q_c active colours to a
site; adjacent active colours must agree and each equal-colour active edge
has weight Q_c. Colour summation gives `a^K Q_c^(C_B+B)`, followed by
the explicit winding projection `Q_c^(-r/2)=m^(-r)`.

Thus the family differs from the standard site random-cluster weight
`a^K Q_c^C_B` by both the **tied local edge factor Q_c^B** and the
**winding factor Q_c^(-r/2)**. Positive weights, FKG for t>=0 and equal
pressure after removing the bounded winding factor do not prove Potts
universality of this particular tied family. Whether its critical line
has the Potts value of Q equal to this microscopic Q_c is a hypothesis,
not a normalization convention that can be imposed after observing a slope.

For a source X let V_X be its first response of this same U, including the
pooled root and denominator, in the same baseline law. A K tilt is a common
homogeneous logit shift, so V_K=0 exactly. Consequently

```text
V_Sstar/2 = V_CB+V_B-(1/2)V_r,
Q_c partial_(Q_c) U = V_CB+V_B-(1/2)V_r.
```

The first equality uses t as the source coordinate; the second uses Q_c.
Both retain the local edge and winding terms. In particular V_CB alone
is not the derivative of the specified closed family.

## 2. Why the thermal derivative and moving root do not change the power

Let `A_N=N^(13/8)/2`, `Q_N=mean(q)`, `Y_N=P4(E)`, and
`U_N=A_N Y_N,p/Q_N,p` at its own pooled Q_N=0 root. Compare geometrically
similar pairs at N and cN: same modulus, orientation projector and typed
topological observer/source. N is area and L is proportional to sqrt(N).

Assume a critical continuation with a relevant thermal coordinate
`z=a(t)[p-p_c(t)]L^(y_t(t))`, a limiting simple pooled root z_*(t), and
a single nonzero angular correction from a field of dimension x(Q_c):

```text
Q_N = mathcal Q(z,t)+smaller terms,
Y_N = u(p,t) L^(2-x(Q_c)) mathcal Y(z,t)+smaller terms,
mathcal Q(z_*,t)=0,  mathcal Q_z(z_*,t)!=0,
u(p_c(t),t) mathcal Y_z(z_*,t)!=0.
```

The source and finite-size remainders must obey this expansion uniformly
to first order in t near0. Pointwise dominance at t=0 alone does not
ensure dominance of its Q derivative; Section4 gives a counterexample.

The leading thermal derivatives are

```text
Q_N,p ~ a(t)L^(y_t) mathcal Q_z,
Y_N,p ~ a(t)u(t)L^(2-x+y_t) mathcal Y_z.
```

The derivative of the smooth microscopic coupling u contributes a term
smaller by L^(-y_t), provided the displayed leading slope overlap is
nonzero. The moving root changes the limiting scaling functions and
amplitude through z_*(t), not this power. Hence

```text
U_N(t)=mathcal A(t) N^kappa(t)[1+epsilon_N(t)],
kappa(t)=13/8-(x(Q_c(t))-2)/2,
epsilon_N, partial_t epsilon_N -> 0.
```

In logarithmic derivatives the `y_t'(t) log L` contribution to Y_p
cancels the identical contribution to Q_p. There is no extra thermal
dimension or velocity to append. Equivalently U=A_N partial_Q Y on the
thermal quotient, so a common thermal reparameterization cancels exactly.

This argument fails if the leading angular thermal slope is zero, the
asymptotic Q slope vanishes, a different sector enters, or the correction
is nonlinear in the supposed leading field. For example, with
mathcal Y_z(z_*)=0, the u_p term or root-shift corrections can dominate
with a different power. Strict finite-volume Q_p>0 does not alone prove
the nonzero limiting scaling slope required here.

## 3. Parameter-free regular targets and their observer typing

For `R_N=V_Sstar,N/U_N=partial_t log|U_N|`, the preceding expansion gives

```text
R_N(t)=partial_t log|mathcal A(t)|
       -Q_c x'_Q(Q_c) log N + o(1).
```

At t=0, `Q_c=1, dQ_c/dt=2`. Thus

```text
S_N(c)=[R_(cN)-R_N]/log c -> -x'_Q(1).
```

The area factor1/2 in the exponent and colour derivative2 cancel.
Using Q rather than t as the derivative coordinate would halve the stated
target; using a *length* ratio instead of area ratio would double it.
The factor N^(13/8) is kept fixed and must not be retuned with Q.

The [existing velocity dictionary](q-velocity-spin4-spectroscopy.md) supplies
the following conditional values; it was not recomputed here:

| Field family | x(1) | U_N power at t=0 | S_N(c) target |
|---|---:|---:|---:|
| thermal Q4 epsilon, regular nonzero slope overlap | 21/4 | 0 | `9 sqrt(3)/(16 pi)=0.310122503674758...` |
| four-leg V_(2,+/-2), **typed nonzero-overlap control only** | 17/4 | 1/2 | `5 sqrt(3)/(16 pi)=0.172290279819310...` |

Their target gap is `sqrt(3)/(4 pi)=0.137832223855448...`; it has no
amplitude parameter. At an area dilation c, the corresponding differences
R_(cN)-R_N are the displayed constants times log c. Their baseline U ratios
would tend respectively to1 and sqrt(c).

The second row is **not reopened as an ordinary endpoint competitor**.
The later [P275 double-projector selection](https://github.com/LightChainr/Matching-One/blob/ddf41aa5b02bd223329bee5e5d08beac10a04f1e/notes/p275-double-projector-selection-staircase.md),
building on `9320649`, already annihilates the regular unlabelled ordinary
four-leg one-insertion overlap. The old preflight listed its generic-loop
spectrum, not a surviving ordinary lattice coupling. Its constant velocity
is useful only for a genuinely different charged/twisted or otherwise typed
observer with an established nonzero endpoint overlap.

## 4. The live Q-activated four-leg adversary predicts a power, not its velocity

Keep the ordinary endpoint selection, but consider the named Q-derivative
loophole

```text
a_F(1)=0,  a_F'(1)!=0,
U_N(Q)=a_T(Q)N^kappa_T(Q)+a_F(Q)N^kappa_F(Q)+... .
```

At Q=1, kappa_T=0 and kappa_F=1/2. The baseline can therefore remain
thermal and of order1 while its t tangent is

```text
R_N=alpha+gamma_T log N+beta sqrt(N)+...,
gamma_T=9 sqrt(3)/(16 pi),
beta=2 a_F'(1)/a_T(1).
```

There is no first-order four-leg velocity term multiplying this activated
amplitude: differentiating its exponent also multiplies a_F(1)=0. The
velocity `5 sqrt(3)/(16 pi)` is therefore **not** the prediction for this
ordinary Q-activation mechanism. Its first response measures the overlap
derivative and the existing dimension gap1.

This is a conditional Q-residue hypothesis, not a claim that every positive
colour family activates the non-singlet. If the regular overlap theorem
extends as `a_F(Q) identically0` in the relevant analytic family, then
the activation is zero as well. A marked/singular Q continuation must be
identified rather than inferred solely from the availability of a colour
fugacity. A scalar derivative of a regular analytic zero remains zero.

Nevertheless the hypothesis has an amplitude-free three-size shape. Put

```text
T_N=R_N-gamma_T log N,
Delta_N=T_(cN)-T_N.
```

For a resolved nonzero activated contribution,

```text
Delta_(cN)/Delta_N -> sqrt(c),
R_(c^2 N)-2R_(cN)+R_N
                  ~ beta (sqrt(c)-1)^2 sqrt(N).
```

For c=4 the increment ratio is2; for c=2 it is sqrt(2). No beta or new
exponent is fitted. The ratio is undefined in the pure thermal leading
model, where both increments vanish; it must not be formed from unresolved
near-zero denominators. An equivalent division-free leading constraint is
`Delta_(cN)-sqrt(c)Delta_N=o(sqrt(N))`,
together with a separately resolved nonzero Delta. These are all contrasts
of the same R data, not independent evidence.

## 5. Jordan contamination can curve the slope, or imitate a shifted constant

For the minimal fixed-rank logarithmic ansatz, with ell=log N, write

```text
U_N=N^kappa(t)[a(t)+b(t)ell],
R_N=kappa_t ell+(a_t+b_t ell)/(a+b ell).
```

Here b is the coefficient of log N, not log L. Direct subtraction gives

```text
S_N(c)=kappa_t+
       (a b_t-b a_t)/[(a+b log N)(a+b log(cN))].
```

For nonzero b and regular coefficients this generally produces finite-size
curvature; its extra slope is of order1/(log N)^2 at large N. However the
especially important case `b(0)=0, b_t(0)!=0` gives

```text
R_N=a_t/a+[kappa_t+b_t/a]log N,
S_N(c)=kappa_t+b_t/a.
```

Thus a source-activated Jordan partner can give a **constant but shifted**
slope. Constancy alone does not establish the ordinary thermal velocity.
If a Jordan collision at Q=1 resolves into distinct generic-Q dimensions
with singular normalization, its tangent can also contain log(N)^2 terms;
the regular a,b ansatz does not exhaust that case. Fixed-degree logarithmic
responses with resolved nonzero increments have equal-log-scale increment
ratios tending to1, rather than
the sqrt(c) of a nonzero power activation, away from their zeros.

An ordinary mixture of fields with nonzero endpoint amplitudes also gives
scale-dependent weights and possible poles in R when U crosses zero. Such
behavior is not permission to choose whichever constant velocity is closer.

## 6. The winding-projector derivative is explicit and cannot be discarded

For each geometry write `mu_g=E[q], e_g=E[q^2]`. Since q^3=q and r=q+1,
the response to a separate r tilt has the exact unmarked interface

```text
j_q,g^(r)=Cov(q,r)=e_g-mu_g^2,
j_E,g^(r)=Cov(E,r)=mu_g(1-e_g).
```

Let `D=Q_p, b=Y_p/D`. The complete moving-root response is

```text
V_r=A_N/D partial_p{
       P4[mu_g(1-e_g)]-b mean[e_g-mu_g^2] } at Q=0.
```

This includes the p derivative of b, root shift and denominator motion.
The actual m^(-r) factor contributes **-V_r** to V_Sstar, or
`-V_r/2` to the Q derivative at Q=1. Its inputs are already q/E profiles;
it is not an unknown normalization constant. Pooled mean(mu_g)=0 does
not set each mu_g to zero or remove its rank variance.

Under a proved single-field continuation, a size-independent analytic
projector normalization can be absorbed into mathcal A(t), and that
amplitude derivative cancels across sizes. But deleting m^(-r) before
establishing this can change the selected field, the root scaling section,
or its logarithmic/activated residue. Equality of bulk pressure does not
establish equality of these topological derivatives. The local Q_c^B
term likewise remains part of the specified continuum-continuation question.

## 7. What can be excluded, and the limit order

The active ordinary hypotheses are now: regular thermal logarithmic
velocity; a Q-activated four-leg sqrt(N) tangent; and a specified Jordan
logarithmic mixture. The old constant four-leg velocity remains only a
typed nonzero-overlap control.

With a common critical continuation, nonzero baseline U, and controlled
first-derivative finite-size corrections, an asymptotic S_N(c) separated
from the thermal target excludes the **regular single-thermal-field**
description of this source/readout family. A nonzero corrected increment
ratio separated from sqrt(c) excludes the declared dominant activation
shape. A different typed observer with established ordinary four-leg
overlap can separately test its constant target. Agreement is compatibility,
not a proof of universality or a physical-field identification. A finite
two-size discrepancy without a correction bound is not an exclusion of an
asymptotic law, and the measurements retain their common dependence block.

The missing physical conditions are concrete: the tied positive occupation
model must reach the claimed Potts critical branch with parameter Q_c;
the pooled root must remain in its scaling window; the selected derivative
overlap and its regularity must be known; and the explicit topological
projection and field normalization must be continued consistently.

Finally Q_c=exp(2t) exceeds4 for t>log2. The real dense/FK spectral branch
`Q=4 cos^2(pi u)` used by the old velocity dictionary cannot extend there
with real u. In the existing m=2,4,8,16 readout, m=2 is already at Q_c=4,
and the other points lie beyond it. The exact positive model and its
fixed-N two-state/negative-tail results remain valid; they do not measure
the weak-Q critical velocity. The present question takes the derivative
at t=0 and then studies growing matched sizes. It must not be confused
with the earlier limit of fixed N followed by t->infinity.

The [numerical constant packet](../results/closed-source-weak-colour-predictions.json)
records both conditional velocities, their gap and c=2,4 contrasts. Its
four-leg entry assumes nonzero overlap and is subject to the typing rule
above; it is not a restored ordinary-endpoint proposal.
