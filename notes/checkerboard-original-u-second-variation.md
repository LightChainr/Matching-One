# A genuine Bernoulli coupling path: exact second transmission to original U

This is a finite-model derivation, with no sampling, source search, fitting,
fresh P154 readout or change to the lag1 source's stopping decision.
It starts from `4bb851cb`. The related older
[`c4-self-matching-tangent-improvement.md`](c4-self-matching-tangent-improvement.md)
already uses two sublattice probabilities on a **different graph**, the
degree8/degree4 checkerboard triangulation. There a unit translation does
not exchange equivalent sublattices and a staggered first derivative can
be nonzero. Here the graph is ordinary degree4 square site percolation;
the occupation law changes while that graph and the rank readouts stay fixed.

## 1. Legal quotients and the even transverse chart

Both period vectors must have even coordinate sum. Then vertex parity
descends to the torus, a unit translation exchanges colors A/B, and each
color has N/2 sites. Merely taking even N without this period condition
is insufficient. Use the same legal condition in every orientation of P4.

First write p_A=p+h, p_B=p-h, with |h|<min(p,1-p). For each geometry the
configuration readouts q=rank-1 and E=1{rank=0 or2} are translation
invariant. Therefore their expectations F(p,h) obey F(p,h)=F(p,-h).
At h=0, F_h=F_ph=0, including both the pooled Q=mean(q) and Y=P4(E).

Let p0 be Q(p0,0)=0 and, throughout this note,

```
A_N=N^(13/8)/2, D=Q_p, T=Q_pp, B=Y_p, r=B/D,
U=A_N Y_p/D evaluated on the corresponding pooled-root branch.
```

P4 and mean have their original fixed geometry weights. An h subscript
means an ordinary derivative, not the coefficient of h² in a Taylor series.

## 2. The correct Bernoulli second score

For a configuration set K=K_A+K_B, Delta=K_A-K_B, u=p(1-p), and
V=K-Np. At h=0 the log-density derivatives and density second score are

```
ell_h = Delta/u,
ell_hh = -K/p² -(N-K)/(1-p)²,
H = ell_hh+ell_h²
  = [Delta² -(1-2p)K -Np²]/[p²(1-p)²]
  = [Delta² -Np(1-p) -(1-2p)(K-Np)]/u².
```

The final K-dependent centering term must not be omitted away from p=1/2.
At p=1/2, H=16Delta²-4N. Exactly E[H]=0 and, for any fixed configuration
observable O, F_hh=E[O H]. Translation makes the first derivative zero;
it does not force this even second response to vanish.

Thermal derivatives must also differentiate H. The needed mixed score is

```
S_p=V/u,
H_p=2V/u² -2(1-2p)H/u,
J=H S_p+H_p,
F_phh=E[O J].
```

Thus the complete observable interface is
Q_hh=mean E[qH], Y_hh=P4 E[EH], Q_phh=mean E[qJ], Y_phh=P4 E[EJ].
H and J can also be used as centered covariance weights because their
unweighted expectations vanish. No unknown coupling coefficient enters.
Ordinary K-only q/E histograms do not in general determine the necessary
Delta² correlations; this note does not invent them from existing summaries.

## 3. Moving root and slope: the full transverse second derivative

Evenness and D>0 give

```
p_h(0)=0,  gamma:=p_hh(0)=-Q_hh/D,
(Y_p along the root)_hh=Y_phh+gamma Y_pp,
(Q_p along the root)_hh=Q_phh+gamma T.
```

Both first root-comoving derivatives vanish, so the complete ratio rule is

```
U_hh = A_N/D [Y_phh+gamma Y_pp-r(Q_phh+gamma T)]
     = A_N/D partial_p [Y_hh-r Q_hh] at p0.
```

The p derivative acts on r as well. Equivalently, if
L(j_Q,j_Y)=A_N/D partial_p(j_Y-r j_Q), then U_hh=L(Q_hh,Y_hh).
This is the same original global U, including root and denominator motion;
it is not an unnormalized angular second moment.

## 4. The endpoint-compatible black-saturation family

For the finite path use **p_A=s+(1-s)p, p_B=p**, 0<=s<=1.
At s=1, partial_p p_A=0 and partial_p p_B=1: the thermal direction is
solely the remaining white-color probability. Fixed-h differentiation at
the p+h=1 edge would not have this property and cannot justify a child-U map.

Put m=(p_A+p_B)/2=p+s(1-p)/2. In the inverse thermal chart,

```
h = [s/(2-s)](1-m),
F_tilde(m,s)=F(m,h(m,s)),
F_tilde_s|0=0,
F_tilde_ss|0=c(m) F_hh(m,0),  c(m)=(1-m)²/4.
```

The thermal Jacobian dm/dp=1-s/2 multiplies numerator and denominator
of U equally, so this chart leaves original U exactly unchanged. Consequently

```
U_s(0)=0,
U_ss(0)=L(c Q_hh,c Y_hh)
       = A_N/D {c[Y_phh-r Q_phh-r'Q_hh]
                 +c'[Y_hh-r Q_hh]},
c=(1-p0)²/4, c'=-(1-p0)/2.
```

In particular **U_ss=c U_hh+(A_N/D)c'(Y_hh-rQ_hh)**, not c U_hh alone.
The temperature dependence of the physical saturation dose is part of the
model. Even a transverse change invisible to U_hh at this root can be
visible to this particular s path through its nonzero transverse value.

The directly measurable saturation weights simplify to

```
W=cH=[Delta²-(1-2p)K-Np²]/(4p²),
Z=W S_p+partial_p W=W S_p+V/(2p²)-2W/p.
```

Set j_Q=mean E[qW], j_Y=P4 E[EW] and replace W by Z for their p
derivatives. Then U_ss=A_N/D*(j_Y,p-r j_Q,p-r'j_Q). These translation-even
weights already remove the vanishing odd part of the transformed-family
second score and retain the c' term exactly; no finite difference is needed.

For an equivalent check in the original p coordinate, write w=(1-p)/2.
At s=0, F_s=wF_p and F_ss=w²(F_pp+F_hh). The original pooled root obeys

```
p_s=-w,
p_ss=2ww'-w²Q_hh/D=-(1-p)/2-c Q_hh/D.
```

The w² F_pp, first-order root motion and first-order slope motion cancel
as a common thermal reparameterization; the remaining expression is the
boxed c-dependent second transmission above. Thus a nonzero first-order
source score on configurations must not be reported as nonzero U_s.

## 5. Exact null, finite endpoint, and what their combination implies

A thermal-only hypothesis has every q_g,E_g on a common temperature orbit
F_s(p)=F_0(phi_s(p)); it predicts U(s) constant exactly. Infinitesimally
Y_hh-rQ_hh vanishes if the second transverse jet is only b(p)F_p.
A nonzero measured U_ss excludes that thermal-only second-jet explanation.
A zero U_ss does not prove it: a transverse profile can have a stationary
slope at the root or cancellations in the displayed formula. There is no
claim of a new continuum field or a universal sign.

Given the independently proved endpoint topology, write r_child=1-p at
s=1. Its complement and quarter-diagonal rotation imply

```
q_parent(p,1)=-q_child(1-p), E_parent(p,1)=E_child(1-p),
delta_cos4_child=-delta_cos4_parent,
U_parent(1)=2^(13/8) U_child(N/2).
```

The complement reverses the E thermal derivative; rotation reverses the
angular denominator. Those signs cancel. This note supplies the thermal
direction and normalization check, not a replacement proof of endpoint
topology. Comparisons to a native parent U must retain the original shared
covariance when parent/child anchors are dependent.

For these proper square quotients the saturation path also has a regular
root branch. At p=0 all occupied vertices are in A and have no edges, so
Q=-1; at p=1 the torus has rank2, so Q=1. The rank readout is increasing in
occupation. Conditional on all A vertices occupied (positive probability
for p>0), it is a nonconstant increasing function of B, hence at least one
B-site influence is positive under the interior product measure. Since
partial_p p_B=1 and partial_p p_A>=0, Q_p>0 throughout 0<p<1, including
s=1. Polynomial dependence and the implicit-function theorem give the
unique analytic root branch and smooth U(s) on this finite interval.

Since U_s(0)=0, exactly

```
Delta U=U(1)-U(0)=integral_0^1 (1-s) U_ss(s) ds.
```

Thus a positive finite endpoint contrast and zero initial first derivative
are fully compatible; a positive contrast forces positive curvature
somewhere, not necessarily at s=0. It does not justify quadratic
interpolation, a guessed all-s curve, or a prediction that U_ss(0)>0.
The finite coupling family is independent of the P154 lag1 intervention;
its endpoint information does not rescind that source's recorded stop line.
