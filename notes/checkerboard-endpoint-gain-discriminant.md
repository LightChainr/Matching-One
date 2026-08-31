# Does the closed source acquire only a geometric gain near saturation?

Pure finite-model derivation from `0d19179f`; no new numerical outcome,
sampling, defect enumeration, gain fit or F4 re-score is used here.
The source remains exactly `S_star=C+F+Bv`. The probabilities are
`p_A=s+(1-s)p, p_B=p`, and every thermal derivative below holds `(s,t)`
fixed in the positive law proportional to `P_(p,s) exp(t S_star)`.

## 1. The scalar hypothesis and its coefficient-free discriminator

Normalize each geometry's expectation separately, then use the original
fixed mean and P4 weights to form `Q=mean(q)` and `Y=P4(E)`. Write

```text
Q(p0(s,t),s,t)=0,  D=Q_p,  A_N=N^(13/8)/2,
U(s,t)=A_N Y_p/D at p0(s,t),  V=U_t.
```

Use a simple-root branch with `D!=0` and `U!=0` near `(s,t)=(1,0)`.
Finite partition functions are analytic there; the physical s derivative
at1 is from the left. No regularity of every finite-t root is assumed.
The established endpoint map gives, at the **same** source strength t,

```text
U(1,t)=gamma U_child(t),  gamma=2^(13/8).
```

The new hypothesis is stronger: a source-independent gain explains the
first departure from saturation,

```text
H_gain: U(1+delta,t)
       =gamma U_child(t)+delta g1 U_child(t)+O(delta^2).
```

It predicts the determinant

```text
R = U U_st-U_s U_t = U^2 partial_s(V/U) = 0                 (s=1).
  = gamma [U_child U_st-V_child U_s].
```

Thus only the relative source transmission `V/U` is being compared.
The unknown g1 cancels: no gain is fitted, and no division by V is needed.
A nonzero R at t=0 rejects this multiplicative-gain explanation of the
mixed first jet of the original U. It does not undo the exact endpoint
identity or the configurationwise closure of S_star.

Zero R at one t is only a necessary first-jet condition. For example
`U=gamma u(t)+delta[a u(t)+b t^2]`, with `u(0)!=0, b!=0`, passes R(1,0)=0
but fails H_gain as a function of t. Even an exact scalar U factorization
would not establish equality of the complete q/E profiles or measures.

## 2. One-defect partition data supply the entire required s jet

Let `M=N/2`. For one geometry and fixed observable O, let `H0_O(p,t)`
be its unnormalized moment with all A sites occupied. Let `Hd_O(p,t)`
be the average of the corresponding moments with exactly one uniformly
chosen A site vacant. In both cases B is independent Bernoulli(p), the
observable is evaluated on the **parent**, and the weight is the actual
parent `exp(t S_star)`. In particular, a defective parent's source must
not be replaced by the no-defect child source using endpoint closure.
Translation-invariant O and source permit an equivalent representative
A defect when translation equivalence has been established.

With `epsilon=1-s`, the A-vacancy probability is `epsilon(1-p)`, so

```text
H_O(p,s,t)=H0_O+epsilon M(1-p)(Hd_O-H0_O)+O(epsilon^2).
Z=H_1,  F_O=H_O/Z.
```

No two-defect sector is needed for U_st. For each O in `{1,q,E}`, only
the following four ordinary functions of p are required:

```text
P_O = H0_O(p,0),                   T_O = H0_(O S_star)(p,0),
A_O = M(1-p)[Hd_O(p,0)-P_O],
B_O = M(1-p)[Hd_(O S_star)(p,0)-T_O].
```

Their two-parameter jet is `H_O=P_O+t T_O+epsilon A_O+epsilon t B_O`.
Here `P_1=1, A_1=0`, but generally `T_1!=0, B_1!=0`. In particular
`Z_t=T_1` and `Z_st=-B_1` are not zero. Exact division gives

```text
F_O   =P_O,
F_O,t =T_O-P_O T_1,
F_O,s =-A_O,
F_O,st=-B_O+A_O T_1+P_O B_1                         at s=1,t=0.
```

Differentiate these complete functions in p for the thermal jets. The
derivative also acts on the explicit `(1-p)` factor, T_1 and B_1. More
generally, if `Delta H=Hd-H0`,

```text
partial_p^i H_s=-M[(1-p)partial_p^i Delta H
                         -i partial_p^(i-1) Delta H],
```

with the second term omitted for i=0; the same rule holds for H_st
after replacing H by its first source moment.

A useful interpretation of the mixed normalization correction is

```text
F_O,st = -M(1-p) [ Cov_d(O,S_star)-Cov_0(O,S_star)
                  +(E_d O-E_0 O)(E_d S_star-E_0 S_star) ].
```

The defect expectation here includes the uniform choice of missing A
site. A bare difference of the two normalized covariances omits the final
product: the source also changes the relative statistical weight of the
defect sector. Each geometry has its own partition denominator; pooling
H_q or H_E before division changes the model.

## 3. Complete moving-root and slope mixed derivatives

All partial derivatives in this section are those of the **normalized**
Q and Y before substitution of the root. Put `r=Y_p/D`. For a=s,t,

```text
p_a=-Q_a/D,
p_st=-(Q_st+Q_ps p_t+Q_pt p_s+Q_pp p_s p_t)/D.
```

For either `L=Y` or `L=Q`, define the root-comoving slope jets

```text
J_L   = L_p,
J_L,a = L_pa+L_pp p_a,
J_L,st= L_pst+L_pps p_t+L_ppt p_s
                    +L_ppp p_s p_t+L_pp p_st.
```

Thus J_Q=D and its comma derivatives include root motion. The complete
ratio rules are

```text
U_a = A_N/D [J_Y,a-r J_Q,a],

U_st= A_N/D [J_Y,st-r J_Q,st
       -(J_Q,s/D)(J_Y,t-r J_Q,t)
       -(J_Q,t/D)(J_Y,s-r J_Q,s)].
```

These include denominator motion, root p_st and the third thermal jets.
Substitute them directly into R; there is no assumption of zero U_s,
zero partition-source mean, or a rigid translation of the thermal root.

For a minimal exact/arbitrary-precision calculator, the required p orders
of Q/Y are: base through3, s and t through2, and st through1. Equivalently:

1. Form H and Z as jets in `delta=s-1`, t and `xi=p-p00`, then divide.
   The signs of the defect terms are `-delta A_O-delta t B_O`.
2. Pool normalized Q/Y, solve `Q(p00,1,0)=0`, and insert
   `p=p00+p_s delta+p_t t+p_st delta t` into Q. Set its three nonconstant
   coefficients to zero. Mixed Taylor coefficients have no factor1/2.
3. Differentiate Q/Y in p **before** that root substitution, divide their
   slopes, and multiply by A_N. The delta, t, delta*t coefficients are
   U_s,U_t,U_st. Return `U0*U_st-U_s*U_t`.

Only first source moments are required: neither S_star^2 nor a second
defect sector enters this mixed first jet. At a nonzero source reference
t0 the same scheme uses moments weighted by exp(t0 S_star), with general
partition division rather than assuming Z=1 or Z_s=0.

## 4. Quotient meaning and reparameterization boundary

Use the common thermal coordinate `x=Q(p,s,t)` and write
`mathcal Y(x,s,t)=Y(p(x,s,t),s,t)`. Then

```text
U=A_N partial_x mathcal Y(0,s,t),
zeta_a=Y_a-r Q_a,
U_a=A_N/D partial_p zeta_a at the root.
```

This proves invariance under any common invertible thermal relabeling
`p=phi(p_tilde,s,t)` that leaves s and t fixed. The relabeling may depend
on s,t, and its Jacobian may be negative; both thermal slopes acquire
the same Jacobian. It must be common to all orientations and observables.

The zero of R also survives separate nonsingular changes of s and t:
R acquires their two chain-rule Jacobians. Multiplying U by separate
nonzero factors `a(s)b(t)` multiplies R by `a(s)^2 b(t)^2`.
It is **not** invariant under mixing the two physical parameters. For
example `t=a(s) tau`, with `a(1)=1`, gives at tau=0

```text
R_new=R+a'(1) U U_t.
```

Consequently an s-dependent rescaling of the source can manufacture or
remove a mixed gain residual. The exact dictionary fixes the coefficient
of the named S_star to the same t, so that freedom is not part of H_gain.
Likewise, a saturation coordinate that depends on p changes the fixed-s
thermal direction and cannot be treated as a harmless common p clock.

## 5. Density-gauge caution specific to this two-sublattice family

Dropping the `-4K` term from the action was a common-logit relabeling in
the homogeneous family, and remains one on the saturated endpoint. It
does not automatically preserve the present endpoint **mixed** derivative.
For the source change `S_star -> S_star+cK`, let
`f_l(z)=logistic(logit(z)+l)`. The exact measure relabeling is

```text
p'=f_(ct)(p),
s'=[f_(ct)(p+s(1-p))-f_(ct)(p)]/[1-f_(ct)(p)],
p'=p+ct p(1-p)+O(t^2),
s'=s+ct s(1-s)(1-p)+O(t^2).
```

Although s'=1 at s=1, this relabeling mixes saturation and temperature
nearby. Its first source tangent on any normalized observable is

```text
j_F^(cK)=c p(1-p) F_p+b F_s,  b=c s(1-s)(1-p).
```

The first term is a common thermal clock and cancels in U. At the endpoint
b=b_p=0, so adding cK leaves U_t unchanged there. However b_s=-c(1-p)
and b_ps=c. With `zeta_s=Y_s-r Q_s`, the exact t=0 differences are

```text
Delta U_st = c [A_N zeta_s/D-(1-p)U_s],
Delta R    = c U [A_N zeta_s/D-(1-p)U_s]                 at s=1.
```

Removing `-4K` means c=4; the constant2N has no normalized effect.
Thus the endpoint equality of U(t) does not license exchanging these
two sources for this R while retaining the same linear saturation family.
No new observable is introduced by this caution: zeta_s and U_s are
already supplied by the one-defect/root jets above.

## Scientific scope

The new target is a one-defect, first-source-jet obstruction to a
source-independent geometric gain. Its calculation has a finite named
input, all normalization terms and an explicit nuisance-free determinant.
A nonzero R would identify coupling-changing behavior of this observer
immediately off saturation, beyond its exact endpoint gain. Zero leaves
higher source jets and other profile directions unresolved. Neither
outcome alone identifies a continuum field, an RG flow, or asymptotic
amplitudes, and neither changes the completed F4-only stopping decision.
