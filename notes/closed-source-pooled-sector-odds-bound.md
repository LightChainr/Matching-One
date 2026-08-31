# A sharp sector-odds bound for the original pooled denominator

**Result.** At the original equal-weight two-geometry pooled root, let
`b=(P1_f+P1_s)/2` and suppose the two log odds of rank2 versus rank0 differ
by at most delta. Then the within-geometry denominator obeys the sharp bound

```text
kappa >= b(1-b)+(1-b)^2 sech^2(delta/4).                         (1)
```

The bound is sharp given only b and delta, even with nonzero rank-one
probability. A sharper exact formula below uses the two individual P1
values. This converts the denominator problem into a specified partition
cross-ratio; it does not assert a new bound on that cross-ratio for the
original tilted lattice. This note starts from `e17b286b` and uses no
sampling, numerical re-score or cross-geometry mixture covariance.

## 1. Exact sector parametrization on the pooled root

For geometry g=f,s let Z_j,g be its positive restricted partition function
in ambient rank j=0,1,2 under the same source and thermal parameter.
Configuration-independent normalization factors do not matter. Define

```text
Pj_g=Z_j,g/(Z_0,g+Z_1,g+Z_2,g),
w_g=1-P1_g,  ell_g=log(Z_2,g/Z_0,g),
q_g=P2_g-P0_g=w_g tanh(ell_g/2).
```

At the equal-weight pooled root, q_f=a and q_s=-a. Hence

```text
kappa=(Var_f(q)+Var_s(q))/2
     =(w_f+w_s)/2-a^2=1-b-a^2.                                (2)
```

For a!=0 the two ell values have opposite signs. Therefore, putting
d=|ell_f-ell_s| and u=|a|,

```text
d=2 atanh(u/w_f)+2 atanh(u/w_s).                               (3)
```

If a=0 then ell_f=ell_s=0 and the same relation holds. All finite positive
laws have w_f,w_s>0 and u<min(w_f,w_s). Boundary distributions can be
handled by limits.

The observable quantity controlling this mismatch is exactly

```text
Xi = (Z_2,f Z_0,s)/(Z_0,f Z_2,s),  d=|log Xi|.                 (4)
```

It is evaluated at the common pooled root. It involves separately
restricted partition functions, not either geometry's total pressure.

## 2. Exact inversion, and the optimal bound with the individual P1 values

Set `W=w_f+w_s, P=w_f w_s, rho=tanh(d/2)`. Applying the addition formula
to (3) gives `rho=u W/(P+u^2)`. The physical solution is

```text
u = 2 rho P/[W+sqrt(W^2-4 rho^2 P)],
kappa = W/2-u^2.                                              (5)
```

This form is nonsingular at d=0, where its numerator is zero. If only
`d<=delta` is known, use `rho_delta=tanh(delta/2)` in (5) to obtain an
upper bound u_delta on u and the lower bound `W/2-u_delta^2` on kappa.
The inverse in (3) is monotone, so this is the best possible bound at
fixed w_f,w_s and delta among normalized sector distributions.

To attain it, choose `a=u_delta` and

```text
(P0_f,P1_f,P2_f)=((w_f-a)/2,1-w_f,(w_f+a)/2),
(P0_s,P1_s,P2_s)=((w_s+a)/2,1-w_s,(w_s-a)/2).
```

They have the prescribed pooled zero and odds mismatch delta. Sharpness
of the probability bound is not a claim that every extremizer is realized
by the closed source on a square quotient.

## 3. Sharp mean-only bound, including nonzero rank-one mass

For fixed u>0, the function `w -> atanh(u/w)` is convex on w>u, since
its second derivative is `2u w/(w^2-u^2)^2`. Jensen's inequality in (3)
therefore gives, with w=(w_f+w_s)/2=1-b,

```text
delta>=d>=4 atanh(u/w),  u<=w tanh(delta/4).
```

Substitution into (2) proves

```text
kappa >= w-w^2 tanh^2(delta/4)
       = b(1-b)+(1-b)^2 sech^2(delta/4).
```

Equality is attained by w_f=w_s=w and ell_f=delta/2, ell_s=-delta/2,
so no larger uniform lower bound follows from b and delta alone.
With the actual cross-ratio Xi, the same result is directly computable as

```text
kappa >= b(1-b)+(1-b)^2 *4 sqrt(Xi)/(1+sqrt(Xi))^2.             (6)
```

In particular d=0 makes each geometry individually balanced and kappa=1-b.
For b=0, equation(1) becomes the exact two-state formula
`kappa=sech^2(d/4)`. If b->0 while d stays bounded, kappa is bounded away
from zero; if both b,d->0, kappa->1.

If only `0<=b<=epsilon<1` is available, one must not simply replace b by
epsilon in (1), since its right side is not monotone in b. Writing
`sigma=sech^2(delta/4)`, the sharp uniform bound over that interval is

```text
kappa >= min{ sigma,
              (1-epsilon)[epsilon+(1-epsilon)sigma] }.         (7)
```

This follows because `w-(1-sigma)w^2` is concave on
w in[1-epsilon,1]. The convenient weaker bound `(1-epsilon)sigma` is
always valid.

## 4. Direct consequence for the existing original-U bound

The [fixed-coupling denominator reduction](closed-source-fixed-coupling-peierls.md)
already gives, with the original nonzero angular denominator Delta,

```text
|U/A_N| <= 4N b/(|Delta| kappa),  A_N=N^(13/8)/2.
```

Equation(1) makes this explicit without replacing a within-geometry
variance by a mixture variance:

```text
|U/A_N| <= 4N b/
          {|Delta|(1-b)[b+(1-b)sech^2(delta/4)]}
       <= 4N b cosh^2(delta/4)/{|Delta|(1-b)}.                 (8)
```

Thus a sufficient condition for the full U, including its area factor,
to vanish along a size sequence is

```text
N^(21/8) b cosh^2(delta/4)/[|Delta|(1-b)] -> 0.                (9)
```

For example, if N is proportional to L^2, Delta stays nonzero, and
`b<=exp[-tau L+o(L)]`, a separately proved `delta=o(L)` closes (9).
More generally `delta<=gamma L+o(L)` is sufficient when gamma<2tau.
The factor1/2 in this competition follows from
`cosh^2(delta/4)<=exp(delta/2)`; it must not be replaced by an unexplained
claim that any surface-order mismatch is harmless.

## 5. What is still not implied by same area or equal pressure

Even identical total partition functions do not bound d. For any
`0<=b<1` and d>=0, consider positive sector weights

```text
Z0_f=exp(-d/4),  Z2_f=exp(d/4),
Z0_s=exp(d/4),   Z2_s=exp(-d/4),
Z1_f=Z1_s=[2b/(1-b)]cosh(d/4).
```

Both total partition functions are exactly `2 cosh(d/4)/(1-b)`, the
pooled q is zero, and kappa attains (1) with mismatch d. For b=0 interpret
the rank-one weight by its limiting value; b>0 makes all sectors positive.
This is a logical counterexample to a pressure-only inference, not a
claimed realization by the specified lattice action. Choosing d of order
L already makes kappa potentially exponentially small while the two
pressure densities remain exactly equal.

An actual rank-preserving, weight-preserving bijection between the two
geometries would give Xi=1. Equal area and the same named source do not
supply such a bijection for the original axis/tilted pair. Nor do the
current action identities prove a sublinear bound on their restricted
sector free-energy mismatch. One needs a separate estimate of (4), or
equivalent control of the rank0 and rank2 partition ratios.

The mathematical closure delivered here is therefore exact and limited:
the original pooled denominator is now a sharp function of rank-one mass
and a single specified partition cross-ratio. The tilted winding bound
and the needed fixed-t control of that cross-ratio are not asserted to
have been proved by this note. Any independently established regime with
b,d->0 immediately yields kappa->1 through the same formula.

The separate [Poisson double-scaling theorem](closed-source-poisson-double-scaling.md)
now supplies such a regime for the actual two fixed lattice laws:
`N/m² -> zeta < infinity` with growing systole gives b,d->0, including
the oblique companion. The fixed-m problem remains distinct.
