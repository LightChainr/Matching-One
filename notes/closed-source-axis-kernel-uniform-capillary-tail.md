# Uniform capillary tail of the exact directed determinant

## Result

Let

```text
r=m^-1,
c=Lr=L/m -> infinity,
alpha=L^2 r^3=L^2/m^3.                                      (1)
```

For the exact directed kernel

```text
t_r(z)=sum_(d in Z) r^|d| z^d
      =(1-r^2)/[(1-rz)(1-r/z)],
K_n=[z^n] t_r(z)^L,                                          (2)
```

define the two-boundary endpoint determinant

```text
J_(1,L)=K_0^2-K_1^2.                                         (3)
```

In the joint window

```text
c -> infinity,       alpha -> 0,                              (4)
```

equivalently

```text
L^(2/3) << m << L,                                            (5)
```

the determinant has the uniform asymptotic

```text
J_(1,L)
 =exp(4L/m)/[8pi (L/m)^2]
   *[1+O(m/L)+O(L^2/m^3)].                                   (6)
```

The key point is that `K_0` and `K_1` may each have a relative Bessel error
of order `L/m^2`, while their determinant is only an order `1/c` fraction
of `K_0^2`. The cancellation-sensitive relative error is therefore

```text
c(L/m^2)=L^2/m^3=alpha.                                      (7)
```

Thus `alpha`, rather than the ordinary rare-jump mass `L/m^2`, is a clean
sufficient power gate for replacing the exact determinant by the Bessel
determinant. This note does not claim that (4) is necessary; exploiting
additional symmetry of the compound remainder may weaken it.

Under the same uniform two-cloud and sector-odds assumptions used for the
bounded-c law, (6) gives

```text
Ustar/A_N
 =-[exp(4L/m)/(8pi Delta)] m^(-(2L-1))[1+o(1)] <0.            (8)
```

Hence the original directed two-cloud response remains negative throughout
the window (5), while acquiring the explicit `exp(4L/m)` capillary factor.

## 1. Exact compound-Poisson factorization

The logarithm of (2) is

```text
L log t_r(z)
 =L log(1-r^2)
  +sum_(k>=1) lambda_k(z^k+z^(-k)),
lambda_k=L r^k/k.                                             (9)
```

Separate the `k=1` Skellam core and normalize every higher harmonic:

```text
t_r(z)^L
 =A_L exp[c(z+z^-1)] G_R(z),                                 (10)

G_R(z)
 =exp{sum_(k>=2) lambda_k[z^k+z^-k-2]},

A_L
 =(1-r^2)^L exp(2 sum_(k>=2)lambda_k).
```

`G_R` is the probability generating function of a symmetric compound
Poisson integer R: for every `k>=2`, independent jumps `+k` and `-k` have
rate `lambda_k`. Its total event intensity is

```text
Lambda
 =2 sum_(k>=2)lambda_k
 =2L[-log(1-r)-r]
 =Lr^2[1+O(r)].                                               (11)
```

The common amplitude in (10) is much closer to one:

```text
log A_L
 =L[log((1+r)/(1-r))-2r]
 =2L(atanh(r)-r)
 =(2/3)Lr^3+O(Lr^5).                                         (12)
```

Since

```text
[z^n] exp[c(z+z^-1)]=I_n(2c),
```

equation (10) yields the exact representation

```text
K_n=A_L E[I_(n-R)(2c)].                                      (13)
```

This representation isolates both possible errors: a common amplitude
`A_L`, and a higher-jump convolution R.

## 2. Individual coefficients versus the determinant

For integer j and positive c,

```text
0<=I_j(2c)<=I_0(2c).
```

Also `P(R!=0)<=Lambda`. Consequently, for `n=0,1`,

```text
K_n=A_L[I_n(2c)+e_n],
|e_n|<=2 Lambda I_0(2c).                                     (14)
```

For `c>=1`, `I_1(2c)` is uniformly comparable with `I_0(2c)`, so (14)
also gives the familiar individual relative estimate

```text
K_n/[A_L I_n(2c)]=1+O(L/m^2),       n=0,1.                   (15)
```

That estimate is not yet sufficient for a cancellation-sensitive
determinant. Substitution into (3) gives

```text
J_(1,L)
 =A_L^2{I_0(2c)^2-I_1(2c)^2
          +O[Lambda I_0(2c)^2]}.                             (16)
```

The standard large-c Bessel ratio is

```text
1-[I_1(2c)/I_0(2c)]^2
 =1/(2c)+O(c^-2).                                             (17)
```

Therefore the error in (16), relative to the Bessel determinant, is

```text
O(c Lambda)
 =O[c(L/m^2)]
 =O(alpha).                                                   (18)
```

This proves that `alpha->0` makes the Bessel approximation relative-error
accurate at the level of `J_(1,L)`, not merely at the level of each K.
Equation (12) contributes only

```text
A_L^2=1+O(L/m^3),                                             (19)
```

which is already `1+o(1)` under (4).

The estimate (18) is deliberately robust: it uses only positivity and the
total higher-jump intensity, without assuming cancellation between `e_0`
and `e_1`. It is why `m>>L^(2/3)` is a transparent sufficient gate.

## 3. Large-c evaluation and exact sign

The Bessel asymptotics give

```text
I_0(2c)^2-I_1(2c)^2
 =exp(4c)/(8pi c^2)[1+O(c^-1)].                               (20)
```

Combining (16)-(20) proves (6).

The exact directed determinant is positive even before taking the joint
limit. Indeed, using the Fourier coefficient representation of (2),

```text
K_0-K_1
 =(1/2pi) integral_(-pi)^pi t_r(e^(i theta))^L
                         [1-cos(theta)]dtheta >0,

K_0+K_1
 =(1/2pi) integral_(-pi)^pi t_r(e^(i theta))^L
                         [1+cos(theta)]dtheta >0.              (21)
```

Thus `J_(1,L)>0` for every `0<r<1` and finite L. The role of the gate (4)
is not to manufacture positivity; it makes the explicit Bessel magnitude
in (6) uniform enough to dominate the cancellation scale and to propagate
the fixed-L original-U expansion.

## 3a. Fixed-m local CLT and the directed proliferation threshold

The same exact kernel has a distinct fixed-m limit which cannot be obtained
by sending c to infinity in (6) while discarding r. Normalize one step by

```text
phi_r(z)=t_r(z)/t_r(1),
t_r(1)=(1+r)/(1-r).
```

Its mean is zero and its variance is

```text
sigma_r^2=2r/(1-r)^2.                                        (21a)
```

If `p_L(n)` is the mass at n of the L-step normalized walk, then

```text
K_n=t_r(1)^L p_L(n).
```

The analytic lattice saddle at theta=0 gives, for fixed n and fixed
`0<r<1`,

```text
p_L(0)=[2pi L sigma_r^2]^-1/2[1+O(L^-1)],
p_L(1)/p_L(0)
 =1-[2L sigma_r^2]^-1+O(L^-2).                               (21b)
```

The ratio estimate in (21b), one order sharper than an ordinary local CLT,
is needed because the determinant subtracts two neighboring masses. It gives

```text
J_(1,L)
 =t_r(1)^(2L)/[2pi L^2 sigma_r^4]
   *[1+O(L^-1)]
 =[(1-r)^4/(8pi r^2 L^2)]
   [(1+r)/(1-r)]^(2L)[1+O(L^-1)].                            (21c)
```

Through the two-path nonintersection determinant, (21c) is a positive
noncrossing-path partition, not a formal signed cancellation.

Restore the bare two-boundary winding barrier `m^(-2L)`, with `r=1/m`.
Its exact exponential base is

```text
q_dir(m)
 ={[(1+r)/(1-r)]/m}^2
 =[(m+1)/(m(m-1))]^2.                                       (21d)
```

The unique positive threshold is

```text
q_dir(m_dir)=1,
m_dir^2-2m_dir-1=0,
m_dir=1+sqrt(2).                                              (21e)
```

Consequently:

- if `1<m<m_dir`, this positive directed rank-one endpoint subfamily already
  grows exponentially after the nominal `m^(-2L)` barrier. That is enough to
  invalidate a simple two-phase suppression argument based only on the
  minimum winding cost;
- if `m>m_dir`, this directed subfamily is exponentially suppressed, but the
  conclusion controls neither overhangs nor other non-directed rank-one
  families and is not a full phase-transition theorem;
- at `m=m_dir`, the endpoint determinant with its barrier has the explicit
  polynomial law

```text
m_dir^(-2L) J_(1,L)
 =[(3-2sqrt(2))/(2pi)] L^-2[1+O(L^-1)].                      (21f)
```

Equation (21f) is the single endpoint-channel prefactor. Translation,
orientation and other lattice multiplicities must be restored separately;
they must not be silently absorbed into this local CLT.

## 4. Propagation to original U

The matched two-cloud capillary law has the finite-kernel form

```text
Ustar/A_N
 =-[L^2 m^(-(2L+1))/Delta] J_(1,L)
   *[1+epsilon_sector(L,m)].                                  (22)
```

Here `epsilon_sector` collects the already-declared non-directed,
additional-component, tilted-sector and odds-alignment remainder; it is not
controlled by the one-dimensional kernel alone. If

```text
epsilon_sector=o(1)
```

uniformly in (5), equations (6) and (22) give

```text
L^2/(L/m)^2=m^2
```

and hence (8). Neither root motion nor the positive within-geometry
denominator changes the sign. The relative determinant error is `O(alpha)`,
so it cannot overturn the leading negative term when `alpha->0`.

The logarithm of the magnitude is

```text
-(2L-1)log m+4L/m+O(1),                                      (23)
```

which makes the division of roles explicit: `2L log m` remains the winding
tension, while `4L/m` is the resummed capillary entropy.

## Scientific card

- **Exact bridge:** equation (13) expresses the finite directed kernel as a
  Skellam/Bessel core convolved with an explicit higher-jump compound
  Poisson remainder.
- **Power gate:** the ordinary remainder mass is `L/m^2`, but determinant
  cancellation amplifies it by `c=L/m`; the robust relative-error parameter
  is `alpha=L^2/m^3`.
- **Uniform law:** for `L^(2/3)<<m<<L`, equation (6) supplies the explicit
  `exp(4L/m)/(8pi c^2)` determinant tail.
- **Fixed-m threshold:** the exact directed exponential base is (21d), with
  `m_dir=1+sqrt(2)` and the critical endpoint prefactor (21f).
- **Signed consequence:** conditional on the stated uniform sector
  remainder, original U obeys (8) and remains negative.
- **Boundary:** `alpha->0` is a sufficient positivity-preserving Bessel gate,
  not a proof of its optimality. Likewise, `m>m_dir` suppresses only the
  directed endpoint family and does not control a complete lattice phase.
