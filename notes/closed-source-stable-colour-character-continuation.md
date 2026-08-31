# Stable colour character and a fixed occupation continuation of the seam packet

For the declared one-seam torus trace, the generic `[Q-2,2]` character
coefficient can be obtained exactly without enumerating colours. The
result also states when a finite-Q seam value is outside the stable
range. Under the separate N25 packing input stated in Section 4, the
actual trace at **every integer Q>=4** reduces to two named configuration
classes, with rational weights `(Q-3)/(2Q)` and `(Q-3)/2`. This fixes an
explicit occupation-based continuation to Q=1; no interpolation of a
few colour numbers or additional coupling calculation is used.

Source: `0dda27ba`, in particular
[the finite torus closure](closed-source-finite-torus-pair-closure.md)
and [the S4 rank-one filter](closed-source-s4-rank-one-filter.md).
The seam, primitive deck basis, original q/E, and source Sstar are
unchanged. The packing statement below is proved in the separate
[N25 winding packing note](n25-winding-packing-and-pair-continuation.md).
This note supplies its character-theoretic consequence.

## 1. Finite character coefficient and its stable range

Consider a rank-one configuration with w>=1 essential components of
common primitive winding `(u,v)` and c0 contractible hypergraph blocks.
Under seam permutation pi its colour factor is

```
Q^c0 Fix(pi^u)^w.
```

Let X_d be the number of d-cycles of pi. For n=|u|>=1,

```
Y_n(pi)=Fix(pi^u)=sum_(d divides n) d X_d(pi).
```

At integer Q>=4, the irreducible character and its dimension are

```
chi_[Q-2,2] = binom(X1,2)+X2-X1,
d2(Q)=Q(Q-3)/2.
```

Define the finite coefficient by the class-function inner product

```
kappa_Q(u,w) = E_(pi uniform in S_Q)[chi_[Q-2,2](pi) Y_n(pi)^w]. (1)
```

The projected configuration colour factor is therefore
`d2(Q) kappa_Q(u,w) Q^c0`. For u other than +/-1 the power operation
`pi -> pi^u` need not define a representation action. Thus “multiplicity”
in (1) means a character coefficient and need not be nonnegative at
small Q. This is compatible with a signed central transfer trace.

For nonnegative integers a_d, the exact permutation factorial moments are

```
E prod_d (X_d)_(a_d) = prod_d d^(-a_d), if sum_d d a_d <= Q,
                    = 0,              otherwise.             (2)
```

Here `(X)_a` is a falling factorial. To prove (2), choose the ordered
disjoint marked cycles and permute the remaining labels; the remaining
factorials cancel, leaving one factor 1/d for each marked d-cycle.

The weighted degree of `chi Y_n^w`, assigning degree d to X_d, is at
most `nw+2`. Expanding in falling factorials consequently proves

```
Q >= max(4, |u|w+2)  =>  kappa_Q(u,w)=kappa(u,w),               (3)
```

where the stable expectation replaces the X_d by independent Poisson
variables with means 1/d. This is an exact factorial-moment identity,
not a Poisson approximation to a finite-size error. The bound is
sufficient, not asserted optimal in every case.

If u=0, `Fix(pi^0)=Q` is constant, so orthogonality gives
`kappa_Q(0,w)=0` for every Q>=4. There is no divisor-of-zero convention.

## 2. Closed stable formula

Let

```
Y=sum_(d divides |u|) d P_d,       P_d independent Pois(1/d).
```

For a Poisson variable, marking a factorial pair of 1-cycles adds 2
to Y, and marking one 1-cycle adds 1. Marking one 2-cycle adds 2
exactly when u is even. Applying these shifts to (1) gives

```
u odd:  kappa(u,w) = (1/2) E[(Y+2)^w-2(Y+1)^w+Y^w],
u even, u!=0:
        kappa(u,w) = E[(Y+2)^w-(Y+1)^w],
u=0:    kappa(0,w) = 0.                                      (4)
```

In particular the X2 term for odd u is `(1/2)E[Y^w]`, since X2 is
independent of Y; for even u it is `(1/2)E[(Y+2)^w]`. This accounts
for the different first and second finite differences in (4).

An equivalent exact generating function, useful without any random
variables, is

```
M_n(z)=exp(sum_(d divides n) (exp(dz)-1)/d),
sum_(w>=0) kappa(u,w) z^w/w!
  = (exp(z)-1)^2 M_n(z)/2,          u odd,
  = exp(z)(exp(z)-1) M_n(z),        u even and nonzero.         (5)
```

Consequences needed here are

```
kappa(odd,1)=0,                  kappa(even nonzero,1)=1,
kappa(odd,2)=1,                  kappa(1,3)=6.                 (6)
```

For example `E[Y]=number_of_divisors(n)` and the odd w=3 value is
`3 E[Y]+3`. The stable coefficients in (4) are nonnegative, although
this does not turn a general finite-Q twisted propagation kernel into
a positive operator.

## 3. Two exact finite-Q aliases

At Q=4 the relevant class table is

| Class | Size | chi_[2,2] |
|---|---:|---:|
| identity | 1 | 2 |
| transposition | 6 | 0 |
| double transposition | 3 | 2 |
| three-cycle | 8 | -1 |
| four-cycle | 6 | 0 |

For u=3,w=1, the three contributing fixed-colour values are 4,0,4.
Hence

```
kappa_4(3,1) = [2*4 + 3*2*0 - 8*4]/24 = -1,
kappa(3,1)=0,                    stable range Q>=5.           (7)
```

For u=1,w=3, those values before cubing are 4,0,1, so

```
kappa_4(1,3) = [2*4^3 + 3*2*0^3 - 8*1^3]/24 = 5,
kappa(1,3)=6,                    stable range Q>=5.           (8)
```

These finite answers are genuine character contractions. Calling them
aliases means that the support restriction in (2) has removed terms
present in the stable polynomial; it does not mean the negative Q4
trace in (7) is a numerical or geometric error. A Q4 value cannot be
silently treated as the stable coefficient and continued to Q1.

## 4. Conditional N25 simplification: precisely two surviving classes

Use the following separate geometric packing input for each fixed N25
geometry and its declared first deck seam:

```
w<=2;
w=2 implies u=0 or |u|=1;
w=1 implies |u|<=2.                                          (9)
```

No requirement on the second primitive coordinate v is added here.
Under (9), every nonzero-u case already satisfies the sufficient range
(3) for all Q>=4. The possible coefficients are

| Class | Essential-cluster data | kappa_Q, every integer Q>=4 |
|---|---|---:|
| A | w=2, abs(u)=1 | 1 |
| B | w=1, abs(u)=2 | 1 |
| other allowed rank1 | u=0, or w=1 and abs(u)=1 | 0 |

Rank0 has constant colour character. Rank2 has exactly one essential
block carrying the full saturated image, whose seam factor is Fix(pi).
Their `[Q-2,2]` projections vanish by character orthogonality. Moreover
q=E=0 on rank1, so the **whole direct numerator remains zero**, as in
the existing finite closure.

Relative to the untwisted colour factor `Q^(c0+w)`, the exact central
partition insertion on the two surviving classes is

```
beta_Q(A)=d2(Q)/Q^2=(Q-3)/(2Q),
beta_Q(B)=d2(Q)/Q  =(Q-3)/2,
beta_Q(other)=0.                                             (10)
```

The activity, contractible colours, and original multiplier Q^(-r/2)
cancel in this ratio. At Q4 these are respectively 1/8 and 1/2, but
(10) follows from the all-integer character theorem and packing,
not from fitting those two numbers. On these two N25 geometries the
integer-Q>=4 packet is even nonnegative configurationwise, conditional
on (9); the more general rank1 negative example (7) is outside (9).

## 5. Explicit occupation completion to Q1

Let `Z_A(p,Q),Z_B(p,Q)` be the original closed-source occupation sums
restricted to the geometric classes A,B. They are not independently
renormalized. The fixed continuation of this trace packet is

```
Z_[2](p,Q) = (Q-3)/2 [Z_A(p,Q)/Q + Z_B(p,Q)],
N_q,[2](p,Q)=N_E,[2](p,Q)=0.                                 (11)
```

Equation (11) is defined for every positive real Q by the existing
occupation weights and rational coefficients (10). Under (9) it agrees
with the actual character trace at every integer Q>=4. It is an explicit
choice of completion: agreement at the integers alone would not exclude
arbitrary analytic additions involving, for example, sin(pi Q). Within
the specified rational configuration coefficients and fixed occupation
family, no such additions are made.

At Q1 it has the finite values and derivatives

```
beta_1(A)=beta_1(B)=-1,
partial_logQ beta_Q(A)|1=3/2,
partial_logQ beta_Q(B)|1=1/2.                                (12)
```

In particular the continued packet is **not** a channel whose value
vanishes at Q1 and is first activated by differentiation. It is a signed
rank1 insertion there, with no literal nontrivial S1 representation.
For the separately normalized partition ratio f=Z_[2]/Z, at fixed p,

```
f|1 = -P(A union B),
partial_logQ f|1
  = (3/2)P(A)+(1/2)P(B)
    -(1/2)Cov(1_(A union B),Sstar).                          (13)
```

The covariance includes the original measure derivative; it cannot be
replaced by differentiating beta alone. Thermal/root transmission to
original U must still use the complete normalized functional for each
geometry. Neither the sign of (12) nor the dimension polynomial alone
fixes that response. There is no inference of a four-leg continuum
amplitude, a Jordan collision, or a regular-endpoint activation.

## Scientific card

The new result is an exact stable-character formula with an explicit
finite-Q validity threshold, followed by the conditional two-class N25
occupation completion. The finite aliases explain why a fixed-Q
isotypic trace could not previously be differentiated without specifying
its continuation. The remaining geometric input is exactly (9), supplied
by the separate packing proof. No new configuration populations, source
scores, stochastic samples, or coupling grid were computed; (7)-(8)
are hand-derived finite character sums.
