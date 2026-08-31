# Where lifetime re-enters the complete topology profile

The cancellation of W=K2-K1 from the *unweighted* A_top integral does not
remove lifetime information from the whole profile. It specifies the first
thermal moment that can recover it, and the center fluctuation that must be
separated from it.

For a fixed permutation let
`A(p)=P(Bin(N,p)>=K1)+P(Bin(N,p)>=K2)-1`, and define
`J_q=integral_0^1 p^q A(p) dp` for integer q>=0. With the rising factorial
`(x)_[m]=x(x+1)...(x+m-1)`, an exact identity is

\[
 \boxed{J_q=\frac{1}{q+1}\left[
 1-\frac{(K_1)_{[q+1]}+(K_2)_{[q+1]}}
 {(N+1)_{[q+1]}}\right].}
\]

To derive it, regard the binomial tail as the CDF of the K-th order statistic
U_K of N independent uniforms. Then
`integral p^q 1{U_K<=p} dp=(1-U_K^(q+1))/(q+1)`.
The beta integral gives
`E[U_K^m]=(K)_[m]/(N+1)_[m]`. Apply this to the two thresholds and subtract
the integral of one. This is an identity of the canonical kernels, not an
additional random simulation or smoothing assumption.

## The first two levels separate clock center and width squared

Write C=(K1+K2)/2 and W=K2-K1. The first two identities are

\[
 J_0=1-\frac{2C}{N+1},\qquad
 J_1=\frac12-\frac{C^2+C+W^2/4}{(N+1)(N+2)}.
\]

Thus the first p-weighted moment contains *two* new terms beyond the mean
center: its second moment and the squared lifetime. An apparent lifetime
signal in J1 cannot be identified from J1 alone. The same-path clock-center
coordinate supplies an exact subtraction:

\[
 \boxed{J_1-\frac12+\frac{C^2+C}{(N+1)(N+2)}
       =-\frac{W^2}{4(N+1)(N+2)}.}
\]

For the existing orientation contrast, the corresponding relation is

\[
 \Delta J_1=-\frac{\Delta(C^2+C)+\Delta(W^2)/4}{(N+1)(N+2)}.
\]

The same fixed cos(4theta) normalization may be applied to every term.
Because C and W come from the same two births, their covariance is part of
the signal; this is not a prescription to combine independent error bars.
All necessary baseline quantities are already present in the completed
9c495ab1 counter archive. No longer-horizon simulation is needed to obtain
this particular lifetime-sensitive observable.

## The finite polynomial selection rule

Every numerator `(C-W/2)_[q+1]+(C+W/2)_[q+1]` is even in W. Consequently
the complete A_top thermal-moment hierarchy contains only even lifetime
powers: none in J0, W^2 first in J1, then higher even powers as the degree
increases. The plateau profile F1-F2 instead uses their difference and is
odd in W; its unweighted integral is W/(N+1).

This is a birth-exchange algebraic rule, not a claim that these powers are
independent fields, asymptotic correction exponents or new random sources.
It also shows why completing the global observable need not end the lifetime
route: it moves that route to a specified thermal moment, with an explicit
center-moment subtraction.

Scientific card: exact finite-N canonical identity, same two-birth source
and same orientation pair. New discriminator is J1 together with C^2+C,
not another fit family or a scan over moment orders. The present note derives
the relation only; it does not report a new numerical lifetime anisotropy.
