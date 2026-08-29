# The all-order Arguin homology Q-jet

## Exact identity

Let

\[
H(Q,\tau)=\pi_Q(\mathrm{cross};\tau)-\pi_Q(\mathrm{trivial};\tau),
\qquad P(Q,\tau)=\pi_Q(\mathrm{trivial};\tau).
\]

The Arguin homology relation is

\[
H(Q,\tau)=(Q-1)P(Q,\tau).
\]

It fixes every ordinary derivative at `Q=1`:

\[
\partial_Q^n H|_1=n\,\partial_Q^{n-1}P|_1.
\]

For the critical-manifold score coordinate `t=log Q`, write
`D=partial_t=Q partial_Q`.  Directly differentiating
`H(e^t)=(e^t-1)P(e^t)` gives the cleaner all-order form

\[
\boxed{D^nH|_1=\sum_{k=0}^{n-1}\binom nk D^kP|_1.}
\]

Thus the logarithmic-Q jet is the truncated Pascal triangle.

## A necessary basis correction

At third order the two valid expressions are

```text
D^3 H|_1 = P + 6 partial_Q P + 3 partial_Q^2 P
          = P + 3 D P         + 3 D^2 P.
```

The coefficient `6` belongs to the ordinary `partial_Q` basis.  It must not be
used in front of `D P`.  The oracle derives the change of basis with signed and
ordinary Stirling numbers and recovers the Pascal coefficients through order
10.  It also checks every row on the polynomial basis `P=(Q-1)^m`.

## Consequence for Q-score tomography

The exact critical-manifold score modes estimate `D^k P` from a single
`Q=1` stream.  Therefore the topology-forced part of the order-`n` contrast
score is

\[
\sum_{k<n}\binom nk\,\mathbb E[I_{\rm trivial}H_k(T)],
\]

where `H_k(T)` is the normalized order-`k` score polynomial and `H_0=1`.

Define the typed residual

\[
R_n(\tau)=D^nH|_1-\sum_{k<n}\binom nk D^kP|_1.
\]

For the pure continuum homology identity, `R_n` vanishes at every modulus.
A nonzero residual is useful only after its semantics are typed: it can come
from an explicit `Q` derivative of the field/projector, a singular collision
residue, or a finite-lattice violation of the continuum identity.

The forced Pascal jet must be removed before a higher score order is called a
new rank-3 or higher logarithmic state.  Higher derivative order alone is not
new representation content.

## Evidence boundary

The derivative identities and coefficient conversion are exact consequences
of `H=(Q-1)P`.  They do not prove that a finite square-site observable has the
same generic-`Q` continuation, nor do they identify any residual with a
specific LCFT field.
