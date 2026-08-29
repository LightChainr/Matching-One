# Cross-microscopic amplitude invariants: an exact gauge boundary and one weight-12 null

Status: exact normalization bookkeeping plus a conditional, zero-amplitude modular prediction for Issue #118.

## 1. Why one bare H4 amplitude divided by one bare E6 amplitude cannot work

Let a square/matching realization supply a spin-4 response and a triangular
realization supply a spin-6 response at matched continuum modulus:

```text
A4^S(tau)=g4^S Z4^S mu_S^4 F4(tau),
A6^T(tau)=g6^T Z6^T mu_T^6 F6(tau).
```

The irrelevant couplings `g4^S,g6^T`, field normalizations `Z4^S,Z6^T`, and
constant coordinate metrics `mu_S,mu_T` are independent.  Under the two
coupling gauges, a monomial `(A4^S)^p(A6^T)^q` has charge `(p,q)`.  The charge
matrix is the identity, so its kernel is zero.  This is an exact no-go: no
nonconstant single-modulus monomial of the two bare amplitudes is universal.

Matching weights by cubing and squaring removes coordinate dimension but not
the two independent microscopic couplings.  Thus `(A4^S)^3/(A6^T)^2` is still
not a universal number.

## 2. Two shapes close both gauges

Use the same typed observable in each microscopic model at two matched moduli
`tau_j,tau_k`.  The weight-12 double ratio

```text
U46(j,k) = [A4^S(tau_j)/A4^S(tau_k)]^3
           /[A6^T(tau_j)/A6^T(tau_k)]^2                 (1)
```

cancels `g4^S,Z4^S,mu_S` separately from `g6^T,Z6^T,mu_T`.  It never equates
the two microscopic normalizations.  Its covariance-friendly form is the
polynomial null

```text
(A4_j^S)^3 (A6_k^T)^2 - (A4_k^S)^3 (A6_j^T)^2 = 0.     (2)
```

The powers are not fitted: `3*4=2*6=12`.

Under the conditional torus bridge

```text
F4(tau)=E4hat(tau),  F6(tau)=E6hat(tau),
```

the three degree-2 children of the hexagonal point obey

```text
E4hat_j = A*(1,zeta,zeta^2)_j,
E6hat_j = B*(1,1,1)_j.
```

Therefore, for every pair of children,

```text
U46(j,k)=1,                                           (3)
```

because `zeta^3=1`.  This is a genuinely cross-microscopic dimensionless
prediction: either lattice may have any nonzero constant coupling, field
normalization, or metric.  A failure of (3) rejects the joint pure-E4/pure-E6
bridge; it does not say which side failed.

The cheapest experiment is six common-frame response estimates: the complex
H4 coefficient of the square/matching model and the complex E6 center-score
coefficient of the triangular model on the same three registered child tori.
Fit no amplitudes.  Score the two independent complex polynomial nulls with
the full covariance.  Real-only H4 data are insufficient because they alias
the two conjugate C3 characters.

## 3. The log-pair mixed pairing supplies a second normalization closure

The Issue #234 same-flow archive already measures the symmetric Gram triple

```text
LL=<L L>,  LD=<L D>,  DD=<D D>.
```

Independent field rescalings act as

```text
(LL,LD,DD) -> (a^2 LL, a b LD, b^2 DD).
```

Consequently

```text
K_J=LL*DD/LD^2,       J_234=K_J-1                   (4)
```

are exactly normalization-free.  The mixed pairing `LD` is essential: using
only `LL` and `DD` leaves two independent field gauges and admits no
nonconstant continuous invariant.

There is one important boundary.  At finite cutoff, the allowed Jordan shear
`D -> D+alpha L` changes (4).  In the continuum null-bottom limit

```text
LL -> 0,  LD -> nonzero,
```

the prediction becomes

```text
K_J -> 0  (equivalently J_234 -> -1),               (5)
```

and is then shear-independent as well as rescaling-independent.  Thus the
existing triangular archive can score approach to (5), but its noisy
finite-cutoff `J` values are not universal constants.  A cross-microscopic
Jordan test needs only the same three coflow statistics in one independent
realization; no absolute energy-field normalization is required.

## Claim boundary

The gauge counting, cancellation in (1), C3 closure in (3), and Gram
transformation laws are exact.  The lattice identification with pure E4/E6
responses and the continuum Jordan limit are hypotheses.  No bare amplitude
and no finite-cutoff log-pair number is promoted to a universal constant.
