# A positive microscopic source closes after two saturated decimations

## Named model and exact prediction

On an ordinary square-site torus define three **configuration counts**:

```text
C = occupied NN components + vacant matching components,
F = fully occupied elementary faces,
Bv = NN edges whose two endpoints are vacant.
```

Saturate one checkerboard sublattice, contract the remaining coset by1+i,
and complement it to an ordinary square child. The proposed exact dictionary
is

```text
C_parent = C_child + F_child,
F_parent = Bv_child,
Bv_parent = 0.
```

The first line was proved in the existing checkerboard work. The second
comes from the face-to-edge bijection: the two unsaturated corners of a
parent face become adjacent child sites, and their occupations are reversed
by complement. The last line holds because every parent NN edge touches
the saturated sublattice. Use honest nonalias quotients and the same lifted
winding map; full companion proofs accompany this calculation.

It follows that the **positive count source**

`S_star = C + F + Bv`

is preserved configuration by configuration. It is not a regression-weighted
combination and no coefficient is inferred from the F4 or lag1 data. This
supplies a concrete finite microscopic law

`P_(p,t)(omega) proportional to P_p(omega) exp(t S_star(omega))`.

For its saturated endpoint, the normalized entire child law is the same
law at `(1-p,t)`. Consequently on corresponding simple-root branches,

```text
p0_parent,end(t) = 1 - p0_child(t),
U_parent,end(t) = 2^(13/8) U_child(t).
```

Here U is unchanged: the same pooled homology-balance root, area factor,
DeltaCos4 projector, and thermal-derivative ratio. This identity holds at
finite source strength wherever that root branch is defined, not only for
the source derivative at t=0. It does not assert a homogeneous critical
renormalization fixed point or a continuum field identity.

## The source closure is finite, positive, and coefficient-fixed

For couplings g=(c,f,b) of (C,F,Bv), the exact child coupling is

```text
Tg = (c,c,f),
T^2g = (c,c,c),
T^3g = T^2g.
```

Thus a bare cluster source reaches S_star after exactly two such boundary
decimations. Normalizing the C coefficient to1 leaves the fixed coefficient
vector(1,1,1). The initial F and Bv couplings disappear after two steps:
their different twice-decimated normalized configuration laws agree when
their C coupling agrees. This is more than equality of one scalar response.

In particular the earlier two-source family C+F is not closed. The missing
term has the specified vacant-bond shape; it is not another face descriptor
selected from a residual. There is no further local count to guess in this
three-step chain. The nilpotent transient in this finite count transformation
is not a logarithmic-CFT Jordan block.

The ordinary occupied-bond representation is equivalent:

`Bv=2N-4K+Bocc`,

so `S_star=C+F+Bocc-4K+2N`. The -4K term is a common chemical-potential
shift, not an extra unidentified physical source. At the same source strength,
removing that term changes logit(p) by a common amount; root/slope-normalized
U is unchanged by that reparameterization. The manifestly positive count
representation avoids having to hide such a shift in the endpoint map.

## Nested saturation is a concrete spatial environment

Repeated steps require the period lattice to remain checkerboard-compatible
at every intermediate quotient. To descend k times to a child of area M,
take a parent of area2^k M with periods multiplied by(1+i)^k. Short-period
aliases must remain outside the chosen elementary-stencil implementation.

The physical parent pattern is not obtained by occupying the same color
again without complement. Successive layers fix alternating **occupied and
vacant** original sites. Their fractions are

```text
occupied_fixed = (2/3) * (1 - 4^(-ceil(k/2))),
vacant_fixed   = (1/3) * (1 - 4^(-floor(k/2))),
free_fraction = 2^(-k).
```

The residual free occupation is complemented k times. This gives exact
relations, on the declared nested endpoint only,

```text
U_(2^k M),nested(t) = 2^(13k/8) U_M(t),
rho_(2^k M),nested(t) = 2/3 + (-1/2)^k [rho_M(t)-2/3],
S_star,parent(omega) = S_star,child(omega_child).
```

rho is the actual occupation density in the source-tilted law; for t=0 its
child value is the Bernoulli p. The2/3 density is a feature of this explicitly
inhomogeneous forced hierarchy, **not** a square-site critical threshold.
There is no assertion of an interpolating U curve between the ordinary
homogeneous parent and the nested endpoint.

## Direct finite computation specified before its result

The next calculation is restricted to the coefficient-fixed S_star on the
existing honest N25 pair `(5,0),(4,3)`, at its own pooled root and original U.
It reports V_Sstar and, separately, V_Bv required by the second line of the
source dictionary. There is no source search, new lag, Monte Carlo top-up,
or re-scoring of the completed N65/N85/N130/N170 F4 block.

The per-K counts saved by the preceding exact enumeration did not retain
Bv. The existing rollback enumerator can directly collect the missing
statistics from `Bv=2N-4K+Bocc`; no new connectivity engine is needed.
The entire fixed finite population is summed, so this is deterministic
coefficient calculation, not another independent statistical evidence vote.
A nonzero rational enclosure of V_Sstar/A would eliminate a pure common
thermal-clock description of this named model at this finite observer.
It would not identify the continuum H4 field or an all-size exponent.
