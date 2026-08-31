# A positive microscopic source closes after two saturated decimations

## Named model and exact prediction

On an ordinary square-site torus define three **configuration counts**:

```text
C = occupied NN components + vacant matching components,
F = fully occupied elementary faces,
Bv = NN edges whose two endpoints are vacant.
```

Saturate one checkerboard sublattice, contract the remaining coset by1+i,
and complement it to an ordinary square child. The exact dictionary
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

The [finite quotient proof](checkerboard-positive-source-closure.md)
establishes these identities and their winding/partition-function scope.
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

## The birth-channel map is fixed too

Let F1(p,t)=P_(p,t)(r>=1) and F2(p,t)=P_(p,t)(r=2). At t=0 these are the
usual first/second birth CDFs; at nonzero t they denote static rank-event
probabilities, without asserting a monotone process representation. The
rank complement gives

```text
F1_parent(p,t) = 1 - F2_child(1-p,t),
F2_parent(p,t) = 1 - F1_child(1-p,t).
```

Write the two readout contributions on their common root as
`U01=-A P4(F1_p)/D` and `U12=+A P4(F2_p)/D`. The thermal chain rule and
DeltaCos4 rotation then give

```text
U01_parent,end(t) = 2^(13/8) U12_child(t),
U12_parent,end(t) = 2^(13/8) U01_child(t).
```

This same swap holds for their complete t derivatives including the
root/slope motion. Two steps return each readout channel to itself.
These are readout contributions to a common full-source response, not
separately manipulated event-source interventions. Thus the named model
has a fixed microscopic source -> rank-event -> original-U map, rather
than an inferred assignment of entry versus completion after seeing data.

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

The [cluster-gas action](decimation-closed-cluster-gas-action.md) gives a
more physical form using the occupied graph cycle dimension beta1:

`S_star = 2 beta1 - 3K - q + 2N`.

Modulo the common thermal K term and normalization, this is a cycle-weight
source with the **fixed** rank correction -q. Contractible cycles carry
weight2 at fixed K/rank; dropping -q would define a different model. This
action identity gives meaning to the three counts, but does not separate
their individual contributions to the observed global-U response.

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

## The fixed source reaches the original global U

The [pre-calculation contract](../analysis/p337_closed_source_n25_contract.json)
selected only S_star and the dictionary-required Bv companion. Two complete
2^25 configuration sums on `(5,0),(4,3)`, followed by one rational-interval
evaluation, give

```text
p0                  = 0.5926655393282267,
U25                 = 0.8804661569633677,
V25(S_star)         = +0.12616536341416915,
V25(Bv)             = +0.33291070842057197.
```

Here V=dU/dt at t=0, including per-geometry covariance centering, the
pooled-root displacement and denominator/slope motion. The exact rational
enclosure of V25(S_star)/[25^(13/8)/2] is strictly positive: **the closed
source is not a pure common thermal-clock alias at this finite observer**.
See the [result and all rational bounds](../results/p337-closed-source-n25/REPORT.md).
The entire compile/enumerate/score run took3.11seconds locally; no cloud job,
random sample, source fit or old C/F4 response re-score was used.

The preceding saved per-K profiles lacked occupied-edge sufficient
statistics, so this one pass collected S_star/Bv using the inherited
connectivity engine. These are exact finite-population calculations, not
independent statistical votes for a continuum mechanism.

The identities above now imply the following derivatives without any new
parent simulation:

| Specified parent endpoint | Original-U source response |
|---|---:|
| N50, one saturated level, closed S_star source | +0.389147178497717 |
| N50, one saturated level, F4-only source | +1.026836996840865 |
| N100, two nested levels, bare C source | +1.200293982712272 |

The first and third are 2^(13/8)V25(S_star) and4^(13/8)V25(S_star);
the second is2^(13/8)V25(Bv). This turns the microscopic source ->
rank-event -> U dictionary into nonzero, parameter-free finite endpoint
predictions. These forced-endpoint responses are not homogeneous N50/N100
measurements, nor an extrapolation of the stopped larger-N F4 block.

## Scientific card and next discriminant

- **Mechanism change:** exact three-count closure plus a nonzero original-U
  response removes common-thermal invisibility of this coefficient-fixed
  finite model. Bare cluster coupling has no additional count to generate
  after its second saturated reduction.
- **Lifecycle:** algebraically specified; contract55988126 before computation;
  code d2cec6a7; exact result d3d60b44; integrated branch delivery, not main.
- **Observer/sector/source/geometry:** original pooled-root normalized U,
  ordinary homology q/E with normalized cos4 projector, bulk exp(t S_star),
  honest N25 axis Z5xZ5 and tilted Z25 quotients. The Smith classes differ.
- **Dependency:** exact integer profiles in `p337-closed-source-n25`; the
  inherited geometry/scoring code is declared in run provenance. The two
  source responses and all transported endpoints belong to this one finite
  calculation; none is an independent production confirmation.
- **Boundary:** no continuum field identity, asymptotic exponent, interior
  saturation law or universal homogeneous RG fixed point follows. The F4
  random-block `INCONCLUSIVE_STOP_FIXED_BLOCK_WITHOUT_TOP_UP` is unchanged.
- **Subsequent discriminant completed:** the
  [single-defect calculation](checkerboard-single-defect-global-u-result.md)
  retains S_star and the same U, and rejects a source-independent geometric
  gain immediately off saturation: R=U U_st-U_s U_t is strictly positive.
  The exact saturated identity still holds. An interior theory must now
  carry the measured source-dependent defect term; neither a fitted extra
  source nor another endpoint replay supplies that theory.
