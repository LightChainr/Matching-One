# The prescribed Q1 closed trace transmits to original U through normalization

**Primary result.** The fixed continuation
`beta1=-1_A-1_B` has the strictly negative original pooled-root response

\[
\boxed{V_{\beta_1}=-0.001904836180602413.}
\]

The exact rational enclosure of `V_beta1/A25` is

```text
lower = -203814516535388901003038529243400709 / 10000000000000000000000000000000000000000
upper = -50953629133847225250759632310850177 / 2500000000000000000000000000000000000000
```

It excludes zero, with an outward-enclosure width of `1e-40`.
This rejects the single frozen normalization-transmission null at Q1.
It does not reinstate a regular invariant endpoint: the source is the
specified closed-seam trace continuation on the actual lattice.

## Fixed packets and accepted gates

Contract: `964ef2032effbe59f9158c158cf06a2c0844d7ee`.
Scorer committed before gate release: `bc93924c68d91a3299d1a1146bf258e7f7c2d997`.
Packing proof: `e901ba4a7fff8cd5eb644eb69d7a908945c5d69b`.
Generic-character proof: `58d91061267e90b72261d98226a103711251fc0e`.
Both were accepted by the coordinator before the single execution.

The types and continuations remain exactly those in the contract:

- A: `q=0,bad2=1,n_bad3=2`, with `beta_Q=(Q-3)/(2Q)`.
- B: `q=0,bad2=0,n_bad3=1`, with `beta_Q=(Q-3)/2`.
- All other allowed types have zero trace coefficient.

The complete reused histograms contained no unsupported mod6 pattern.
Their exact support counts are:

| Geometry | A configurations | B configurations | Zero-coefficient configurations |
|---|---:|---:|---:|
| axis `(5,0)` | 1,810 | 10 | 33,552,612 |
| tilted `(4,3)` | 850 | 55 | 33,553,527 |

These are configuration counts, not probabilities under the matching-root
measure. Both source packets have identically zero q-weighted and E-weighted
numerator polynomials, because all their support has ambient rank1.

## Three-term transmission and moving root

The old iid p-root, Dp and U/A were imported unchanged from
`results/p337-closed-source-n25/latest.json`. Only the coordinate change
`h=y=p/(1-p)`, `D_h=D_p/(1+h)^2` was applied. The matching root remains
`p0≈0.5926655393282267`, with no root search.

The existing normalization map includes separately normalized geometries,
the pooled root movement and the slope denominator. Its three named terms,
in full original-U units, are:

| Term | Primary `beta1` | Secondary fixed-gauge H |
|---|---:|---:|
| common thermal | `+9.501489779041238e-7` | `-1.3816936403095363e-5` |
| geometric thermal | `-0.0019059910369126577` | `+0.03827406448371185` |
| geometric value, including root movement | `+2.047073323405389e-7` | `+6.949599018239921e-7` |
| **sum** | **`-0.001904836180602413`** | **`+0.03826094250721058`** |

The corresponding h-root tangents are
`-8.95474165887118e-9` and `-3.040040780637895e-8`.
These are h=y derivatives, not p derivatives. The p tangents are obtained
by dividing by `(1+h0)^2`.

For the primary, the mean trace fractions in the two geometries are
`-5.350389311103746e-5` and `-4.9109172730798645e-5`. Their h derivatives
have opposite signs: `+3.737242708923633e-5` and
`-2.019067867477188e-5`. Thus the resolved transmission is carried mainly
by a difference of geometric thermal responses, not by a large root shift.
This is a decomposition of the same exact response, not three independent
mechanism tests.

## The secondary quantity has a different meaning

In the frozen reduced partition `y^K Q^(-(K+g)/2)`, the raw trace
Q-derivative packet is

\[
H=\frac{K+g+3}{2}\mathbf1_A+
  \frac{K+g+1}{2}\mathbf1_B.
\]

Passing this **one additive derivative packet** through the same linear
normalization map gives

\[
J_H=+0.03826094250721058,
\qquad J_H/A_{25}\approx0.0004093861497753176.
\]

Its rational interval is strictly positive and has outward width `9e-40`.
This is fixed-gauge raw-Q attribution. It is not a gauge-invariant share,
not `d_Q d_epsilon U`, and not an independent second evidence block.
Neither packet was changed after reveal.

## Scope and receipt

The finite-lattice trace-to-normalized-U interface is now explicitly
nonzero at Q1 for this named continuation. No continuum H4 assignment,
regular-endpoint activation, asymptotic size law or independent statistical
confirmation follows. These are deterministic uses of the same complete
N25 populations.

[score.json](../results/p337-q1-trace-continuation/score.json) retains all
rational bounds, coefficients, three terms, geometry fractions, input
hashes, gate SHAs and inherited-root provenance. The single managed
research-Python execution took about0.235 seconds. There was no new
enumeration, Monte Carlo, Q/seam scan, old Sstar/Bvac response calculation,
root search or test suite.
