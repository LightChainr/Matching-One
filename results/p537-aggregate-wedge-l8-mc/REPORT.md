# Square-L8 prospective `G4`: the frozen three-way test is unresolved, but a `35/8` law emerges

## Prospective result

The square-L8 production was frozen at commit `30cd7c00` before any output was
inspected: 100 million Bernoulli backgrounds, 100 batches, 16 disjoint shards,
proposal `p*=0.5926`, seed `20260901845`, and the unchanged radius-one landing
source and root-conditioned `G4` coordinate.

At the reweighted matching root,

```text
p8*          = 0.5926379848582818
G4(8)        = 0.0006802040610058551
SE[G4(8)]    = 0.0000126970016013406
95% interval = [0.000655317937867228, 0.000705090184144483]
S            = -0.0133766503649138
```

The full interval crosses the frozen geometric boundary
`0.000676620613854308` between the direct-L6 `L^-9/2` and `L^-4`
continuations. The registered decision is therefore
`UNRESOLVED_MODEL_BOUNDARY`; no samples were added and no coordinate was
changed. Conditional on the L6 anchor, the L8 point lies `-3.69` L8 standard
errors from the `L^-4` continuation and `+3.98` from the `L^-9/2`
continuation. The old three-way model set is too coarse.

## New phenomenon found after opening the L8 result

The exact L4 and L5 values already contain a much sharper exponent:

```text
-log(G4(5)/G4(4))/log(5/4) = 4.374872226062653
35/8                             = 4.375
```

Equivalently,

| L | source | `L^(35/8) G4` |
|---:|---|---:|
| 4 | exact population | 6.0311692933 |
| 5 | exact population | 6.0313412560 |
| 6 | 100M production | 5.8310474546 |
| 8 | prospective 100M production | 6.0765617219 |

The two exact amplitudes differ by only `2.85e-5` relatively. Their mean,
fixed without using L6 or L8, predicts

```text
G4(8) = 0.000675132503964198.
```

The prospective L8 result differs from that value by only `+0.40` of its own
standard error. The global effective powers are also close:
`p(4,8)=4.36418` and `p(5,8)=4.35911`.

This is a numerical discovery, not a fitted continuum identity. In particular,
the L6 amplitude is `-3.13` L6 standard errors below the exact-L4/L5 `35/8`
continuation. Thus the current pattern is a `35/8` backbone with a non-monotone
finite-size or commensurability correction, not a clean four-point power-law
fit. The scientific gain is that L8 rejects the previously suggested monotone
`9/2` continuation and exposes the exact L4/L5 exponent that it had obscured.

## Mechanism candidate

`35/8 = 4 + 3/8`. The `4` is the established canonical pair normalization.
The extra `3/8` is numerically half the square-lattice thermal eigenvalue
`y_t=3/4` in linear-size units. A possible mechanism is therefore a leading
`L^-4` collar term whose root-conditioned coefficient cancels, leaving a
square-root thermal-coordinate mismatch. This square-root step is not yet
derived; it is the precise theoretical gap exposed by the production.

## Boundary

This remains the reduced landing interface, not the full pooled-root
original-`U` response. No new claim is made that it transmits to the global
observable. The next production should therefore change scientific level—to
the full one-defect diagonal-edge gate—rather than add another reduced-`G4`
size.
