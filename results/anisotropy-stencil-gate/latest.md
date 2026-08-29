# Exact spin-4 stencil feasibility gate

## Decision

Reject a blind zero search inside the horizontal/vertical-only exactly-critical square-bond family.
For every nonzero nonnegative axis weighting, `A4` is strictly positive.

## Exact calculations

| stencil/proxy | exact `A4` | zero? |
|---|---:|---|
| axis shell, arbitrary nonnegative weights | `2(w_h+w_v)` | only the empty stencil |
| critical axis family, probability weights | `2` | no |
| critical axis family, Bernoulli-variance weights | `4t(1-t)` | only degenerate endpoints |
| C4 axis + integer-diagonal shells | `4w_axis-16w_diagonal` | `w_axis/w_diagonal=4` |

The unit axis orbit has spin-4 phase `+1`; the integer diagonal has moment `-4` because
its squared length is two and its spin-4 phase is `-1`.

## Geometric gate

With reflection-cancelled imaginary part and nonnegative weights, cancelling a retained positive axis contribution needs at least one orbit with cos(4*theta)<0.
Thus an improved-action candidate with nonnegative microscopic weights must contain an oblique
orbit in a negative spin-4 phase sector. The exact axis/diagonal cancellation certificate is
`w_axis:w_diagonal = 4:1` for the declared integer-step normalization.

## Evidence boundary

This is an exact statement about the declared microscopic fourth-moment proxy, not a derivation
of the renormalized spin-4 coupling. The calculation neither constructs an exactly-critical
mixed axis/diagonal model nor predicts a vanishing measured H4 amplitude. The next admissible
step is to find an exactly-critical star-triangle/isoradial family with a negative-phase orbit
and derive the physical edge weights before simulation.

## Reproduction

```bash
python scripts/anisotropy_stencil_gate.py --format json
python scripts/anisotropy_stencil_gate.py --format markdown
python -m unittest tests.test_anisotropy_stencil_gate
```

Sources: Grimmett--Manolescu `arXiv:1105.5535` and `arXiv:1204.0505`.
