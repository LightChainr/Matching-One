# Paired motif microcanonical certificate

All entries come from complete exact enumeration; no Monte Carlo samples are used.

| Gaussian pair | N | masks | K values | equal multiplicities | maximum fixed-K difference sum |
|---|---:|---:|---:|---:|---:|
| `gaussian-2-1 / gaussian-1-2` | 5 | 32 | 6 | `True` | `0` |
| `gaussian-3-2 / gaussian-2-3` | 13 | 8192 | 14 | `True` | `0` |

The N=5 and N=13 conjugate/swapped controls are configurationwise degenerate under the
deterministic quotient labelling (zero nontrivial masks); they certify the exhaustive counter and
fixed-K algebra, not useful control covariance.

Declared nontrivial same-N gates:

- `gaussian-8-1 / gaussian-7-4` (N=65): equal multiplicities; K=4 witness differences `{'nn_edge': 3, 'diagonal_pair': 2, 'face': 1, 'right_angle': 1}`.
- `gaussian-9-2 / gaussian-7-6` (N=85): equal multiplicities; K=4 witness differences `{'nn_edge': 3, 'diagonal_pair': 2, 'face': 1, 'right_angle': 1}`.
- `gaussian-11-3 / gaussian-9-7` (N=130): equal multiplicities; K=4 witness differences `{'nn_edge': -4, 'diagonal_pair': -2, 'face': -1, 'right_angle': -1}`.
- `gaussian-12-1 / gaussian-9-8` (N=145): equal multiplicities; K=4 witness differences `{'nn_edge': -3, 'diagonal_pair': -2, 'face': -1, 'right_angle': -1}`.
- `gaussian-13-1 / gaussian-11-7` (N=170): equal multiplicities; K=4 witness differences `{'nn_edge': 4, 'diagonal_pair': 2, 'face': 1, 'right_angle': 1}`.

Independent controls:

- `axis` (N=9): 512 masks, formula failures `0`, incremental failures `0`.
- `diamond` (N=8): 256 masks, formula failures `0`, incremental failures `0`.

All six unimodular-basis joint-histogram checks pass exactly.

## Interpretation boundary

No covariance, variance reduction, fitted control coefficient, production sample, matching/Euler identity, or >=2x promotion claim is included.
