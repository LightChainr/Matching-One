# Projective essential-birth N65 smoke

**Status:** local variance/runtime smoke; the pilot mean did not select a model.

- Samples: 20,000 per shape in 20 aligned batches.
- Runtime: 0.159667 wall seconds on the recorded local run.
- Exact crosswalk gates: PASS.
- Direct rank-two paths across both shapes: 785.

## Non-micro projective support

| orientation | primitive lines | distinct chi4 | Var(chi4), path | Var(chi4), birth flux |
|---|---:|---:|---:|---:|
| first | 8 | 4 | 0.473711 | 0.530822 |
| second | 7 | 4 | 0.474137 | 0.517856 |

Both N65 shapes expose more than the quarter-turn-only tiny-control support; the projective mark therefore carries a genuinely varying chi4 value before any continuum model is fitted.

## Marked source/sink at p_ref

| coordinate | mean | batch SE |
|---|---:|---:|
| `first_j4_birth_re` | 2.4450587 | 0.0251 |
| `first_j4_exit_re` | 2.5078637 | 0.0249 |
| `first_j4_activity_re` | 4.9529223 | 0.0465 |
| `second_j4_birth_re` | -1.3827877 | 0.0134 |
| `second_j4_exit_re` | -1.3655043 | 0.00999 |
| `second_j4_activity_re` | -2.748292 | 0.0217 |
| `second_minus_first_j4_birth_re` | -3.8278464 | 0.0236 |
| `second_minus_first_j4_exit_re` | -3.873368 | 0.0223 |
| `second_minus_first_j4_activity_re` | -7.7012143 | 0.0414 |

`A4` is retained only to verify the Issue #156 fixed-p character crosswalk. The nonredundant production coordinates are the ingress and egress fluxes, with the old line retained at the second birth. `DIRECT_RANK2` has no line and enters the unmarked derivative with multiplicity two.
