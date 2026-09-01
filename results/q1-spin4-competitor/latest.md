# Q=1 spin-4 competitor preflight

This is exact generic-loop spectrum arithmetic, not a Potts multiplicity or lattice-overlap claim.

## Exact field data

- beta^2 = 2/3; c = 0
- legs = 4
- (Delta, DeltaBar) = (1/8, 33/8)
- x = 17/4; spin = -4
- x(Q4 epsilon) - x(V_(2,2)) = 1

## Continuum two-field dilation oracle

| area multiplier Q | length dilation | x=17/4 factor | x=21/4 factor | relative Q4/four-leg factor |
|---:|---:|---:|---:|---:|
| 2 | 1.4142136 | 0.22925101 | 0.16210494 | 0.70710678 |
| 4 | 2 | 0.052556026 | 0.026278013 | 0.5 |
| 5 | 2.236068 | 0.032710617 | 0.014628633 | 0.4472136 |
| 10 | 3.1622777 | 0.0074989421 | 0.0023713737 | 0.31622777 |

## Unresolved gates

- potts_q1_multiplicity: unresolved_requires_representation_limit
- global_matching_matrix_element: unresolved_requires_selection_rule
- local_pivotal_overlap: unresolved_requires_character_weighted_control
- normalized_shell_transfer: unresolved_requires_lattice_to_radial_normalization

The relative factor follows only from the exact dimension gap. It must not be applied to the normalized local shell observable until that observable's lattice-to-radial normalization is derived.
