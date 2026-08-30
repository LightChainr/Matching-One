# Exact N=17 multi-orbit flux and N=13 cross-quotient comparison

Status: `exact_cross_quotient_multiorbit_flux`.

The N=17 Gaussian quotient 4+i again resolves rank-one traffic into axis and diagonal primitive-line orbits with exactly opposite chi4. Both orbits carry nonzero source and sink flux and reinforce dA4/dp at p_ref, so the N=13 mechanism is not quotient-specific.

## N=17 orbit gate

- period matrix: `[[4, -1], [1, 4]]`
- subset states: 131,072
- directed addition edges: 1,114,112
- direct 0->2 edges: 8,823
- orbit count: 2
- coefficientwise source/sink identity: `True`

| orbit | primitive lines | chi4 | rank-one states | birth | exit |
|---|---|---|---:|---:|---:|
| axis_orbit | [[0, 1], [1, 0]] | (161/289, 240/289) | 36516 | 150824 | 81600 |
| diagonal_orbit | [[1, -1], [1, 1]] | (-161/289, -240/289) | 2380 | 16218 | 9418 |

## Frozen p_ref signed-share comparison

| quotient | axis share | diagonal share |
|---|---:|---:|
| 3+2i (N=13) | 0.755739917417081006 | 0.244260082582918994 |
| 4+i (N=17) | 0.764844997919214790 | 0.235155002080785210 |
| N17-N13 | 0.00910508050213378321 | -0.00910508050213378321 |

The p_ref axis signed share moves from 75.574% at N=13 to 76.484% at N=17, a +0.911 percentage-point shift despite the different physical chi4 phase.

The character phase rotates between the two Gaussian generators, but the normalized source-minus-sink partition stays close and both line orbits reinforce the total derivative. That favors a transported projective-current mechanism over an N=13-only cancellation accident.

## Boundary

- This is an exact finite-volume source/sink localization, not a continuum-field identification.
- The cross-quotient share stability is a two-geometry mechanism clue, not an asymptotic limit.
- p_ref is inherited unchanged from the N=13 certificate; no parameter was fitted.
- No path enumeration, Monte Carlo sample, Huawei production, new PR, or merge is used.
