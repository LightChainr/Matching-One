# Third-geometry exact falsifier

Status: `exact_third_geometry_falsifier`.

## Flux-blind selection

First HNF in increasing N/lexicographic order with four distinct face corners, no physical quarter-turn lattice symmetry, exactly two projective stabilizer orbits, and distinct constant non-real chi4 on both orbits. Source/sink counts are forbidden during selection.

Selected HNF `[[7, 2], [0, 1]]` at N=7 after 29 geometry-only rows. It has no physical quarter-turn symmetry and is not Gaussian-similar.

Primitive-line orbits: `[[[0, 1]], [[1, -3]]]`.
Characters: `[{'real': '-7/25', 'imag': '24/25'}, {'real': '7/25', 'imag': '24/25'}]`.

## Exact subset-boundary census

- states: 128
- directed edges: 448
- direct 0->2 edges: 14
- coefficientwise continuity: `True`

| orbit | line | chi4 | states | birth | exit | net zero |
|---|---|---|---:|---:|---:|---:|
| orbit_0 | [[0, 1]] | (-7/25, 24/25) | 35 | 84 | 49 | 0.592783237894885157837638571036 |
| orbit_1 | [[1, -3]] | (7/25, 24/25) | 7 | 28 | 21 | 0.571428571428571428571428571429 |

## Frozen score

The zero separation is `0.0213546664663137292662099996079` versus the frozen parent-envelope `0.0409499231130266717212867090127`: close-pair gate **passes**.

The exact character Gram is `527/625` (positive). Therefore the two contributions cancel between the zeros and reinforce outside:

- `(0, 0.571428571428571428571428571429)`: positive/positive; **reinforce**.
- `(0.571428571428571428571428571429, 0.592783237894885157837638571036)`: positive/negative; **cancel**.
- `(0.592783237894885157837638571036, 1)`: negative/negative; **reinforce**.

At `p_ref`, orbit nets are `0.000248142518380236771432082690140` and `-0.0360801677033176642178956766900`; total phase is `-106.4732826` degrees.
The net slopes at their own zeros are `-6.67298485293693831983311735865` and `-1.67930029154518950437317784257`; the first orbit is only `0.000109418866625653742190674952498` of its source-plus-sink activity at `p_ref`.

Frozen verdict: **paired timing survives, but between-zero reinforcement is falsified**.

## Mechanism update

The close pair of source/sink balance times survives the asymmetric HNF quotient, but the Gaussian claim that reinforcement occurs only between them is false. The timing pair and the spin-4 alignment are separate layers: the latter flips with the exact chi4 Gram sign.

Exact two-orbit alignment theorem:

```text
Re[(chi1 J1) conjugate(chi2 J2)] = Re[chi1 conjugate(chi2)] J1 J2 = Gram(chi1,chi2) J1 J2
```

Paired zeros control only `sign(J1 J2)`. Gaussian opposite characters have negative Gram; this HNF has positive Gram, so the same timing window acquires the opposite alignment topology.

Next falsifier: For a fresh two-orbit quotient, the sign of Re(chi_a conjugate(chi_b)) must determine the alignment topology: positive Gram gives cancellation between the simple net zeros and reinforcement outside; negative Gram gives the reverse. Failure rejects the two-scalar-current reduction.

## Boundary

- The HNF candidate was selected from geometry and line-orbit support before boundary flux counts were read.
- This exact N=7 quotient is an asymmetric finite-volume falsifier, not an asymptotic geometry.
- The Gram-controlled alignment rule is a new mechanism classification, not a continuum-field identity.
- No N13/N17 recomputation, Monte Carlo sample, Huawei production, PR, or merge is used.
