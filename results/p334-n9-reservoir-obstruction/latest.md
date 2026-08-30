# P334 N9 reservoir obstruction and two-mark repair

## Exact obstruction

All six deterministic matching/layer-4/Y=0 rows have `(D,M,Y,F)=(216,432,0,72)`. Rows `1` and `3` are the Smith-(3,3) representatives. Because `Y=0`, the `YN` target channel is empty. The existing one-carrier plus one-output-mark reservoir has coarse demand `6912`, reaches `4752` ordinary `MM` targets, and leaves deficiency `2160=5/16`.

On every candidate row the residual minimum cut contains all 768 coarse classes (192 from each source replica). This is an all-site image-capacity obstruction, not a small exceptional family.

## Minimal legal repair

keep both lower bases fixed and release two of the four ordered output-mark slots to arbitrary quotient sites; reapply the fixed-line topology gate.

The output remains an ordinary untagged `M x M` pair: no released-slot, phase, source, or provenance label is retained, and every output face is reclassified on the frozen projective line.

The repair reaches `20736` normalized targets, exactly `M^2/N = 432^2/9 = 20736`. This adds `15984` target orbits (`143856` raw tokens) beyond the old image. Every coarse class has degree `216`, and exact integral max flow is `6912/6912`.

This is minimal along the output-mark axis: the failed existing builder already strictly contains all fixed-base zero/one-mark releases, while two released slots with no base mutation saturate.

## Same-descriptor isomorphism

The strict descriptor (including `3x3` HNF/Smith type but omitting the line label) selects exactly rows 1 and 3. The explicit site permutation induced by `(x,y)->(y,x)` is an involution, swaps lines `(0,1)` and `(1,0)`, maps every D/M/Y/F face family bijectively, and conjugates the translation group to itself. Hence both old and repaired compatibility graphs are exactly isomorphic.

## Complete N9 Y=0 candidate class

| row | HNF | Smith | line | old flow | deficit | two-mark flow | MM orbits |
|---:|---|---|---|---:|---:|---:|---:|
| 1 | `[[3, 0], [0, 3]]` | `[3, 3]` | `[0, 1]` | 4752/6912 | 2160 | 6912/6912 | 20736 |
| 3 | `[[3, 0], [0, 3]]` | `[3, 3]` | `[1, 0]` | 4752/6912 | 2160 | 6912/6912 | 20736 |
| 6 | `[[3, 1], [0, 3]]` | `[1, 9]` | `[1, 0]` | 4752/6912 | 2160 | 6912/6912 | 20736 |
| 9 | `[[3, 2], [0, 3]]` | `[1, 9]` | `[1, 0]` | 4752/6912 | 2160 | 6912/6912 | 20736 |
| 15 | `[[9, 3], [0, 1]]` | `[1, 9]` | `[1, -3]` | 4752/6912 | 2160 | 6912/6912 | 20736 |
| 24 | `[[9, 6], [0, 1]]` | `[1, 9]` | `[2, -3]` | 4752/6912 | 2160 | 6912/6912 | 20736 |

There are two translation-equivariant classes. Rows 1/3 have group `Z3 x Z3` and are D4-isomorphic; rows 6/9/15/24 have group `Z9` and explicit D4 maps from row 6. The classes cannot be translation-equivariantly isomorphic because their group exponents are 3 and 9. Nevertheless all six have the same exact old and repaired flow signature.

## Scientific card

- **Question:** Why does the corrected combined reservoir first fail at N9, and what is the smallest local repair?
- **Answer:** Y=0 removes synergy and the one-mark MM image has a 5/16 all-site Hall deficit.
- **Repair:** A second output-mark release, with bases fixed and no decorated capacity, reaches every MM orbit and saturates.
- **New capacity:** 15984 previously unreachable MM orbits (143856 raw tokens), all already present in M^2.
- **Boundary:** Exact for all six N9 matching/layer4/Y=0 rows; not an arbitrary-HNF theorem.
