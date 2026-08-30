# N65 flat-twist constraint-sector score

**Archive sufficiency: YES.** No new field or sample is required for F2/F3 flat-twist tomography.

For each fixed-p configuration the archive distinguishes rank zero, rank one with primitive `ell mod q`, and rank two. This determines every `T_alpha` exactly.

All aggregate/source-inversion gates pass with maximum residual `2.66e-15`.

## Same-modulus orientation contrasts

| projector | contrast | batch SE | |z| |
|---|---:|---:|---:|
| `F3_H4_axis_diag` | 0.0011029842 | 0.00151 | 0.732 |
| `raw_chi4_parallel` | 0.0021934716 | 0.00301 | 0.729 |
| `F2_H4_axis_diag` | 0.00090614408 | 0.00164 | 0.553 |

Within the same H4 sector, **`F3_H4_axis_diag`** is nominally sharpest at `|z|=0.732`, effectively tied with raw chi4 rather than a material variance improvement.

The balanced F3 projector is the minimal axes-versus-diagonals character on `P1(F3)`. It cancels the common rank-zero contribution automatically; unlike raw physical `chi4`, it uses unit finite-field orbit weights and retains only the twist constraint information.

## Additional twist-only sector

The sharper new contrast is **`F3_diagonal_odd` = 0.0020775613 +/- 0.00104** (`|z|=2`). The joint F3 axis-odd/diagonal-odd diagnostic is `5.903 / 2 df`.

This does not replace the H4 score: it is reflection-odd/projective and therefore a different finite-twist sector. Its value is that ordinary chi4 collapses this modular line information, while T_alpha retains it.

## Selected sector values

| value | mean | SE |
|---|---:|---:|
| `first_F2_S` | 2.3085302 | 0.00469 |
| `first_F3_S` | 4.2405839 | 0.0127 |
| `second_F2_S` | 2.3167377 | 0.00565 |
| `second_F3_S` | 4.2610875 | 0.0155 |
| `first_F3_H4_axis_diag` | 0.14958178 | 0.00126 |
| `second_F3_H4_axis_diag` | 0.15068476 | 0.0011 |
| `second_minus_first_F3_H4_axis_diag` | 0.0011029842 | 0.00151 |

The complete F2/F3 sector vector, all characters and their cross-orientation/cross-character covariance are stored in the JSON artifact.
