# Exact axis L=5 matching polynomial and frozen zero-map score

The C++ kernel exhausted all `2^25 = 33,554,432` configurations. Its
L=1..4 coefficients were first regressed exactly against the Python oracle.
PR #78's three predictions were read unchanged from its committed JSON.

## Frozen score

| metric | frozen prediction | observed | observed - predicted | relative error |
|:---|---:|---:|---:|---:|
| physical_root_0_1 | 0.592597902243523 | 0.5919882565183338 | -0.0006096457252 | -0.0010287679 |
| imaginary_rms | 0.5332328861951509 | 0.3401748476358213 | -0.1930580386 | -0.36205201 |
| nonreal_fraction | 0.7885714285714286 | 0.8 | 0.01142857143 | 0.014492754 |

The cheap two-point `a+b/N` rule is surprisingly close for the physical root
and gets the discrete nonreal fraction nearly right, but badly overpredicts the
imaginary RMS. The latter is a clean falsification of that frozen cloud-scale rule.
No alternative model was fit after seeing L=5.

## Exact and numerical audit

- Exact power-basis coefficients (ascending): `[-1, 0, 0, 0, 0, 10, 0, 100, -200, 300, 240, -3050, 2850, 1950, 6850, -26220, 20450, 3450, -2100, -18650, 22970, -10800, 1650, 300, -100, 2]`.
- Degree / real / nonreal roots: `25` / `5` / `20`.
- Physical root: `0.5919882565183338446109686802119288790479`.
- Maximum normalized polynomial residual: `5.33153e-120`.
- Maximum 100/130-digit root shift: `4.5593e-90`.
- Maximum conjugate-pair error: `0.0`.
- Maximum exact-partner `z -> 1-z` pairing error: `5.0e-120`.

This is a finite-size exact result only; the root cloud is not assigned a CFT or
Lee-Yang interpretation.
