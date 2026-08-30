# N65 charged-source reveal from the existing projective archive

No new permutations were generated. The frozen `q_A/q_D` sources were applied to the existing 20k aligned archive.

| orientation | channel | W | SE | z | birth | exit |
|---|---|---:|---:|---:|---:|---:|
| first | A | 0.336204267896 | 0.00205 | 164.127 | 3.42110897046 | 3.47899626478 |
| first | D | 0.0370407065872 | 0.00104 | 35.633 | 0.640526567485 | 0.627038212802 |
| second | A | 0.338344544086 | 0.0019 | 178.364 | 3.48303234796 | 3.44390381009 |
| second | D | 0.0369750143646 | 0.000645 | 57.366 | 0.628804925877 | 0.625908683408 |

Both A and D charged activations are resolved under the frozen marginal-z gate.
Their F3 one-points have the parameter-free phase `O_C(omega)=(omega-omega^2)W_C/2`; both lie on the +i ray.

The A-D cross response is statewise zero. H-to-A/D cross response reduces to the unweighted R-odd control; it is not forced to vanish batchwise.

- first unweighted `(A,D)` null: quadratic 4.51 / 2 df.
- second unweighted `(A,D)` null: quadratic 3.207 / 2 df.

The activation question and the orientation-modulation question are different:

- A `(W,birth,exit)` orientation score: `12.15 / 3 df`; `W_A` alone is only `0.900 sigma`.
- D `(W,birth,exit)` orientation score: `1.509 / 3 df`.
- Joint frozen six-vector: `15.53 / 6 df`.

Thus the charged sources themselves are precisely measured, while only the A current triplet carries a visible same-N orientation modulation in this engineering block.

Internal T-shear relabeling residuals are `4.16e-17` for the H/A/D vector and `4.16e-17` for the response matrix.
This is an exact representation check, not an independent identity/shear experiment: the two N65 orientations are different microscopic quotients.

The full aligned-batch covariance of both orientations, the frozen six-vectors and their contrasts is stored in `latest.json`.
