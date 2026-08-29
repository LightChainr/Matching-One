# P154 norm-4 fourth-generation pilot

The frozen N520/N680 acquisition and score are complete. Both targets contain
100,000,000 samples in 100 batches and use the pre-declared fourth-generation
geometries, seeds, replica-counter interval and runner commit from
`experiments/norm4_generation4_pilot_20260829.yaml`.

## Frozen score

| secondary eigenvalue | view | chi-square / df | p-value |
|---:|---|---:|---:|
| 1/2 (primary) | scalar U | 1.31420 / 2 | 0.51835 |
| 1/2 (primary) | thermal jet r2--r6 | 9.29810 / 10 | 0.50407 |
| 0 (fixed diagnostic) | scalar U | 1.23881 / 2 | 0.53827 |
| 0 (fixed diagnostic) | thermal jet r2--r6 | 9.06007 / 10 | 0.52641 |
| 1 (fixed diagnostic) | scalar U | 1.40714 / 2 | 0.49482 |
| 1 (fixed diagnostic) | thermal jet r2--r6 | 9.55161 / 10 | 0.48067 |

The two marginal scalar residuals for the frozen lambda=1/2 row are `+1.003`
and `-0.554` standard errors. All ten thermal-jet marginal residuals have
absolute value below `1.07` standard errors. Scalar U and the jet reuse the
same curves and are correlated views, not additive evidence.

## Mechanism reading

The preregistered Jordan-plus-one-even-mode recurrence at lambda=1/2 survives
its first fourth-generation target. The larger conclusion is narrower: this
100M pilot does not identify the secondary eigenvalue. Pure Jordan
(lambda=0), the frozen analytic even mode (lambda=1/2), and persistent
curvature (lambda=1) all receive nearly identical full-covariance scores.

Thus generation four removes the earlier visible tension but does not by
itself select a unique radial transfer law. A free post-reveal eigenvalue fit
would mostly profile a flat direction and is not promoted as new evidence.
The economical one-secondary-mode closure remains live; identifying its
eigenvalue requires either a more discriminating geometry/operator coordinate
or substantially sharper fourth-generation jets, not a relabeling of this
pilot.

## Numerical edge and provenance

The N680 bracket exposed endpoint underflow in the historical binomial-tail
recurrence: `(1-p)^680` vanishes in binary64 at the upper bracket although the
mass near the binomial mode is ordinary. The scorer now re-centers the
recurrence at the binomial mode only when that endpoint value is exactly zero.
The established calculation is bit-for-bit unchanged wherever the original
endpoint is representable, and a large-N deterministic-threshold regression
test covers the fallback.

The acquisition ran on Huawei DevEnv
`f550f3cb1f774374b6842aa648fda796` with two four-thread lanes. Accepted binary
SHA-256 is
`ee9010f524935099ba22f1820fddc05a79dd309d98784dfd5ba7da28129b6856`
from runner commit `bfab0330f5f56ca4d746b45d737f1607e3d229a0`. The run completed with status
`0,0`; stderr is empty. Remote and local SHA-256 digests agree for all six raw
artifacts.

Machine-readable results are in `analysis/score.json`; the exact reveal
command is recorded in `commands.txt`, and `checksums.sha256` covers the
submitted result bundle. Nineteen focused retrospective, production, thermal
jet, transfer and generation-four tests pass.
