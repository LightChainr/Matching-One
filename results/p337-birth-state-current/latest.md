# P337 birth-state current prism

This is a zero-new-sample reconstruction from immutable joint birth archives. Each size is one dependency block; directions share one delete-one batch unit.

## Pooled-root state-current decomposition

| run | N | p_bar | E-dot | J01 | J12 | J02 | age1 completion | age beta | K_A activity |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| N85 | 85 | 0.592739048 | 0.00322521 +/- 0.00152 | -0.00250343 +/- 0.00136 | 0.000721779 +/- 0.00117 | 0.000157961 +/- 0.000378 | 0.000417532 +/- 0.000109 | -0.00235069 +/- 0.00569 | -0.0160511 +/- 0.00497 |
| N170 | 170 | 0.592732466 | 0.00236323 +/- 0.000353 | -0.00121874 +/- 0.000265 | 0.0011445 +/- 0.000278 | -1.31653e-05 +/- 6.04e-05 | 0.000222204 +/- 2.18e-05 | 0.00152379 +/- 0.000958 | -0.0111119 +/- 0.000976 |
| N340 | 340 | 0.592747082 | 0.000953495 +/- 0.000426 | -0.000287702 +/- 0.0003 | 0.000665792 +/- 0.000294 | -7.8025e-06 +/- 5.34e-05 | 6.73631e-05 +/- 1.76e-05 | 0.000171009 +/- 0.000639 | -0.00485732 +/- 0.00125 |
| N680 | 680 | 0.592746229 | 0.000371404 +/- 0.000207 | -0.000215755 +/- 0.000158 | 0.000155649 +/- 0.000151 | -9.35484e-06 +/- 1.52e-05 | 3.44134e-05 +/- 6.38e-06 | -3.90571e-06 +/- 0.000163 | -0.00216756 +/- 0.000557 |
| N130-control | 130 | 0.593023248 | -0.00761472 +/- 0.00733 | 9.04928e-05 +/- 0.00671 | -0.00752423 +/- 0.00496 | -0.00192944 +/- 0.00153 | -0.000101347 +/- 0.000524 | 0.019382 +/- 0.0244 | 0.0142445 +/- 0.0177 |

## Exact interpretation

For line paths, `J01` is first birth and `J12` is second completion; `J02` is the line-free direct `0->2` current. The scorer verifies

```text
F1' = J01 + J02
F2' = J12 + J02
E_top' = J12 - J01
A_top' = J01 + J12 + 2 J02.
```

Thus collision cancels from fixed-p E exactly. Its mass and current remain correlated coordinates, not a third projective-line component.

Across the four primary sizes the directional collision mass and `J02` current are unresolved. At N170 the resolved E response is assembled by a negative first-birth contrast and a positive second-completion contrast, so the two line currents reinforce rather than cancel. At N340 the resolved piece is concentrated in second completion; N85 and N680 do not individually resolve both pieces.

The completion-age first moment is directionally nonzero in all four primary blocks, while the conditional age-hazard coefficient is unresolved in this coarser archive. This distinguishes a robust age-weighted completion current from evidence that age itself causes the directional response. The richer N325/N425 P334 analysis separately rejects coarse-state age independence.

The age coefficient is a line/current-layer fixed-effect diagnostic. A nonzero value rejects sufficiency of the retained coarse state; it cannot distinguish intrinsic time memory from unrecorded current geometry.

## Dependency and boundary

N85/N170/N340/N680 are independent seed blocks. N130 is a separate cross-lineage control. These outputs do not share sample covariance with the ten-size E_top or Euler archives. Completion winding, complement line, transporter, ambiguity, microscopic state and path order are not scoreable.
