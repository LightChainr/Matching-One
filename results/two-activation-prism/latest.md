# P205 K1/K2 quotient-prism reanalysis

## Outcome

The frozen quotient prism selects `H4/H4`: both the first ambient-homology activation (`K1=K_minus`) and the second (`K2=K_plus`) are compatible with the same H4 character line and reinforce at the fixed P205 probability.

This is a retrospective reuse of branch-only archives. It generated no Monte Carlo samples, kept `p_ref=0.59274605079` and `N^-13/8` fixed, and fitted only one amplitude per activation component.

| component model | fitted amplitude | SE |
|:---|---:|---:|
| K1 H4 | 0.5085811186 | 0.0079461494 |
| K2 H4 | 0.2952662594 | 0.0087394862 |

K2 supplies 36.732% of the signed fitted H4 amplitude. The joint score is chi-square=2.585155 on 4 df (p=0.629455). The runner-up is `H4/H12` at delta chi-square=8.918835.

## Activation-resolved contrasts

| N | Delta F1 | Delta F2 | Delta M | K1/K2 correlation | original closure |
|---:|---:|---:|---:|---:|---:|
| 25 | +5.017698803e-03 | +2.941457190e-03 | +7.959155993e-03 | +0.229 | -7.29e-17 |
| 50 | +1.589802045e-03 | +8.145610342e-04 | +2.404363080e-03 | +0.377 | -8.41e-17 |
| 125 | +2.194732386e-04 | +1.798596265e-05 | +2.374592012e-04 | +0.563 | +1.19e-16 |

Every size reinforces: `DeltaF1` and `DeltaF2` have the same sign. The complete 6x6 covariance in JSON retains each within-size K1/K2 correlation; cross-N blocks are zero only because the archived counter streams are distinct.

## Frozen nine-pair score

| K1 line | K2 line | chi-square / 4 df | p-value | delta chi-square |
|:---|:---|---:|---:|---:|
| H4 | H4 | 2.585155 | 0.629455 | 0.000000 |
| H4 | H8 | 278.037854 | 5.90247e-59 | 275.452699 |
| H4 | H12 | 11.503990 | 0.0214473 | 8.918835 |
| H8 | H4 | 1091.570540 | 5.08494e-235 | 1088.985385 |
| H8 | H8 | 960.920452 | 1.05035e-206 | 958.335298 |
| H8 | H12 | 1110.549141 | 3.91388e-239 | 1107.963986 |
| H12 | H4 | 68.315288 | 5.14687e-14 | 65.730133 |
| H12 | H8 | 354.372032 | 1.9951e-75 | 351.786878 |
| H12 | H12 | 41.970730 | 1.69163e-08 | 39.385575 |

## Boundary

The total P205 prism was frozen before reveal, but assigning separate character lines to K1 and K2 is a retrospective analysis. The source archives remain `branch_only`, all three sizes are deliberately small, and the result identifies neither a continuum operator nor an asymptotic theorem. It does show that the ordinary matching H4 signal on this exact quotient code is not solely a first-activation effect: the second activation supplies a resolved reinforcing share.

The JSON records the immutable source commit and input hashes, verifies `DeltaF1+DeltaF2` against the original P205 `DeltaM`, and documents aligned delete-one dependency groups.
