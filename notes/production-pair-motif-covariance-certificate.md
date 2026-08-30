# Exact fixed-K covariance on declared motif-control pairs

This certificate applies the generic overlap-union oracle to the five declared same-order Gaussian geometry pairs from Issue #40. It is deterministic and geometry-only: no Monte Carlo or production target values are read.

## Exact result

For each pair, the four controls are the first-geometry minus second-geometry counts of `nn_edge`, `diagonal_pair`, `face`, and `right_angle`. Each control has conditional mean zero for every feasible occupation count `K`.

| Gaussian pair | N | checked K values | rank-deficient K values |
|---|---:|---:|---|
| (8,1) / (7,4) | 65 | 66 | 0:0, 1:0, 2:2, 3:3, 63:3, 64:0, 65:0 |
| (9,2) / (7,6) | 85 | 86 | 0:0, 1:0, 2:2, 3:3, 83:3, 84:0, 85:0 |
| (11,3) / (9,7) | 130 | 131 | 0:0, 1:0, 2:2, 3:3, 128:3, 129:0, 130:0 |
| (12,1) / (9,8) | 145 | 146 | 0:0, 1:0, 2:2, 3:3, 143:3, 144:0, 145:0 |
| (13,1) / (11,7) | 170 | 171 | 0:0, 1:0, 2:2, 3:3, 168:3, 169:0, 170:0 |

Across all 600 fixed-K cases:

- conditional-mean failures: 0;
- covariance-symmetry failures: 0;
- negative principal-minor failures: 0;
- compact-histogram versus generic-oracle failures at the 15 representative K values: 0.

The covariance matrix is full rank for every `4 <= K <= N-3`. At `K=2` it has rank 2, and at `K=3` and `K=N-2` it has rank 3. It vanishes at `K=0,1,N-1,N`.

For the N=130 and N=170 pairs, the signed union histogram for `nn_edge` versus `diagonal_pair` cancels identically, so those two controls have exact covariance zero for every K. The other three pairs do not share that cancellation.

The JSON artifact stores each upper-triangular signed union histogram. Every covariance entry at arbitrary K is therefore reconstructed exactly as a short rational sum of `C(K,u)/C(N,u)`; representative K=4, midpoint, and N-4 matrices are included as explicit rationals.

## Boundary

This establishes only the control/control covariance forced by the declared finite geometries. It does not measure covariance with the production orientation target, fit a control coefficient, estimate variance reduction, measure wall time, or satisfy the `>=2x` promotion gate. Issue #40 therefore remains open.

