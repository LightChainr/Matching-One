# Charged activity/net decomposition of the existing N65 archive

The frozen transform is `J_minus=birth-exit` and `J_plus=birth+exit`; no source or data block changed.

| channel | coordinate | contrast | SE | z | marginal quadratic |
|---|---|---:|---:|---:|---:|
| A | W | 0.00214027619098 | 0.00238 | 0.900 | 0.8096 |
| A | J_minus | 0.0970158321947 | 0.0286 | 3.390 | 11.49 |
| A | J_plus | 0.0268309228011 | 0.0502 | 0.535 | 0.2859 |
| D | W | -6.56922225705e-05 | 0.00137 | -0.048 | 0.0023 |
| D | J_minus | -0.0105921122137 | 0.00904 | -1.172 | 1.374 |
| D | J_plus | -0.012851171001 | 0.038 | -0.338 | 0.1145 |

A full triplet remains `12.15 / 3 df`; net alone carries `11.49`, while activity alone carries `0.2859`.
Conditioned on `(W,net)`, activity adds only `0.2361`. A is therefore a net-timing counterflow, not a common-activity amplitude.

D gives `1.509 / 3 df`; neither net nor activity resolves an orientation response.

## Covariance eigenmodes

A, standardized `(W,J_minus,J_plus)` correlation modes:
- lambda=0.5447, vector=(+0.717, -0.669, -0.197), quadratic contribution=5.477.
- lambda=0.9675, vector=(-0.076, -0.356, +0.931), quadratic contribution=0.628.
- lambda=1.488, vector=(+0.693, +0.653, +0.307), quadratic contribution=6.048.

D, standardized `(W,J_minus,J_plus)` correlation modes:
- lambda=0.1574, vector=(-0.676, -0.117, +0.727), quadratic contribution=0.0374.
- lambda=0.7634, vector=(-0.390, +0.895, -0.219), quadratic contribution=1.197.
- lambda=2.079, vector=(+0.625, +0.431, +0.651), quadratic contribution=0.2748.

## #334 crosswalk

`J_plus` is the common positive activity and `J_minus` the source-sink residual. The A result matches the counterflow morphology: the large common activity cancels between orientations and information remains in the smaller net timing current.

The sectors and dependency groups remain distinct. This A/B1 coordinate is intra-axis; the exact #334 N13/N17 result is axis-versus-diagonal H4 orbit composition. The exact tables provide a mechanism crosswalk only. The conditional N65 line-sorting score reuses this same 1714141 block and is not additive evidence.
