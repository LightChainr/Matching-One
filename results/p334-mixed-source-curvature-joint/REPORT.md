# Mixed-source curvature: one original-20-batch covariance join

All ± values are paired delete-one-original-batch SE. No new sampling,
raw fork replay, weight calculation or model fitting is performed.

## N325

### New64 original00 tensor

| Receiver / observable | Own pure | Other pure | Mixed | Own − other | Own − mixed |
|---|---:|---:|---:|---:|---:|
| first / A_ref | 5.403858e-07 ± 1.248e-07 | 1.373258e-08 ± 7.727e-08 | -2.346678e-09 ± 9.494e-08 | 5.266532e-07 ± 1.475e-07 | 5.427324e-07 ± 1.765e-07 |
| first / E_ref | 7.008067e-08 ± 6.388e-08 | 4.863277e-08 ± 5.794e-08 | 1.997273e-08 ± 5.543e-08 | 2.14479e-08 ± 7.061e-08 | 5.010794e-08 ± 9.793e-08 |
| first / C | -6.222058e-08 ± 1.208e-08 | -9.81551e-09 ± 9.589e-09 | -3.109713e-09 ± 9.699e-09 | -5.240507e-08 ± 1.596e-08 | -5.911087e-08 ± 1.903e-08 |
| first / W | -5.970987e-08 ± 1.552e-08 | -1.491592e-08 ± 1.389e-08 | -5.40004e-09 ± 8.432e-09 | -4.479395e-08 ± 1.913e-08 | -5.430983e-08 ± 1.863e-08 |
| second / A_ref | 5.875932e-07 ± 1.294e-07 | -6.905566e-09 ± 9.04e-08 | 1.711465e-07 ± 9.721e-08 | 5.944988e-07 ± 1.429e-07 | 4.164467e-07 ± 1.509e-07 |
| second / E_ref | -2.709592e-08 ± 7.792e-08 | 8.626937e-08 ± 7.303e-08 | -2.700691e-08 ± 6.119e-08 | -1.133653e-07 ± 9.052e-08 | -8.901468e-11 ± 8.775e-08 |
| second / C | -5.985114e-08 ± 1.148e-08 | 4.884941e-09 ± 9.975e-09 | -2.228667e-08 ± 9.236e-09 | -6.473608e-08 ± 1.472e-08 | -3.756447e-08 ± 1.489e-08 |
| second / W | -4.745642e-08 ± 1.515e-08 | -1.195056e-08 ± 1.734e-08 | -1.123673e-08 ± 1.093e-08 | -3.550586e-08 ± 2.353e-08 | -3.621969e-08 ± 1.888e-08 |

### Paired finite-vs-zero-source and conditional-stream differences

| Observer / observable | old8 rectangle − old8 Hfs(0) | new64 Hfs(0) − old8 Hfs(0) |
|---|---:|---:|
| first / A_ref | 6.630661e-11 ± 3.167e-10 | -1.928476e-07 ± 2.192e-07 |
| first / E_ref | 9.253861e-11 ± 2.18e-10 | 1.049758e-07 ± 2.007e-07 |
| first / C | -3.62238e-12 ± 3.73e-11 | 2.658858e-08 ± 2.659e-08 |
| first / W | -1.698776e-11 ± 5.375e-11 | 8.183546e-09 ± 4.118e-08 |
| second / A_ref | 3.135924e-10 ± 3.633e-10 | 7.299927e-07 ± 2.909e-07 |
| second / E_ref | -1.150917e-10 ± 2.161e-10 | -1.540845e-07 ± 1.428e-07 |
| second / C | -2.504785e-12 ± 3.073e-11 | -6.106842e-08 ± 2.404e-08 |
| second / W | 5.439081e-11 ± 3.866e-11 | 2.065737e-08 ± 2.526e-08 |
| S / A_ref | 1.899495e-10 ± 2.754e-10 | 2.685726e-07 ± 1.903e-07 |
| S / E_ref | -1.127655e-11 ± 1.816e-10 | -2.455434e-08 ± 1.274e-07 |
| S / C | -3.063583e-12 ± 2.765e-11 | -1.723992e-08 ± 1.907e-08 |
| S / W | 1.870152e-11 ± 3.866e-11 | 1.442046e-08 ± 2.784e-08 |
| D / A_ref | 3.239033e-10 ± 5.259e-10 | 1.208767e-06 ± 4.546e-07 |
| D / E_ref | -2.719612e-10 ± 3.113e-10 | -3.393259e-07 ± 3.11e-07 |
| D / C | 1.463864e-12 ± 5.263e-11 | -1.148161e-07 ± 4.373e-08 |
| D / W | 9.349407e-11 ± 6.92e-11 | 1.633863e-08 ± 5.184e-08 |

### Own-source local gain curvature

| Receiver / observable | H1 | H2/H1 | t=1/2 quadratic / linear |
|---|---:|---:|---:|
| first / A_ref | -2.987067e-05 ± 1.777e-06 | -0.01809085 ± 0.003538 | -0.004522712 ± 0.0008844 |
| first / C | 3.027966e-06 ± 1.871e-07 | -0.02054864 ± 0.003524 | -0.005137161 ± 0.0008811 |
| second / A_ref | -2.959075e-05 ± 1.752e-06 | -0.01985733 ± 0.003926 | -0.004964331 ± 0.0009815 |
| second / C | 2.868602e-06 ± 1.633e-07 | -0.02086422 ± 0.003389 | -0.005216054 ± 0.0008472 |

## N425

### New64 original00 tensor

| Receiver / observable | Own pure | Other pure | Mixed | Own − other | Own − mixed |
|---|---:|---:|---:|---:|---:|
| first / A_ref | 1.023171e-06 ± 1.386e-07 | 1.792355e-07 ± 8.777e-08 | 3.084432e-08 ± 1.192e-07 | 8.439357e-07 ± 1.858e-07 | 9.92327e-07 ± 1.892e-07 |
| first / E_ref | -3.46965e-07 ± 8.98e-08 | -1.101381e-08 ± 8.796e-08 | 6.968934e-10 ± 5.754e-08 | -3.359512e-07 ± 1.132e-07 | -3.476619e-07 ± 1.109e-07 |
| first / C | -8.493126e-08 ± 1.014e-08 | -8.184502e-09 ± 7.148e-09 | 8.331824e-09 ± 8.956e-09 | -7.674676e-08 ± 1.32e-08 | -9.326309e-08 ± 1.408e-08 |
| first / W | -2.393913e-08 ± 1.279e-08 | -8.751517e-09 ± 1.378e-08 | 8.933935e-09 ± 1.048e-08 | -1.518761e-08 ± 1.695e-08 | -3.287306e-08 ± 1.349e-08 |
| second / A_ref | 6.770544e-07 ± 9.73e-08 | -1.519392e-07 ± 1.01e-07 | 1.854838e-08 ± 1.187e-07 | 8.289936e-07 ± 1.385e-07 | 6.58506e-07 ± 1.575e-07 |
| second / E_ref | -1.454548e-07 ± 7.241e-08 | -2.833857e-09 ± 8.61e-08 | 3.942134e-08 ± 7.012e-08 | -1.42621e-07 ± 1.217e-07 | -1.848762e-07 ± 1.027e-07 |
| second / C | -6.51988e-08 ± 8.364e-09 | 8.023462e-09 ± 1.041e-08 | -4.182909e-09 ± 1.015e-08 | -7.322226e-08 ± 1.247e-08 | -6.101589e-08 ± 1.483e-08 |
| second / W | -2.614029e-08 ± 1.511e-08 | 1.348621e-08 ± 1.5e-08 | -1.573615e-08 ± 1.301e-08 | -3.96265e-08 ± 2.107e-08 | -1.040414e-08 ± 2.163e-08 |

### Paired finite-vs-zero-source and conditional-stream differences

| Observer / observable | old8 rectangle − old8 Hfs(0) | new64 Hfs(0) − old8 Hfs(0) |
|---|---:|---:|
| first / A_ref | 4.42856e-10 ± 4.192e-10 | 6.25846e-07 ± 3.851e-07 |
| first / E_ref | -1.928801e-10 ± 2.912e-10 | -1.98917e-07 ± 2.517e-07 |
| first / C | -1.28247e-11 ± 4.344e-11 | -1.864118e-08 ± 2.959e-08 |
| first / W | 2.470113e-11 ± 5.181e-11 | 3.645003e-08 ± 3.287e-08 |
| second / A_ref | 3.562159e-11 ± 4.165e-10 | 4.176704e-07 ± 3.82e-07 |
| second / E_ref | 2.775724e-10 ± 2.287e-10 | 2.787675e-07 ± 2.08e-07 |
| second / C | 1.679363e-11 ± 3.327e-11 | -5.068044e-08 ± 2.661e-08 |
| second / W | -1.228816e-11 ± 3.175e-11 | -7.855807e-08 ± 3.279e-08 |
| S / A_ref | 2.392388e-10 ± 2.769e-10 | 5.217582e-07 ± 2.571e-07 |
| S / E_ref | 4.234612e-11 ± 1.969e-10 | 3.992526e-08 ± 1.539e-07 |
| S / C | 1.984462e-12 ± 2.293e-11 | -3.466081e-08 ± 1.759e-08 |
| S / W | 6.206486e-12 ± 2.56e-11 | -2.105402e-08 ± 1.843e-08 |
| D / A_ref | -4.560808e-10 ± 7.011e-10 | -2.331456e-07 ± 6.376e-07 |
| D / E_ref | 5.268817e-10 ± 3.867e-10 | 5.349813e-07 ± 3.855e-07 |
| D / C | 3.317095e-11 ± 6.981e-11 | -3.588225e-08 ± 4.921e-08 |
| D / W | -4.142604e-11 ± 7.731e-11 | -1.288029e-07 ± 6.085e-08 |

### Own-source local gain curvature

| Receiver / observable | H1 | H2/H1 | t=1/2 quadratic / linear |
|---|---:|---:|---:|
| first / A_ref | -3.591292e-05 ± 1.498e-06 | -0.02849034 ± 0.003371 | -0.007122585 ± 0.0008428 |
| first / C | 3.061755e-06 ± 1.356e-07 | -0.02773941 ± 0.003142 | -0.006934852 ± 0.0007854 |
| second / A_ref | -3.053273e-05 ± 1.904e-06 | -0.02217471 ± 0.002637 | -0.005543678 ± 0.0006593 |
| second / C | 2.837708e-06 ± 1.661e-07 | -0.02297587 ± 0.002314 | -0.005743967 ± 0.0005786 |

## Boundary

Fixed commuting physical-source coordinates and normalization. Pure/mixed contrasts and ratios use new64 original00, zero padded to the original population; old8 finite rectangle minus old8 Hfs is paired within its original stream. New64-minus-old8 keeps shared prefixes/batches and is not independent population replication. The side-one rectangle is an integral of Hfs over [-1/2,1/2]^2, not identically Hfs(0). H2/H1 is a local fractional gain derivative; one-quarter of it is the ratio of specified second/first Taylor contributions at t=1/2, not a finite-response prediction. No ratios are formed for cross responses, E, or W. All old/new matrices and first means are appended jointly to172fbeb1 with consistently signed original20 LOO factors. No fit, finite weights, determinant, shape test, inverse or new MC/DP.
