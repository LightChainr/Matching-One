# Conditional noise weight in the fixed 147-prefix mixture

| Readout | Mean | Within-prefix suffix variance | Between-prefix variance | Total | Removable fraction |
|---|---:|---:|---:|---:|---:|
| canonical g(p_ref) | 0.1706204126 | 0.02092217261 | 0.003543006813 | 0.02446517943 | 85.5182% |
| integrated g | 0.3744119018 | 0.000577549865 | 0.0001104170921 | 0.0006879669571 | 83.9502% |

The uniform finite-mixture variance uses denominator147, not a sample ddof correction. Baseline readouts are binomial tails, not Bernoulli events. All two-readout covariance is retained.

This is the conditional-stratum noise weight for one fresh uniform suffix per uniformly drawn member of these 147 fixed prefixes. It is not global production variance reduction, CPU speedup, an independent new random block, or a distribution over coupled orientations. Between-prefix variance remains after averaging.
