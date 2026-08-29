# K1/K2 activation-curve node spectrum

This is a retrospective full-curve analysis of existing threshold-rank histograms. It generates no Monte Carlo samples and fits no exponent.

The exact orientation-normalized curves are `Di(p)=(Fi_first-Fi_second)/Delta cos(4 theta)`. Their areas obey `Ai=-Delta E[Ki]/((N+1) Delta cos(4 theta))`.

| N | A1 | A2 | D2(p_bar) | D2'(p_bar) | nearest stable K2 node | K2 critical branches | local reading |
|---:|---:|---:|---:|---:|---:|:---|:---|
| 65 | 4.584615e-04 | 2.790013e-04 | 2.756718e-04 | 8.931350e-03 | none | scoreable (0) | no_local_negative_with_resolved_nearby_upper_node |
| 85 | 2.819234e-04 | 1.868931e-04 | 1.443149e-04 | 6.851718e-03 | 0.561916408 | scoreable (1) | no_local_negative_with_resolved_nearby_upper_node |
| 130 | 1.452813e-04 | 9.865125e-05 | 7.812570e-05 | 4.698347e-03 | 0.569509936 | scoreable (1) | no_local_negative_with_resolved_nearby_upper_node |
| 145 | 1.274025e-04 | 8.446455e-05 | 7.886284e-05 | 4.137386e-03 | 0.561979384 | scoreable (1) | no_local_negative_with_resolved_nearby_upper_node |
| 170 | 1.070236e-04 | 6.643965e-05 | 7.851370e-05 | 4.036258e-03 | 0.562761956 | scoreable (1) | no_local_negative_with_resolved_nearby_upper_node |
| 185 | 8.155994e-05 | 4.998713e-05 | 6.455317e-06 | 2.873397e-03 | 0.590401161 | scoreable (1) | no_local_negative_with_resolved_nearby_upper_node |
| 265 | 4.843081e-05 | 2.816318e-05 | -3.885980e-06 | 1.980366e-03 | 0.594623297 | scoreable (1) | local_negative_lobe_before_nearby_zero_not_integrated_reversal |
| 290 | 4.366401e-05 | 3.053172e-05 | 5.164857e-05 | 2.274984e-03 | none | scoreable (0) | no_local_negative_with_resolved_nearby_upper_node |
| 325 | 3.223647e-05 | 1.654196e-05 | -5.808382e-06 | 1.504608e-03 | 0.596348667 | scoreable (1) | local_negative_lobe_before_nearby_zero_not_integrated_reversal |
| 425 | 2.260872e-05 | 1.224101e-05 | -4.431895e-06 | 8.846312e-04 | 0.596942957 | scoreable (1) | local_negative_lobe_before_nearby_zero_not_integrated_reversal |

## What the full curves add

All 10 scoreable archives have the same integrated K2 direction, `positive`. The point value `D2(p_bar)` is negative only at N=265,325,425. At each of those sizes the stable matching-root-window branch has a point estimate just above `p_bar`; the negative point is therefore a local lobe next to a zero crossing, not a reversal of the integrated K2 response. Node-position standard errors remain in the JSON, so this does not claim significant ordering of the node and `p_bar`.

The JSON stores every Bernstein coefficient, so `D1(p)` and `D2(p)` can be reconstructed over the entire unit interval. Endpoint zeros are structural and excluded from the node score. Full-domain and critical-window branch stability are reported separately; whenever a delete-one replicate changes the branch count, that node spectrum is marked `not_scoreable` instead of silently aligning different roots.

## Center-width identity

With `C=(K1+K2)/(2(N+1))` and `W=(K2-K1)/(N+1)`, the exact identity `K2/(N+1)=C+W/2` splits every A2 area into center and width contributions. The per-size residuals of this identity and of the rank-area formula are serialized under `identity_audit`.

## Dependence boundary

Both directions are deleted as one aligned batch. N=65,85,130,170 also share a counter stream, so the JSON contains their cross-N covariance for `A2`, `D2(p_bar)` and `D2'(p_bar)`. Those four views are one dependency block and must not be added as independent evidence.

This finite-archive node map fits no free exponent and does not by itself identify a continuum operator.
