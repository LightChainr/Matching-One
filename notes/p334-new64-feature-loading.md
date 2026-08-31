# New64 on original00: contact structure captures four fifths of clock loading

The separately retained New64 stream gives a consistent approximately **80%**
captured signed loading under the previously named four-feature model; source
energy alone captures about one half. Safe contact degree and joint-safe
abundance remain related to the own-source center response after linear removal
of the two baseline mean clocks. Only committed prefix NPZ statistics are used.

## Original00 stays separate from the original8 frozen result

The population is the original1502 N325 and1551 N425 double-R0 prefixes, zero
padded to each size's original20000-prefix denominator. This differs from the
three-cell receiver-R0 pool in root's broader original8 projection. No streams
are pooled here.

|N / receiver|Old8 contact share|New64 energy-only share|New64 contact share|New64 residual loading ×10⁸|
|---|---:|---:|---:|---:|
|325 first|68.86% ±12.43%|52.89% ±5.05%|80.36% ±5.66%|0.4976 ±0.1727|
|325 second|64.03% ±24.04%|46.35% ±7.43%|78.20% ±5.84%|0.5190 ±0.1851|
|425 first|75.80% ±13.82%|55.40% ±5.16%|79.49% ±6.24%|0.4651 ±0.1676|
|425 second|64.00% ±13.41%|49.71% ±6.77%|79.66% ±6.17%|0.3885 ±0.1351|

Every error is one original20-batch delete-one SE, refitting the same named
projection. New64 observed own-source loading is respectively `2.5340±0.2630`,
`2.3806±0.3002`, `2.2677±0.2526`, `1.9097±0.2253`, all ×10⁻⁸.
New64-minus-old8 loading changes are within one paired SE in all four rows.
The sharper approximately80% pattern is a precision improvement on shared
prefixes, not independent population replication or a detected stream shift.

Contact-minus-energy captured loading is respectively `0.6962±0.0996`,
`0.7582±0.1654`, `0.5464±0.1309`, `0.5718±0.1272`, all ×10⁻⁸. The errors use
the common joint factor. The positive residual roughly20% remains in the
result; complete four-feature closure is not established.

## Geometry after removal of the two baseline clocks

The following own-source center-response entries are
`V_FY - K_F,clock K_clock,clock^(-1) V_clock,Y`. All cross-source and lifetime
entries are retained in the artifact with the same joint factor.

|N / receiver|Joint-safe mass ×10⁸|Own energy ×10¹⁰|Own safe degree ×10⁸|Own safe loop ×10⁸|
|---|---:|---:|---:|---:|
|325 first|2.837 ±0.668|7.066 ±1.127|9.826 ±2.018|1.100 ±1.037|
|325 second|2.967 ±0.913|4.924 ±1.662|11.847 ±3.160|4.915 ±1.551|
|425 first|3.005 ±0.659|7.666 ±1.555|9.694 ±2.116|2.473 ±1.449|
|425 second|3.396 ±0.639|5.064 ±1.588|9.762 ±1.981|3.349 ±0.836|

Safe degree and joint-safe abundance are the consistent geometry relations.
Raw loop mean is less uniform; lifetime-response partial moments are more
mixed and are not promoted into the same common pattern. These are partial
moments after the two clocks only, not each feature's unique effect after all
other geometry features. Correlated descriptors do not receive causal shares.

## Fixed estimands and conditional measurement error

The model sets are unchanged from `011f50e3`: own physical source energy alone;
or `joint_safe_mass,own_score_energy,own_safe_degree,own_safe_loop` together.
Coefficients are separately estimated for old8 and New64 on original00: the
fixed item is the feature set, not an out-of-sample coefficient.

Clocks are old8 conditional means of `C=(K1+K2)/(2(N+1))` and
`W=(K2-K1)/(N+1)`. New64's physical-source integral-A response divided by `-2`
is `H_C`; integral-E divided by `-1` is `H_W`. No further normalization enters.
Signed loading is `2 Cov_00(mu_C,H_C)-0.5 Cov_00(mu_W,H_W)`, weighted by00
prevalence. The captured ratio is **not explained variance**, field count,
or closure probability.

Old8 latent clock Gram and old8 clock-response moments use ordered distinct
quartets. Old-clock mean times New64 response mean is already a conditionally
independent-stream product, with no same-quartet subtraction. Exact census
features have no conditional tail measurement error. Within-cell centering
uses distinct-prefix products. Finite-prefix corrections and all coefficients
are recomputed after deleting each original1000-prefix batch. No ridge, PSD
repair, model search, or New64 latent response-variance claim is introduced.

## Sources and output

- Old8: `375cd3a12b2b7a87d79148a59f62b95898f9e471`.
- Descriptors: `1cfa4ae892a2f7f4168e9a71690efd7a5560d4cd`.
- New64: `8ad30617b0a3076a5c01a208eb213096d8879b32`,
  `experiments/p334-mechanism-response-20260831/results-extension/prefix_statistics_N{N}.npz`.
- Allocation `93ee4e98`; model `011f50e3`; reader `5c4c9b45`.

`results/p334-new64-feature-loading/score.json` preserves zero-padded20batch
sufficient rows, all receiver/source estimates, LOO/factors, clock-partial
moments, paired stream differences and input SHA256. The isolated readout took
0.43 seconds using the managed local Python environment. It does not append to
root's frozen result. No sampling, fork gzip, determinant or finite-policy
computation was performed.
