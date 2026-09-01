# Prefix-local rank of a common Euler-invisible response

Frozen before this analysis is run on Huawei HZsCM6, 2026-08-31.

The completed `73608ba9` result establishes rank of an ensemble-mean response map. Its subsequent interpretation at `d74dbbf1` explicitly leaves prefix-local rank open: mixing one-active-source prefixes in cells01/10 can produce a full-rank ensemble mean. The completed common-label signed-histogram analysis at `4db356e1` and covariance alignment at `b582015e` are not repeated as new discoveries.

We use only the existing40k original prefixes (20k per N), eight independent next-label quartets per prefix from `e32a8593`, and their saved contact coordinates from `959a7fa2`. The generator uses distinct per-prefix/quartet/group/replica streams. Labels U,V are sampled independently with replacement; each has two independent suffixes. Both orientations share each suffix. No new Monte Carlo, label census, fitting, or source-centering estimate is introduced.

The inputs are `L_first=1[rank_first=0](e_first-c_first)` and the analogous `L_second`; `g_plus=(L_first+L_second)/2` and `g_minus=(L_first-L_second)/2`. We use the original joint-safe, equal-joint-degree class policy, with fixed class mass and exponential tilt `exp(t*pi_class*L_source)`.

For each prefix Z and quartet q define

```
X_q[o,s] = 1[same joint safe/degree class]
           * (L_s(U)-L_s(V)) * (mean_tail Y_o(U)-mean_tail Y_o(V))/2.
```

This is unbiased for `J(Z)[o,s] = sum_a pi_a² Cov(L_s,E[Y_o|Z,u]|a)`. Pair differencing cancels the exact unknown class mean; no shared empirical mean across quartets threatens conditional independence.

Writing X as [[a,b],[c,d]], use

```
D(q,r)=(a_q*d_r+a_r*d_q-b_q*c_r-b_r*c_q)/2.
```

The mean of all28 distinct-quartet pairs is unbiased for `det J(Z)`. Averaging all70 four-subsets and each subset's three disjoint pairings `D(i,j)D(k,l)` is unbiased for `[det J(Z)]²`. The second target avoids cancellation of opposite determinant signs across prefixes. It is not `[det E J]²`.

The designated readouts are A/E at the original `p_ref=.59274605079` and over the uniform p integral; no curve scan or new threshold selection. The source Gram matrix is a simultaneous support/control readout. All nine rank cells remain in the original denominator, with their additive contributions retained; local determinants vanish outside00 because at most one source column is active. Individual prefix values, all20 batch vectors, delete-one values, and joint covariance factors are retained. The analysis must reproduce the previous per-batch orientation-basis J values from the saved S/D vectors before any output is accepted.

The nonnegative population determinant-square has an unbiased estimator that may be negative. No clipping, nominal Gaussian p-value, field count, or local-rank-one conclusion follows from a nonresolved estimate. Especially, if four distinct informative quartets rarely coexist within a prefix, zero sample values and zero batch SE do not imply zero population response. Report active-quarter and noncollinear-pair counts, usable four-quartet pairings, and their distribution across original batches before drawing a mechanism conclusion.
