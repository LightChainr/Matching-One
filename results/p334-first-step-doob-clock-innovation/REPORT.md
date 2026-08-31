# P334: next direct absorption explains 7%-8% of integrated clock noise, 19%-22% canonically

This zero-DP result consumes the full forty-batch clock archive at
`0d1e586dafbade5e7d1f9bfc598170d0c881e337`. It is a new deterministic
readout of the same sampled prefixes, not independent production.

For p=h/d and direct-event terminal value a, the exact binary-event
innovation is B=p/(1-p)(a-m)^2. It lower-bounds the innovation from revealing
the complete next label. Integrated B=h/(d-h)(mu-1)^2/(N+1)^2.

| N / orientation | Mean integrated B | Mean integrated V | B/V | Canonical B/V |
|---|---:|---:|---:|---:|
| 325 first | .0000534903 | .0006474735 | 8.2614% | 22.3434% |
| 325 second | .0000536956 | .0006423861 | 8.3588% | 22.4428% |
| 425 first | .0000377827 | .0005450511 | 6.9320% | 19.3629% |
| 425 second | .0000377538 | .0005501719 | 6.8622% | 19.3264% |

The fraction is sum(B)/sum(V), not mean(B/V). Batch SEs for the integrated
fractions are .0272/.0362/.0276/.0337 percentage points; canonical SEs are
.0351/.0403/.0357/.0456 points, in table order.

Interpretation: the original direct absorbing set is a meaningful but
minor clock-information channel. Most conditional noise remains within
the first-safe event. That remainder may be resolved by safe next-label
identity and/or by subsequent insertions; it is not necessarily all noise
surviving observation of the complete next label.

Coverage: 35694 complete orientation clocks out of 35954 rank-one source
rows. Unresolved clocks are missing, not zero variance. Fifty completed
first-orientation clocks inside whole-pair fallbacks are included only in
these marginal diagnostics. All 40 original batches and their counter
intervals are retained; no child, network, or DP was reconstructed.

Boundaries: h=0 gives a zero direct-event floor; h=d makes T=1 and V=B=0,
with no prefix ratio. The archive has no solved h=d prefix. No orientation
bounds are summed into an H4 bound: that requires their common next-label
intersection and covariance. This is also not a statement about total
population variance including rank-stratum prevalence.
