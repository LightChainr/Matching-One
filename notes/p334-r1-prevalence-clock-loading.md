# R1 prevalence and conditional-clock loading: one exact symmetric split

This entry point is prepared **before consuming the final producer artifact**.
Only producer code `5b81e2f9`, its manifest, and the supplied schema were read;
no partial batch or result was inspected.

For each orientation `i` in a fixed size, define `r_i=mean(R1_i)` and
`Y_i=mean(R1-weighted hybrid clock readout)`. Missing R1 rows contribute zero
to the original all-counter denominator; they are not relabeled rank two.
For positive `r_i`, put `m_i=Y_i/r_i`. Each of the two endpoints is treated
identically: the canonical second-birth binomial tail at `p_ref`, and its
integral over `p`.

Let `delta=delta_cos4` from the final producer score. The named decomposition is

\[
D=\frac{Y_f-Y_s}{\delta}
=\underbrace{\frac{(r_f-r_s)(m_f+m_s)}{2\delta}}_{C:\ \mathrm{prevalence}}
+\underbrace{\frac{(r_f+r_s)(m_f-m_s)}{2\delta}}_{L:\ \mathrm{conditional\ clock}}.
\]

This is the symmetric product-difference identity, not a fitted model.
`C` holds the average conditional clock loading against the difference in R1
frequency; `L` holds the average R1 frequency against the difference in
within-R1 clock loading. It avoids choosing either orientation as an arbitrary
reference. Swapping orientations and the sign of `delta` leaves both normalized
terms unchanged. `C` and `L` may cancel, so signed values and covariance are
reported rather than attribution percentages.

The producer's `Y` is hybrid: exact conditional means replace solved whole
pairs; whole-pair fallbacks retain the original suffix observations. Therefore
`m_i` is an empirical conditional **hybrid** readout, not a claim that every
suffix was analytically averaged. This is only the named R1 contribution to
second birth, not full `F2`, `A_top`, or a causal attribution.

The 20 original equal-size batches retain the two paired orientations and
both endpoints. Each leave-one-batch-out replicate recomputes risks, ratios,
`C`, `L`, and `D` from common sufficient sums. Full covariance includes risks,
unconditional `Y`, conditional `m`, and both endpoint decompositions. In
particular `Var(D)=Var(C)+Var(L)+2Cov(C,L)` is preserved. Each size is reported
separately; this analysis does not assert cross-size independence or create
an independent random block.

After the producer owner supplies `FINAL_SHA`, run:

```sh
python3 scripts/p334_r1_prevalence_clock_loading.py \
  --source-commit FINAL_SHA \
  --source-directory results/p334-paired-clock-loading \
  --output results/p334-r1-prevalence-clock-loading
```

The CLI reads immutable git blobs only. A final `score.json` plus all 20
`batches/N*.batch*.json.gz` per declared size are required. It consumes `risk`
and `Y=[first_pref,second_pref,first_integral,second_integral]`; `X` is not
substituted. It reports the producer's fallback counts and provenance, and
performs no reliability solve, continuation sample, or old N900 calculation.

Scientific target: distinguish whether the observed R1-layer H4 contrast is
loaded mainly through differing R1 prevalence, differing conditional clock
responses, or their cancellation. No numerical result is claimed before the
completed committed source is supplied.
