# A concrete clock forecast exists, but this U comparison is not yet affordable

For the existing one-step source, the rigid two-birth-clock model predicts
`V=dU/dsource = +0.00980 ± 0.02328` at N260 and
`+0.01753 ± 0.02817` at N340. Errors include the **paired uncertainty of the
calibration and baseline response curves**, previously missing from the
point-forecast budget. Both prediction intervals include angular screening
`V=0`. These are discovery-calibrated predictions, not new production results.

## Specified model and conversion to original U

The execution team's [finite-model derivation, `bde1a51c`](https://github.com/LightChainr/Matching-One/blob/bde1a51ca95c74448265b670ba0d9a0d87915479/notes/p154-prospective-birth-clock-transmission-decision.md)
sets `e_g=−tau1 F1'_g`, `x_g=−tau2 F2'_g`, with geometry-common shifts
constant in p near the common matching root. Write `delta=tau2−tau1`,
`D=bar(q')`, `c=bar(E')/D`, `B=P4(E')`, `H=P4(E'')`,
`T=bar(q'')` and `R=d〈rank1〉/dsource` at the moving root. Then

```text
delta = 2R / [D(1−c²)]
K_rel = N^(13/8)/(4D) * [cH−P4(q'')−(B/D)(cT−bar(E''))]
V_clock = delta K_rel
```

The common shift cancels from root-comoving U. All quantities use the same
source-block root: no high-statistics root or curvature is spliced into it.

## Joint uncertainty and prospective budget

Each N uses its original1M permutations (100 aligned batches, each containing
the original1000 plus9000 increment). For each deletion, the baseline jets,
saved lagged-source R and V, and saved root use the **same omitted batch**.
The full18-coordinate covariance and100 deletion vectors/factors for each of
the two independent N groups are in [latest.json](latest.json). This is the
same discovery population; recalibration does not create confirmation data.

For a reference **8M new permutations per N**, unchanged estimator efficiency
would give marginal `3SE` resolutions of `0.766` and `1.254`, much larger than
the point forecasts. Ignoring calibration uncertainty, three-standard-error
separation from zero would require about **48.9 billion / 41.0 billion** new
permutations. Including fixed discovery-calibration uncertainty gives

```text
M_new = 10^6 SE_old(V)^2 / [(V_clock/3)^2 − SE(V_clock)^2]
```

Neither denominator is positive. The JSON therefore records a null budget
with its reason, not an infinite-looking count or a production recommendation.
The corresponding three-SE prediction intervals are
`[−0.0600,+0.0796]` and `[−0.0670,+0.1020]`.

This planning calculation assumes independent future data, fixed efficiency
and Gaussian error scaling. It is not an exact power certificate. The old
target/calibration covariance is saved but not credited to an independent
future block. Better calibration, variance reduction or a different explicit
transmission mechanism could change the conclusion; #154 remains active.

## Decision and execution

Do not launch an8M block promising to distinguish this tiny U forecast from
zero. The useful next deliverable is a discriminating, affordable prediction
for the original U, not another lag/descriptor scan. The companion
[two fixed clock laws](../p154-fixed-clock-models/REPORT.md) uses existing
entry/completion readouts to remove one concrete angular law before production.

Executed once locally at code commit `a93ae4ee`, Python3.11/ARM64, NumPy as
recorded in JSON, **0.615 seconds**. No new random samples, configuration
replay, root solver, server action or test suite. Contract and input hashes,
full code SHA and execution commit are recorded in the output.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 scripts/analyze_p154_clock_transmission_budget.py
```

The script preserves an existing result; reproduction requires a separate
checkout/output copy. This result is in Draft PR #267, not main.
