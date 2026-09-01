# One fixed angular transmission law fails before new production

The existing one-step source already distinguishes two explicit finite-clock
laws. **A pure cos(4θ) relative clock shift is incompatible with the measured
entry/completion direction contrasts:** paired residuals are about31–35 SE
at N260/N340. The direction-independent clock law has contrast residuals
between0.19 and1.21 SE in absolute value and is not excluded by these readouts.
This removes a specific proposed response law; it does not identify a CFT field
or establish the surviving law as the global norm4 mechanism.

## Exactly two models, calibrated without directional target fitting

Let `e_g` and `x_g` be the existing source derivatives of first activation and
second completion at the pooled root, `f1g=F1'_g`, `f2g=F2'_g`.
Both models have the same form and two scalar calibration coefficients per N:

```text
e_g = −(m−delta w_g/2) f1g
x_g = −(m+delta w_g/2) f2g
scalar:   w_g=1
fourfold: w_g=cos(4theta_g)
```

Both shifts are constant in p near the root. The fourfold model is a concrete
hypothesis, not a consequence of lattice symmetry. The exact rational angular
weights are saved in JSON; no angular gain, intercept, exponent or lag is fitted.
Only the pooled moving-root rank1 response R and root displacement pdot
calibrate m and delta. Directional contrasts and U remain predictions:

```text
c_w = bar(w E')/D
S_w = bar(w q')−c_w bar(E')
delta = 2R/S_w
m = pdot−delta c_w/2
K_w = N^(13/8)/(4D) * [c_w H−P4(w q'')−(B/D)(c_w T−bar(w E''))]
V_pred = delta K_w
```

Here `D=bar(q')`, `B=P4(E')`, `H=P4(E'')`, `T=bar(q'')`; P4 divides the
first-minus-second difference by the exact cos4 difference. Source derivatives
and all baseline jets use the same source-block root. Entry includes01+02;
completion includes12+02, so their sum gives the matching-source derivative.

## Results that remove the angular law

For the **fourfold** law, observed-minus-predicted entry and completion
contrasts are:

- N260: `+0.012220 ± 0.000349` and `−0.012240 ± 0.000373`.
- N340: `−0.015162 ± 0.000485` and `+0.015349 ± 0.000500`.

Its original-U predictions are `+2.77684 ± 0.08936` / `−4.08127 ± 0.15735`.
Observed-minus-predicted U is `−2.39787 ± 0.72680` /
`+5.35780 ± 1.18263`. The entry/completion mismatch is the much sharper
diagnostic. These correlated residuals are not added as independent votes.

For the **scalar** law, entry/completion contrast residuals are
`−0.0000458 ± 0.0001388`, `+0.0000309 ± 0.0001210` at N260 and
`+0.0000263 ± 0.0001372`, `+0.0001613 ± 0.0001329` at N340.
Its U predictions are only `+0.00980 ± 0.02328` / `+0.01753 ± 0.02817`;
their independent-production budget is analyzed in the
[companion report](../p154-clock-transmission-budget/REPORT.md).
Compatibility of two root-level contrasts does not prove the full p-window
law, its cross-size continuation, or angular screening of global U.

## Dependence, execution and next action

This is a **discovery-data model comparison**, not a prospective validation.
The two model definitions and scoring script were committed before this
numerical comparison (`7f399019`), but the research program was already informed
by these old archives. Each original aligned batch is deleted from source
centering, baseline, saved matching root and model calibration together.
The full joint covariance and all deletion vectors/factors are in
[latest.json](latest.json); there are no new random prefixes or permutations.
SE ratios describe the mismatch scale and are not adaptive-search-adjusted
significance certificates.

Do not promote or independently re-run this pure-fourfold clock candidate as
a plausible P0 rival, and do not rescue it by fitting another angular
coefficient on the same archive. Keep #154 focused on an explicit transmission
law with measurable original-U consequences. The scalar law supplies a useful
low-dimensional baseline, but its present near-zero U forecast does not justify
a huge run merely to distinguish it from zero. No Issue or research route is
closed by this result.

Executed once locally, **0.938 seconds**, code/contract/input SHA256 and Python/
ARM64 environment saved in JSON. No Monte Carlo, configuration replay, root
solver, test suite or server action. Reproduce in a separate checkout/output
copy with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 scripts/analyze_p154_fixed_clock_models.py
```

Result integration status: open Draft PR #267, not main.
