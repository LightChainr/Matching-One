# C07 Stage A: leakage-safe polynomial/Padé audit

Model and fitting-window selection used only rolling folds whose test widths
ended at or before 18. Widths 19--21 were evaluated only after the candidate
was frozen. Candidates required at least two valid rolling folds, and rational
fits were rejected when a denominator pole lay within one quarter of the
observed/extrapolation `x=n^-2` interval.

Both preregistered selection rules (median and worst rolling RMSE) selected the
degree-4 polynomial correction `F(x)` with `n_min=9`. Its frozen signed errors
at widths 19, 20, and 21 were

```text
+2.0916e-12, +4.8114e-12, +8.3517e-12
```

with held-out RMSE `5.6943e-12`. It improves the earlier `n_min=8` result, but
all three residuals remain positive and grow with width. The systematic drift
is therefore reduced, not removed.

The best admissible rational candidate was Padé `[2/2]`, `n_min=8`. Its
validation median RMSE was `4.9049e-12` and its target RMSE was `8.2919e-12`,
both worse than the selected polynomial. Its target signed errors were also
positive and increasing:

```text
+3.6429e-12, +7.3335e-12, +1.1799e-11
```

Several other Padé/window combinations failed the nonlinear descent, lacked
two valid folds, or were rejected by the predeclared pole rule. These failures
are retained in `stage_a_audit.json`; they are evidence of numerical
instability, not evidence against every possible rational representation.

Conclusion: within the preregistered Stage-A families and available widths,
Padé corrections do not cure the same-sign drift more predictively than the
polynomial baseline. The selected intercept remains a deterministic model
output, not an uncertainty-calibrated estimate of the infinite threshold.
