# P337 N680 heldout reveal

Preregistration `ba4ca6f` froze four point forecasts, their source covariance,
the N680 exact geometry and a 120M/shape production contract before any N680
data existed.  The design used N340 variance and projected 3.03 standard errors
between the two-mode and free-single forecasts.

Huawei HZsCM6 completed the 80 aligned batches in 1026.0 seconds.  Every
metadata, period-matrix, Smith `(2,340)`, seed/counter and projective-birth gate
passed; stderr is empty.  The 302 MB exact raw archive remains on the server and
is committed as a hash-verified 24 MB lossless gzip.

The heldout amplitude is `-0.00216756 +/- 0.00055693`.  It is closest to the
two-mode recurrence (`-0.586` measurement SE), followed by free-single
(`+1.563`), frozen nominal single-mode (`-2.048`) and scale-neutral (`+12.110`).
Source-uncertainty-aware residuals are `-0.295`, `+0.747`, `-2.024`, and
`+7.159` predictive SE.

Actual variance was 1.99 times the frozen projection, leaving only 2.15 SE
between the first two point forecasts.  The correct reading is therefore:
N680 independently favors the recurrence shape and strongly excludes no decay,
but it does not decisively separate recurrence from a free single transfer.
The projective scalar remains null at `-0.397 sigma`.
