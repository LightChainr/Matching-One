# P3: covariance support repair and N580 recovery without production

Date: 2026-09-12. Companion to the exact-foundations correction (#702).

## Mathematical correction

For nonsingular S and a fixed full-column-rank signed-real design V,
D=min_a (y-Va)^T S^{-1}(y-Va) has the ordinary Gaussian-reference residual
interpretation. Two coordinates/one line give exactly Fieller's contrast squared.
Neither this algebra nor correct Fieller inference creates new information.
The three-versus-two-rung comparison changes the data used, not the two-rung test.

For singular S whose nullspace is KNOWN as an exact deterministic constraint,
feasibility requires U0^T(y-Va)=0, where U0 spans ker S. Find a particular
solution a0 and a null basis Z for U0^T V. On the stochastic support whiten
with W=Lambda_+^{-1/2} U_+^T and fit W V Z against W(y-Va0). The correct df is
rank(S)-rank(W V Z), not blindly rank(S)-number_of_original_amplitudes.
If the support constraint is infeasible, there is no compatible model mean.

Concrete regressions: S=diag(1,1,0), y=(0,0,1), v=(1,0,0) was D=0,p=1 in
blob `1b6270238db586c5c48849bd607376bf3c7d14fa`; it now raises an explicit
incompatibility in exact_support mode. For v=(1,0,1), the deterministic
coordinate forces a=1, giving D=1,df=2. A matrix with eigenvalue -0.1 is
rejected rather than clipped into an apparently valid covariance.

Default `covariance_mode='strict'` refuses rank-deficient/truncated inputs.
`exact_support` is an explicit scientific assertion, NOT something inferred
from too few jackknife batches. Small positive variances cut off numerically
are not silently made deterministic. Exact-support computation is numerical
at the recorded precision/tolerance; it is not an exact rational certificate.
For empirical low rank, explicitly choose/justify a reduced observation space
or acquire adequate covariance information before using a reference p-value.

A signed line {a*v:a in R} is not a positive ray. Cone/bounded-amplitude tests
need their own reference law. Estimated covariance likewise does not supply
exact finite-sample chi-square coverage. These distinctions also correct the
broader claims in the P3 manuscript; this patch does not promote its readiness.

## N580: recover before replaying

The main-branch files `results/aspect-ladder-n580/shards/rung_r{1,2,4}.json`
retain 100 aligned `_deleted` rows, plus seed, replica offset and sample count.
The missing cov(r2,r4) therefore need not be set to zero or commissioned as a
new run. `scripts/recover_n580_covariance.py` reads these arrays, checks the
#577 alignment metadata and previously saved diagonal/cross entries, then
reports the complete covariance, raw curvature and the existing eight ray tests.
It writes only stdout: historical results and freezes remain untouched.

An additional labelled sensitivity fits a common-modulus H8 nuisance plane and
bounded common ratios |rho|<=0.2 and <=1. The bounded results are MINIMUM
RESIDUALS ONLY, not chi-square-calibrated tests. This does not identify H8,
validate a common ratio across moduli, or turn the reanalysis into new evidence.

## Validation and limits

Locally: 14 projective checks and 2 covariance-reader controls passed. The
full-rank implementation exactly reproduced all returned fields in 80
fixed-seed comparisons against the saved original module. No production data
were changed. The real-shard end-to-end command is also put in a bounded,
read-only, path-triggered GitHub Actions job (three-minute cap), so it can run
where a complete repository checkout is available. Its actual result must be
read before any numerical N580 verdict is claimed. Full repository CI is a
separate check; it was not run in this local connector-only checkout.
