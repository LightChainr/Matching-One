# P3 correction addendum — 2026-09-12

The historical manuscript/tables are retained, but are NOT submission-ready.
This addendum supersedes the incompatible claims listed below; it is not a
claim that the entire manuscript has been rewritten or externally reviewed.

1. **Signed line versus ray.** The existing a in R model is a signed line.
A cone or a bounded nuisance image has different inference geometry. A general
one-parameter family of directions need not be a linear subspace.

2. **Singular covariance.** A pseudoinverse alone is not a test. In exact-support
mode enforce U0^T(y-Va)=0, then use df=rank(S)-rank(W V Z), with Z the null
basis of the deterministic amplitude constraints. Empirical low rank does not
establish deterministic constraints. #703 now validates PSD and refuses an
undeclared singular/truncated policy. Saturated df=0 models have exact residual
zero after numerical consistency checking; roundoff is not a rejection.

3. **Fieller equivalence.** On the same two coordinates the GLS line statistic
is Fieller's contrast squared. The three-rung gain is from additional data and
changed nuisance assumptions, not from defeating correct Fieller inference or
recovering information destroyed by the ratio parameterization itself.
A shared random stream is not required for a joint covariance: independent
runs have zero cross-covariance when independence is justified; coupled runs
require the actual covariance. Use source provenance, not matching batch IDs.

4. **Recovered covariance.** #703 recovered cov(r2,r4)=5.5621485097801135e-9
from existing shards. Bare aspect is no longer 'undetermined because covariance
is missing': the retrospective pure-line score is D=10.8644599/2,
p=.00437333, equivalent 2.84990 sigma, NOT rejected at the nominal 3-sigma
cutoff. The original frozen two-rung underpowered verdict is unchanged.
No new replay is needed for this entry.

5. **Curvature and sign.** Raw curvature z=-3.11062 is a one-coordinate contrast,
not the two-df model test. The sign of a convex v is reversed by negative a;
'negative curvature excludes every signed convex-shape model' is invalid.
An annihilator removes an algebraic nuisance; it does not guarantee adequate
power or immunity from covariance uncertainty.

6. **Common bound versus common ratio.** |rho_i|<=B for every rung does not
imply rho_1=rho_2=rho_4. Common-ratio results are labelled sensitivities.
The new rungwise bounded analysis in notes/n580-rungwise-leakage-20260912.md
removes that equality assumption. Even at B=1 seven of eight candidate images
miss a common nominal 99.73% three-dimensional mean ellipsoid. Do not assign
an ordinary chi-square law with guessed residual df to the bounded cones.

7. **Novelty and literature.** The linear algebra is classical. Do not claim
ratio testing is ubiquitous, or that no multivariate Fieller analogue exists,
from an unverified gap. Use #701/#695's retrieved sources only at their actual
PRIMARY_TEXT_READ/ABSTRACT_ONLY level. Frame the contribution around this
worked observation/design failure and its corrected inference, not a newly
invented general statistical theory.

8. **Design.** N650 does not cleanly identify angular H8 because angle and
Smith class move together (#589). Three rows cannot identify four unconstrained
columns (constant, H4, H8, Smith nuisance). Also, two unanchored model cones
both contain zero: their minimum pairwise distance is zero. A maximin sample
allocation needs a justified nonzero signal range, a fixed-SNR criterion, or
a source-anchored design objective. Do not hide this by assuming a strong
amplitude or calling an arbitrary normalization physical.

All recovered and new bounded scores are same-block C2 analyses with estimated
covariance. They do not increase the count of independent experiments, identify
a continuum field, prove a threshold law, or turn this draft into a submission.
