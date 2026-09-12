# N580: seven exclusions survive independently varying bounded H8 leakage

2026-09-12. Direct analysis following #703; same-block retrospective C2, not
new sampling, not a revised freeze, and not a continuum identification.

## Question and answer

A uniform bound on |A8/A4| does not imply the SAME ratio on different moduli.
The common-ratio sensitivity in #703 and the older point-reconciliation column
therefore leave a genuine question: do the exclusions survive rung-specific
ratios? Yes, for B=0.2 and even B=1, against the explicit simultaneous mean
region below. Bare aspect remains compatible, already at B=0.

The input is #703's recovered mean and full covariance, as committed at
b924b2af94ebb99a058e252df6884594c51da381. Local computation used the published
vector/matrix/rays, not a new local shard replay. All eight B=0 distances
reproduce that artifact. The script reads the committed input directly in a
complete checkout. Historical production data, scores and freezes are unchanged.

## Model: weaken the unnecessary assumption, not the observation definition

For each existing positive shape v, let

    mu_i = a v_i (1 + lambda_i rho_i),  |rho_i| <= B,
    lambda = (1148/21025) * (-1,+1,-1),  a in R.

The three rho_i may differ. A4 still obeys the nominated shape a*v; the only
allowed additional angular contribution here is the declared H8 leakage.
For 0 <= B < 21025/1148, the amplitude sign is preserved at every rung.
For a>=0, the model image is the cone over the eight vertices
v_i*(1+lambda_i*B*t_i), t_i in {-1,+1}; a<=0 gives its negative. An arbitrary
convex combination of vertices gives exactly all independently bounded rho_i.
Thus both amplitude signs are covered, without treating a signed line as a
positive ray. Opposite leakage signs do not reduce this larger set: independent
symmetric bounds absorb each sign. Equality of the ratios would be extra physics.

The nearest point on either cone can be represented with at most three
linearly independent generating rays in R^3: a conic representation with more
than three can be reduced along a linear dependence until a coefficient is
zero. At the optimum, the positive coefficients satisfy the least-squares
normal equations on their active span. Enumerating supports of sizes 1,2,3
and the origin therefore includes the global optimum. Both sign cones are
solved, not just whichever sign looks favourable in the observed data.

The independent numerical optimality check re-solves each selected support at
60 dps and verifies q=S^-1(y-mu), g_j^T q<=0 on ALL eight rays, and mu^T q=0,
for both signs. These are sufficient convex projection conditions. This is a
high-precision numerical check, not an exact rational certificate.

## A reference region, not an invented cone chi-square law

Use one three-dimensional Gaussian mean ellipsoid

    E(y) = {mu: (y-mu)^T S^-1 (y-mu) <= 14.156413609126687}.

Its nominal tail alpha is erfc(3/sqrt(2))=0.002699796063260189. The cutoff is
chi-square with 3 df because this is the full MEAN confidence region, not the
residual df of a fitted cone. With known Gaussian covariance it has the stated
coverage. With this estimated jackknife covariance it is nominal/asymptotic;
no exact finite-sample coverage is asserted.

The simultaneous statement needs no guessed cone reference distribution:
on the event that E contains the true mean, every model image containing that
mean intersects E. Hence declaring any image disjoint from this SAME region
cannot falsely exclude a true image on that event. This applies simultaneously
to all eight models and nested bounds. It is conservative relative to an
individual 2-df line test, and is not its replacement or a change to the freeze.

## Executed result

| Shape | min D at B=0 | min D at B=1 | first B intersecting E |
|---|---:|---:|---:|
| bare aspect | 10.86446 | 5.72524 | 0 |
| no modulus dependence | 89.85851 | 73.37075 | 6.55584 |
| area scaling | 57.51249 | 50.29088 | 5.71612 |
| Q4 weight-4 shape | 55.58771 | 48.31480 | 5.49417 |
| weight-12 E12 | 132.26732 | 131.64755 | 17.33571 |
| weight-12 E4 cubed | 132.22612 | 131.60161 | 17.30638 |
| weight-12 delta | 513.65922 | 511.91617 | 18.31398 |
| weight-8 E8 | 115.80438 | 113.36508 | 14.29176 |

The artifact also reports B=0.2. Seven images remain disjoint even at B=1.
In particular, under the stated shape-plus-H8 model a Q4 explanation needs
at least one |rho_i| of about 5.49 before it can reach this nominal region;
area scaling needs about 5.72. This is not a measurement of H8, does not
exclude other angular terms or observer errors, and does not refute the Q4
module as a possible contributor to a different observable.

Do not confuse these uncertainty-aware intersection bounds with fitting the
observed point exactly. For positive z_i=y_i/v_i and common |lambda_i|=lambda,
point interpolation first occurs at

    B_point = (max z - min z)/(lambda*(max z + min z)).

This follows by choosing a=(max z+min z)/2 and equalizing the two extreme
relative residuals. Bare aspect has B_point=4.30327, yet B=0 already intersects
E: the point estimate's mismatch is not a rejection at this coverage level.
The earlier common-ratio point requirement must not be promoted into a lower
confidence bound on contamination.

## Validation and reproduction

Three inexpensive mathematical controls passed: an identity-covariance toy
with exact projection (42/11,14/11,14/11) and D=2/11; B=0 versus GLS; sign
reversal and nesting. All 24 fixed-bound fits and both endpoints of each
positive B threshold were checked at 60 dps on both cones. Fixed-bound maximum
relative KKT/complementarity residuals were 1.76e-61 and 5.10e-61. Threshold
brackets are numerical 38-step bisections, not exact intervals.

    python scripts/n580_rungwise_leakage.py --summary --source-revision b924b2af94ebb99a058e252df6884594c51da381
    python -m unittest discover -s tests -p 'test_n580_rungwise_leakage.py'

Omit --summary for fitted means, nuisance witnesses and branch distances.
--output refuses to overwrite an existing file. NumPy/mpmath are sufficient.
The compact generated result is
results/research-control-20260912/n580-rungwise-leakage-summary.json.

## Decision

No new N580 sampling or covariance replay is justified by this question.
The seven exclusions are not merely a common-ratio artefact. Bare aspect is
not established; it remains a surviving finite model at the declared evidence
resolution. Further angular acquisition belongs to #589's identifiability
question, not an attempt to push this same block across a preferred cutoff.
The independent next analyses remain #622's within-model shape motion and
#681's restricted finite-torus/cylinder dictionary. Neither needs a GPU.
