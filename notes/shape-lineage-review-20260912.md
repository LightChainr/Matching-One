# Review of #706 and a smooth-coordinate diagnostic

2026-09-12. This note supersedes the width standard errors and interval-to-interval norm errors of `notes/issue-622-quantile-shape-lineage-20260912.md`. The v1 source artifact is preserved. The reviewed reader is `scripts/shape_lineage_review.py`; it needs only stored pooled Q and covariance, not another production or histogram replay. Actual numerical output must be read before a result is claimed.

## Two corrections

For B=100 delete-one widths, the jackknife variance is (B-1)/B times the centered sum of squares, not the independent-sample-mean variance. The v1 width SE is too small by B-1=99. Because W=Q(b)-Q(a), its corrected variance is exactly S_Q[a,a]+S_Q[b,b]-2*S_Q[a,b], using the stored jackknife Q covariance.

For independent size blocks X1,X2,X3, Delta1=X2-X1 and Delta2=X3-X2 have cross-covariance -S2. Let g_i=Delta_i/||Delta_i||. The variance of ||Delta2||-||Delta1|| is g1^T S1 g1+(g1+g2)^T S2(g1+g2)+g2^T S3 g2. Use this for Z and A; the individual adjacent differences and per-size A norms remain valid.

Full inverse chi-square <= the diagonal sum is NOT an inverse-robustness certificate. Marginal standardized coordinates can independently show nonzero response under their own assumptions, but they do not validate a full inverse's rank/calibration. The reviewed reader uses a 60-decimal solve and labels the estimated-covariance reference nominal. sqrt(N) here is a Euclidean period length, not an integer count of sites on a winding path. Square-site tiny-rank controls share the F=(1+M)/2 mathematical observable at their own geometry; bond controls are a different model, not a rung of this site lineage.

## A deeper alternative to fitting another exponent

Observed A decreasing with size need not by itself identify a new irrelevant field. Assume, only for the following conditional lemma, Q_N(u)=h(t_N+s_N z_N(u)), with z_N(1-u)=-z_N(u), z_N(a)=-1/2, z_N(b)=1/2; h is a common C4 increasing coordinate map with h' bounded away from zero, bounded grid z_N, and s_N->0. Set X_N(u)=[Z_N(u)-Z_N(1-u)]/2. Taylor expansion gives

A_N(u)/W_N = K_N [X_N(u)^2-1/4] + O(W_N^2),
K_N=h''(t_N)/h'(t_N)^2.

Proof: W=h' s+O(s^3); the reflected numerator is h'' s^2(z^2-1/4)+O(s^4); and X=z+O(s^2). Divide by W^2. Uniformity requires the stated derivative/bounded-grid assumptions. This is a conditional analytic-coordinate lemma, NOT a percolation scaling theorem.

The median eliminates K without fitting a grid direction: K_mid=-4*A(1/2)/W. The remaining three independent coordinates test R=A/W-K_mid*(X^2-1/4). R=O(W^2) is allowed, so rejecting R=0 at finite N does not refute all smooth coordinate explanations. A/W exact equality across sizes is separately a finite model, not a theorem.

The reader also fixes ONE exact quadratic coordinate phi(p)=p+beta*(p-.5)^2 using only N145's median symmetry, then applies it without target refit to N290/N725. All target covariance includes shared source-beta uncertainty via the joint Q Jacobian. This is a stronger finite null than the asymptotic lemma. Rejecting it must not trigger an automatic higher-degree rescue. The experiment was already seen: this is C2, not prospective validation, and the physical occupation parameter is not redefined.

Seven local mathematical checks pass: the factor 99, shared-middle scalar variance 6 versus 4, exact quadratic normal form, order-W^2 smooth cubic remainder, and source-parameter uncertainty propagation. The real-data bounded workflow runs the existing-artifact reader only. The input Q covariance is retained in full; new nonlinear quantities use delta-method propagation, with the original nonlinear-jackknife A covariance shown as a comparison. Equal weighting is a same-block sensitivity, never an independent experiment.
