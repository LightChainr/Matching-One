# Server result intake: pre-registered interpretation rules

Purpose: avoid changing the scientific story after seeing the server output. When new results are pushed, classify them against the following gates before tuning any model.

## A. Gaussian-integer orientation tomography

Preferred raw outputs per geometry `(a,b)`:

- period matrix `[(a,b),(-b,a)]`;
- `N=a^2+b^2`;
- `theta`, exact/rational `cos(4 theta)`;
- matching root and standard error;
- `M(p_ref)` at one common frozen `p_ref`;
- derivative or slope near the root;
- RNG counter range / sample count.

### Strong success

At multiple same-N pairs,

\[
N^{13/8}\frac{M_{\theta_1}(p_{ref})-M_{\theta_2}(p_{ref})}
{\cos4\theta_1-\cos4\theta_2}
\]

is statistically compatible with an N-independent limit or a slowly drifting common asymptote, and fixed-N multi-angle data prefer a `cos(4 theta)` harmonic.

### Partial success

Axis/diamond or same-N differences retain a stable sign but the normalized amplitudes drift substantially. Continue to larger N before claiming spin 4.

### Failure

- sign pattern is inconsistent with `cos(4 theta)`;
- same-N differences are comparable to noise or reverse unpredictably;
- `pi/8` approximants are not suppressed relative to nearby orientations;
- arbitrary orientation labels outperform the harmonic model out of sample.

Do not rescue a failed `cos4` law by fitting many harmonics on the same test set.

## B. Pell axis/diamond pairs

Keep these as an independent cross-check even if Gaussian fixed-N tomography works.

Report:

\[
\Delta_k=p_D-p_A,
\qquad
w_{\rm eff}=-\frac{\log|\Delta_{k+1}/\Delta_k|}{\log(L_{k+1}/L_k)}.
\]

The threshold-free prediction is `w_eff -> 4` if the ordinary root bias is dominated by a sign-changing spin-4 correction.

A simple average is scientifically interesting only if its held-out convergence improves without refitting weights.

## C. kappa_3

Never judge the `-5/3` hypothesis from a raw decimal alone.

For each exact-threshold control and target sequence fit, with choices frozen before the largest sizes:

1. `k(L)=k_inf+a L^-3/2`;
2. `k(L)=k_inf+a L^-3/2+b L^-2`;
3. one free-exponent model trained only on smaller sizes.

Report held-out prediction error and the transformed residual

\[
Y_L=L^{3/2}(\kappa_3(L)+5/3).
\]

Interpretation:

- finite nonzero `Y_L` limit supports both `k_inf=-5/3` and leading `L^-3/2` correction;
- linear/systematic drift in `Y_L` rejects at least one of those two assumptions;
- agreement in square bond but not square site indicates implementation/coordinate/correction issues rather than immediate failure of universality;
- disagreement between same-shape models in the extrapolated limit rejects the claimed universal ratio.

## D. Full universal crossover profile

If microcanonical data are available, preserve them. Reconstruct

\[
\rho(z)=\frac12\mathcal M'(z)
\]

and standardize the horizontal scale by a frozen convention such as unit peak curvature or unit variance when finite.

Compare the whole profile across microscopic models using:

- pointwise residuals on a common standardized grid;
- integrated squared distance;
- `kappa_3`, `kappa_5`;
- tail asymmetry (which should vanish for the matching-odd derivative profile).

A full profile collapse is much stronger evidence than a single rational-looking invariant.

## E. Matching control variates

For channels with common expectation, estimate covariance on a pilot set, freeze weights, and score on a fresh set.

Report variance ratios relative to the best single channel, not just relative to a convenient baseline.

GPU gate: require stable >=2x variance reduction or another clear covariance structure before optimizing a many-replica kernel around it.

## F. Annihilator hierarchy

For every proposed weight vector record:

- cancelled powers/log basis;
- sizes used;
- L1/L2 norms;
- covariance-weighted variance factor;
- conditioning of the constraint matrix;
- worst-case response to one omitted correction basis function.

A formally higher asymptotic order is rejected if its expected mean-square error on withheld sizes is worse.

## G. Blind finite-width predictions

The n=22..24 predictions already frozen in the repository must be scored before importing those values into any model-selection stage. Do not revise them or redefine the scoring metric after the reveal.

## H. Unexpected results

Unexpected results should be preserved before interpretation. If a server run shows a sign flip, anomalous exponent, or a value extremely close to a simple constant:

1. commit raw aggregates and metadata first;
2. reproduce with an independent seed/counter range;
3. run an exactly solved control through the same code path;
4. only then add a new analytic hypothesis.

This is especially important for unusually attractive numbers such as `-5/3`.
