# Cluster fugacity moves the matching root; its angular response remains unresolved

The absolute cluster source produces a precisely resolved **pooled matching-root
shift** at both N65 and N85. After relocating that root and including the full
slope normalization, its response in the original intrinsic readout U remains
unresolved: the single two-size zero score is **chi-square 2.49185 / 2,
nominal p = 0.287675**. Root movement is therefore established on these finite
tori; an additional deformation of the normalized angular thermal response
has not yet been resolved.

These results reuse the original P40 million-counter block at each size. The
source is the common raw density `S=(CB+CW)/N`, without geometry-fitted
counterterms. The saved numbers and full covariance are in
[`latest.json`](latest.json); no new random counters were generated.

## The complete U response, not one selected term

Let `q=rank_black−1`, `E=q²`, and
`P4[X]=(X_first−X_second)/DeltaCos4`. At the root where the mean q of the two
directions vanishes, the original observable is

\[
U=\frac{N^{13/8}}{2}\frac{P_4[\partial_p\langle E\rangle]}
 {\overline{\partial_p\langle q\rangle}}.
\]

The table gives derivatives with respect to the **density-source** coupling
lambda in `exp(lambda*S)`, with one-standard-error uncertainty.

| Contribution to dU/dlambda | N65 | N85 |
|---|---:|---:|
| Direct thermal/source response | −0.00604014 ± 0.00964288 | −0.01615574 ± 0.01285426 |
| Numerator change from root movement | +0.00619012 ± 0.00287995 | +0.00400694 ± 0.00364574 |
| Denominator change from the source | −0.01083285 ± 0.00364533 | −0.00018867 ± 0.00398477 |
| Denominator change from root movement | −0.00020541 ± 0.00006948 | −0.00000323 ± 0.00006810 |
| **Total dU/dlambda** | **−0.01088829 ± 0.00889166** | **−0.01234070 ± 0.01238833** |
| Total z | −1.22455 | −0.99616 |

The four terms are correlated contributions to one derivative, not four
independent mechanism tests. For example, the N65 source-slope contribution
is about three standard errors from zero, while the complete response is
only 1.22 standard errors from zero. The complete sum and its covariance
determine the conclusion.

With `D=mean(q_p)`, `G=P4[E_p]`, `Jf=Cov(f,S)` and
`r=dp0/dlambda=−mean(Jq)/D`, the four contributions are exactly

\[
\frac{N^{13/8}}2\left[
\frac{P_4[J'_E]}D+
\frac{rP_4[E_{pp}]}D-
\frac{G\overline{J'_q}}{D^2}-
\frac{Gr\overline{q_{pp}}}{D^2}\right].
\]

## The common microscopic fugacity moves both roots by about 0.029 per unit log Q

A common cluster fugacity `Q=exp(t)` weights the integer cluster count,
`exp(t*(CB+CW))`. Therefore `d/dt=N*d/dlambda`; the standard errors scale by
the same N. This unit conversion is essential when comparing sizes.

| N | Pooled root p0 | dp0/dlambda, density source | dp0/dt, common cluster fugacity |
|---:|---:|---:|---:|
| 65 | 0.59277973 ± 0.00006722 | 0.000445921 ± 0.000001824 | **0.02898485 ± 0.00011859** |
| 85 | 0.59259822 ± 0.00007026 | 0.000340058 ± 0.000001444 | **0.02890491 ± 0.00012275** |

The root-response z values are 244.42 and 235.47. Positive cluster fugacity
decreases the pooled matching mean at the original p, so a positive p shift
restores its zero. The two bulk-normalized responses are numerically close;
two parent sizes do not establish a universal size limit. These are
lambda/t-zero derivatives of finite-volume pooled roots, not finite-fugacity
shifts or a measured derivative of the infinite-volume critical probability.

For a source `a+bK` common to the two directions, this root relocation removes
the change of Bernoulli log-odds and the numerator/denominator Jacobians
cancel: `dU/dlambda=0`. The present U result is compatible with that necessary
prediction; it does not show that absolute cluster fugacity is entirely a
thermal clock. The derivation and source-coordinate conventions are in
[`The thermal-clock source quotient`](../../notes/p40-thermal-clock-source-quotient.md).

## Original observable, new estimator on old configurations

The input block is pinned at
`291854a518b03eef4293431b89254f0f4429da53`: N65 `(8,1)/(7,4)` and N85
`(9,2)/(7,6)`, each one million configurations in 100 aligned batches at
`p*=0.59274605079`. Reobservation stores counts and first source moments by K.
The exact Bernoulli likelihood ratio supplies self-normalized importance
means and analytic p derivatives, including the normalization denominator.
The pooled root is re-found inside every aligned delete-one calculation.

This retains the original definition of U. It does **not** substitute the
fixed-matching covariance C for U, and it is not a reproduction of the older
threshold-integrated estimator: the fixed-p importance estimate is different
and noisier. First source moments determine the lambda-zero response only,
not general finite-lambda reweighting. At the leave-one-out roots the minimum
importance effective sample sizes are 989999.44 and 989990.80 out of 990000.

Uncertainty uses the complete same-batch direction covariance. The single
two-size nominal score treats the N-separated PRNG domains as independent;
the four contributions and previous views of these same counters are not
additional independent evidence. Input hashes, all derivative vectors and
the analysis environment accompany the saved JSON and
[`analyze_p40_source_thermal.py`](../../scripts/analyze_p40_source_thermal.py).

## Same source quotient on an independent archived stream

The same within-sector source `R=S-E[S|q]-b*(K-E[K|q])`, with one common
within-sector clock coefficient b, has also been evaluated on these million
configurations. No configurations or roots were recomputed; the saved full
and100 delete-one roots were reused. All numbers below are in the common
microscopic fugacity coordinate and retain the joint decomposition covariance.

| N | Raw v | Within-sector W | Three-state remainder v-W |
|---:|---:|---:|---:|
| 65 | −0.70774 ± 0.57796 | −0.31101 ± 0.48405 | −0.39673 ± 0.34815 |
| 85 | −1.04896 ± 1.05301 | −0.97277 ± 0.75961 | −0.07619 ± 0.57649 |

The independent NZ100k source archive gives N85 W=+2.08799±.88163. Its positive
2.37-SE hint does not receive same-direction support here. These are the same
physical readout estimated from different old streams and finite estimators;
the archive comparison was made after the NZ hint was seen, not prospectively.
Neither the three components nor the two estimator views establish a field
identity. [Saved six-coordinate covariance and root mappings](jet-split.json).

## Completed follow-on: the complete two-lineage source response

This report contains only N65/N85 parents, not the complete norm-4 test. The
completed follow-on reuses archived permutations for **65→130→260** and
**85→170→340**; its results are in
[`the complete-lineage source report`](../norm4-source-thermal/REPORT.md)
using100k original production permutations per size, with a distinct RNG
and complete-prefix Binomial estimator. It resolves the same-direction
finite-root response at all six sizes; fixed-source coefficients and
single-generator drift remain unresolved at that subset precision.

That calculation asks whether the fixed transfer coefficients are rigid
under one common microscopic cluster fugacity, or whether one common
generator drift explains the response. It uses the bulk derivatives
`v_N=N*dU/dlambda`, the q2/Jordan source residuals, and the no-division
two-lineage drift contrast defined in the quotient note. These are explicit
source-extension hypotheses: the unperturbed laws do not automatically apply
to a new source, and a compatible derivative cannot certify a failed baseline
law. This provides the next mechanism comparison without treating another
resolved global root shift as the answer to the angular transport question.
