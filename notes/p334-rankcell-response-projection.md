# Within-rank-cell response heterogeneity: projection and interpretable lower bounds

The same-cell prefix transport in dc4bb041 proves that rank-cell identity
does not exhaust the mean-response organization. A low-dimensional latent
projection can quantify what clock means and exact contact descriptors
capture. But there is an important built-in distinction: **including both
clock means reproduces the old transport exactly by projection algebra**.
The genuinely new question is how much *response heterogeneity* is captured,
and whether descriptors explain additional heterogeneity beyond those means.

## Target and fixed-cell centering

Fix one physical orientation, one frozen exact-census source mark, and one
original rank cell G. The policy changes neither G nor the law of the full
ordered prefix Z. Let

```
X=K1/(N+1), Y=K2/(N+1), C=(X+Y)/2, W=Y-X,
mu_C(Z)=E[C|Z], mu_W(Z)=E[W|Z],
r_C(Z)=H E[C|Z], r_W(Z)=H E[W|Z].
```

Here C/W are normalized rank clocks; their conditional first means also
equal the continuous-priority first means, but this does not import any
continuous order-statistic second moments.

Take r to be r_C or r_W, and fix a short feature vector
`A(Z)=(T(Z),mu_C(Z),mu_W(Z))`, or a declared subset. T comprises exact
prefix/census descriptors, not estimated suffix means. All following
covariances are conditional on G:

\[
 K=\operatorname{Cov}_G(A,A),\quad
 v=\operatorname{Cov}_G(A,r),\quad
 \sigma_r^2=\operatorname{Var}_G(r).
\]

The target is the variation of the **latent prefix response**, not the
variance of its eight-quartet estimate. Exact descriptors have no fork
measurement noise; they remain random across prefixes.

## Unbiased latent products from the same eight quartets

For each prefix i and quartet q, represent an exact descriptor by
`a_iq=T_i`, a clock feature by its baseline `b_iq,C` or `b_iq,W`, and the
response by `r_iq=h_iq,C` or `h_iq,W`. Then, conditional on Z_i, their means
are the desired latent quantities. The h values use this round's complete
census score `s=pi_a(g-mean_a g)` and
`h=(s_U-s_V)(f_U-f_V)/2`, without the earlier same-class pair mask.
Different quartets are independent conditional on the prefix; components
within a quartet need not be independent.

For M=8 define

\[
 \Omega_i(a,b)=\mathcal P_M(a_i,b_i)
 ={(\sum_q a_{iq})(\sum_q b_{iq})-\sum_q a_{iq}b_{iq}\over M(M-1)}.
\]

It estimates the latent product without the shared-quartet noise term.
The necessary raw products have a particularly small contract:

| Product | Conditional-unbiased within-prefix readout |
|---|---|
| exact T_a T_b | T_ia T_ib |
| exact T_a times a clock mean | T_ia times bar b_i |
| exact T_a times r | T_ia times bar h_i |
| two clock means | P_M(b_a,b_b) |
| clock mean times r | P_M(b_a,h) |
| r squared | P_M(h,h) |

Let L be the retained count in G and let bars denote quartet averages.
With `P_L` the analogous ordered distinct-prefix product, the
**source-population** covariance estimator for any pair of latent variables is

\[
 \widehat{\operatorname{Cov}}_G(a,b)
 ={1\over L}\sum_i\Omega_i(a,b)-\mathcal P_L(\bar a,\bar b).
\]

This supplies K, v and sigma_r^2, conditional on the cell sample count,
under the original independent-prefix source law. It requires L>=2.
Equivalently,

\[
 \widehat{\operatorname{Cov}}_G(a,b)
 =\operatorname{sampleCov}_i(\bar a_i,\bar b_i)
   -{1\over L M}\sum_i\operatorname{sampleCov}_q(a_{iq},b_{iq}),
\]

where the sample covariances use denominators L-1 and M-1. Thus
`Var(r)` is the between-prefix variance of bar h **minus** its average
quartet estimation variance/M. The same subtraction, including its sign,
is necessary for shared-clock/response measurement covariance. For exact
T the within-quartet correction is identically zero.

### Source population versus the exact empirical prefix mixture

If instead the target holds the realized L prefixes fixed and gives each
weight1/L, the latent global product must include diagonal prefixes with
their corrected product:

\[
 \widehat{\bar a_{\rm latent}\bar b_{\rm latent}}
 ={\sum_{i\ne j}\bar a_i\bar b_j+\sum_i\Omega_i(a,b)\over L^2}.
\]

Subtracting this from mean Omega gives the empirical-mixture covariance.
For this uniform fixed-cell mixture its estimate is exactly
`(L-1)/L` times the source-population covariance estimate above. The same
factor applies to K, v and sigma_r^2. A plug-in beta consequently happens
to agree, while the covariance/variance targets differ. This algebra does
not make the empirical population representative of the source law.
Keep the existing source-population convention for the dc4bb041 transport,
and retain the original cell masses when returning to the full population.

## Population projection and what the inequality actually guarantees

After removing only declared structurally constant or algebraically
redundant coordinates, suppose the population K is positive definite.
The centered best linear predictor and its residual are

\[
 \beta=K^{-1}v,\qquad
 r-E_G r=\beta^T(A-E_G A)+\epsilon,\qquad
 \operatorname{Cov}_G(A,\epsilon)=0.
\]

Therefore

\[
 \boxed{L_A:=v^TK^{-1}v\le\sigma_r^2,\qquad
        \operatorname{Var}_G(\epsilon)=\sigma_r^2-L_A.}
\]

This is both a projection variance and a lower bound on latent response
heterogeneity. With a single declared nonconstant feature,
`Var_G(r)>=Cov_G(A,r)^2/Var_G(A)`. For several responses at once the
matrix counterpart is `Cov_G(r,r) >= V^T K^{-1} V` in the PSD order.

The statements concern population moments. Although the moment estimates
above are unbiased, their inverse/product plug-ins are not unbiased and
are not statistical lower confidence limits. A noisy corrected K may be
indefinite and a corrected sigma_r^2 may be negative. Do not clip these to
PSD, turn a negative variance into zero, or manufacture an R-squared in
[0,1]. If K is not usable on the declared feature basis, report that the
multi-coordinate projection is not resolved; do not select a ridge,
eigenvalue cutoff or a different feature subset to make it work. Exact
structural zero/redundant columns may be removed with their identity stated.
Nonlinearity or an unexplained singularity requires a separately specified
analysis, not an arbitrary regularized answer.

Exact T-only K has no fork-noise subtraction and is an ordinary sample
covariance. Its cross-vector v still has finite-prefix/fork uncertainty.
A nonzero plug-in lower-bound expression is not a confidence statement.
All uncertainty stays with the original20 batch deletions, re-pooling and
re-centering each retained cell before forming nonlinear quantities.

## Reusing the old same-cell transport without overclaiming explanation

The old within-cell prefix contribution is

\[
 B_G=\operatorname{Cov}_G(\mu_X,r_Y)
       +\operatorname{Cov}_G(\mu_Y,r_X)
    =2\operatorname{Cov}_G(\mu_C,r_C)
       -\tfrac12\operatorname{Cov}_G(\mu_W,r_W).
\]

Project r_C and r_W on the same declared A, with coefficients beta_C,beta_W.
Let `k_C=Cov_G(A,mu_C)` and `k_W=Cov_G(A,mu_W)`. The reusable decomposition is

\[
 B_G^{\rm proj}=2k_C^T\beta_C-\tfrac12 k_W^T\beta_W,\qquad
 B_G^{\rm residual}
   =2\operatorname{Cov}_G(\mu_C,\epsilon_C)
       -\tfrac12\operatorname{Cov}_G(\mu_W,\epsilon_W),
\]

with `B_G=B_G^proj+B_G^residual`.

**If A already contains mu_C and mu_W, the residual is exactly zero by
orthogonality.** Clock-only projection therefore captures100% of this
particular signed transport functional by construction, even if it explains
very little of Var(r_C) or Var(r_W). Adding descriptors cannot improve that
already exact transport accounting. Reporting it as newly discovered
mechanistic explanatory power would conflate two different targets.

Two complementary, non-tautological readouts remain:

1. **Descriptor-only transport projection:** choose A=T. The displayed
   B_G^proj then measures the part aligned with those declared descriptors;
   the residual retains clock/response alignment outside their linear span.
   Projecting both the clock and response on T gives the same term, with
   their two residuals supplying the remainder. These are signed covariance
   contributions; their ratios to B_G need not lie in[0,1].
2. **Additional response heterogeneity beyond clock means:** start with
   L=(mu_C,mu_W), then add T. Define

\[
 K_{T\cdot L}=K_{TT}-K_{TL}K_{LL}^{-1}K_{LT},\qquad
 v_{T\cdot L}=v_T-K_{TL}K_{LL}^{-1}v_L.
\]

   On an identified nonredundant residual feature basis the incremental
   projection variance is

\[
 \boxed{\Delta L_{T\mid L}
      =v_{T\cdot L}^TK_{T\cdot L}^{-1}v_{T\cdot L}\ge0.}
\]

   It can be positive even though the increment in explained old B_G is
   identically zero. This is the direct question of additional response
   organization, not another accounting of the same two covariances.

A further parameter-free consequence uses the old transport alone:

\[
 \boxed{\operatorname{Var}_G(r_X)+\operatorname{Var}_G(r_Y)
 \ge {B_G^2\over\operatorname{Var}_G(\mu_X)+\operatorname{Var}_G(\mu_Y)}.}
\]

It follows by Cauchy-Schwarz on the centered vectors `(mu_Y,mu_X)` and
`(r_X,r_Y)`. If the denominator is zero, B_G is necessarily zero and the
ratio is not used. Thus a nonzero same-cell transport already entails some
latent response heterogeneity, even when the direct h-squared variance
estimate is noisy. The finite-data ratio is again not automatically a
statistical lower bound.

## Scope and handoff

All projections and covariance functionals are defined within each physical
orientation and rank cell before applying cell masses and S/D. Projecting
the paired difference directly is a different problem with cross-geometry
features/covariances. Do not multiply pooled S/D means as a substitute.

Exact census descriptors improve measurement, not causal identification:
their association with r does not establish that editing the descriptor
while holding other geometry fixed would produce that response. A finite
linear projection also is not a claim of conditional-mean closure or a
CFT field identification. No predictive test on an independent source is
created by reusing the same archive.

Scientific card: this note turns the dc4bb041 same-rank-cell residual into
latent projection and heterogeneity targets, gives the required finite-M/L
products from03603388, and separates a tautological clock-transport fit
from genuine additional descriptor-associated response variance. It is
theory only: no data pass, ML fit, tuning, validation framework, new MC/DP,
PR or issue comment. The existing quartet extraction and exact descriptor
pass can supply every stated coordinate once.
