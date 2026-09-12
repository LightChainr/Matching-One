# Averaging probabilities is not averaging their matching roots

Date: 2026-09-13. Completed finite identity, perturbative theorem, and exact/tiny controls.

## 1. An exact annealed identity

Let random environment probabilities be p_i=p+sigma*xi_i, with independent xi_i of
mean zero and bounded support. Keep p_i in [0,1]. Conditional on the environment,
site occupations are independent. Separate affinity and independence imply

    E_xi M(p+sigma xi_1,...,p+sigma xi_N)=M(p,...,p)          (1)

EXACTLY for every legal sigma. Indeed the expectation of every product factorizes
and E[p_i]=p. The annealed occupation law itself is homogeneous Bernoulli(p).
Thus the zero of the averaged observable remains exactly p_*.

This statement fails in general for spatially correlated xi: mean-zero marginals
alone do not factorize the environment expectation. It also changes if the noise is
added to log-odds instead of probabilities, since E[logistic(z+sigma xi)] need not
be logistic(z). Neither distinction is optional.

## 2. The quenched finite roots have a different mean

For a transitive torus let q(a) solve M(q(a)+a_i)=0. The implicit-function theorem
and the source-Hessian result give

    grad q=-1/N,
    R_ij=-M_ij/M' + M''/(N^2 M'),
    R*1=0,
    trace R = M''/(N M').                                  (2)

For any bounded mean-zero random vector xi with covariance Sigma, as sigma->0,

    E q(sigma xi)
      =p_*+(sigma^2/2) trace(R Sigma)+O(sigma^3),
    Var q(sigma xi)
      =(sigma^2/N^2) 1^T Sigma 1+O(sigma^3).                 (3)

The remainders improve to O(sigma^4) when the whole vector law is centrally symmetric.
The constants are finite-system local constants; no assertion of uniformity in N is made.

For independent unit-variance fields on all N sites,

    E q = p_* + sigma^2 M''/(2N M') + O(sigma^4),
    Var q = sigma^2/N + O(sigma^4).                         (4)

Equation (1) and equation (4) are compatible. Root extraction is nonlinear and does
not commute with averaging. No new macroscopic disorder transition follows from this.

If only k sites fluctuate independently, the leading bias is

    k sigma^2 M''/(2N^2 M').                               (5)

At the 4x4 root the all-site coefficient in (4) is

    0.020558979708451625... .

For two addressed adjacent sites the coefficient is

    0.002569872463556453... .

## 3. Exact small experiment, no Monte Carlo

Two adjacent sites have independent symmetric offsets +/-sigma, all other sites
have no offset. There are only four environments. The site-occupation probabilities
and annealed equality are checked with rational arithmetic at two p values and two
amplitudes. The roots of all four environments are computed at 65 decimal digits.

| sigma | (mean of the four roots - p_*)/sigma^2 |
|---|---:|
| 1/64 | 0.00257018011526596 |
| 1/128 | 0.00256994938117334 |
| 1/256 | 0.00256989169325377 |
| second-order limit | 0.00256987246355645 |

The root-averaging effect is directly present, although the annealed M is EXACTLY
unchanged. This is not an estimated disorder effect and does not borrow a stochastic
error bar from a deterministic calculation.

## 4. Removing the empirical field mean does not remove the curvature bias

Put xi_c=xi-(1/N)sum xi_i. For originally independent unit-variance xi,
Cov(xi_c)=I-11^T/N. Since R*1=0, trace(R Cov(xi_c))=trace R. Thus the leading mean bias
is unchanged, while the O(sigma) root fluctuation is removed exactly. The variance
of the centered-field root is O(sigma^4) under the stated bounded-support assumptions.

However xi_c has correlated coordinates. The exact annealed identity (1) no longer
applies. Mean centering is not a free operation on the physical ensemble.

For correlated environments supported on one unit-RMS Fourier pattern +/-h, the
second-order mean shift is sigma^2 q''(h)/2. The complete 4x4 Hessian supplies both
signs (stripes down, checkerboard up), so there is no universal convexity or Jensen
sign for arbitrary zero-mean spatial noise.

## 5. Practical interpretation

Distinguish a pooled probability curve, an average of roots measured in different
fixed environments, and a root estimated with random finite samples. The present
calculation compares the first two, not the third. It gives an exact control for a
pipeline that mixes them, without reinterpreting prior frozen analyses as disorder
experiments. Statistical estimator bias requires its own sampling model.

The proof uses only finite multilinearity and implicit differentiation. The formulas
are applications of standard delta-method ideas; no novelty is claimed for them.
All input polynomials and four-environment outputs are included in
results/research-control-20260912/full-site-hessian-independent.json.
