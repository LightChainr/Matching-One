# What the cylinder balance root resolves: an exponentially rare conditional law

2026-09-12. Completed conditional spectral analysis, instantiated by the proved
width-2/3/4 square-site formulas. Connects #636, #276, #337 and #582. It is
separate from the all-fixed-width geometric endpoint theorem in
`thin-torus-two-birth-limits-20260912.md`.

## 1. Local spectral hypotheses, not an all-width assumption

Suppose on a neighborhood of an interior q, at a FIXED width,

    P_0(p;m)=lambda_0(p)^m[1+O(eta^m)],
    P_2(p;m)=lambda_2(p)^m[1+O(eta^m)],                  (1)

uniformly there, with the required derivative bounds, 0<eta<1, analytic positive
simple leading branches, and

    lambda_0(q)=lambda_2(q)=lambda_* in (0,1),
    a_0=(log lambda_0)'(q)<0<a_2=(log lambda_2)'(q).

Leading coefficients are one in these formulas. This is justified locally by
the already proved width-2/3/4 event/spectral representations. It is NOT assumed
at every width. For general nonunit leading coefficients the formulas below
acquire the corresponding amplitude ratio; equality of exponents alone does
not erase it.

At width four the preceding rational spectral certificate bounds the common
subleading block below the two leading eigenvalues. Its relative Perron value
at q is about .773935. This block cancels in M but not in each probability.
Thus it may be used for the hypotheses (1), while the much faster M correction
uses a different observable rate. This distinction survives source tilting.

## 2. A sharp conditional transition inside a flat unconditional CDF

Set p=q+x/m for x in a fixed bounded set. Expanding the logarithms gives

    P_j(q+x/m;m)/lambda_*^m -> exp(a_j*x),
    M(q+x/m;m)/lambda_*^m -> g(x)
                             := exp(a_2*x)-exp(a_0*x).  (2)

The unconditional threshold CDF remains F=(1+M)/2 ->1/2. In contrast,

    P(r=2 | r!=1, p=q+x/m)
      = P_2/(P_0+P_2)
      -> 1/[1+exp(-(a_2-a_0)*x)].                      (3)

Equivalently, the conditional signed contrast M/E tends to
`tanh((a_2-a_0)*x/2)`, where E=P_0+P_2. Formula (3) is a rare-sector conditional
law. It is not the distributional threshold law F and not a claim that most
configurations undergo a critical transition at q.

This gives the precise object that a fixed-width spectral crossing organizes:
the relative costs of the two exponentially rare alternatives to rank one.
The sigmoid and the flat CDF coexist because they refer to different measures.

## 3. Quantile probability window and inverse conditioning

Since a_0<0<a_2, g is strictly increasing from -infinity to +infinity. For fixed
real z define

    u_m = 1/2 + (z/2)*lambda_*^m.

For all sufficiently large m this is an interior probability level. From (2),
monotonicity and locally uniform convergence,

    m [Q_m(u_m)-q] -> g^-1(z).                         (4)

The p window is of order 1/m, but its probability mass is exponentially small,
of order lambda_*^m. Every FIXED nonmedian quantile lies outside it and instead
obeys the endpoint laws proved in the companion note.

Let h'=a_2-a_0. At the finite matching median p_m=q+o(1/m),

    M'(p_m) ~ m*h'*lambda_*^m,
    Q_m'(1/2) ~ 2/[m*h'*lambda_*^m].                   (5)

Thus the inverse CDF becomes exponentially ill-conditioned at its median even
though that median's deterministic spectral bias can be exponentially small.
Calling a shrinking finite-root bias 'increased statistical resolution' would
be wrong here.

At width four, the preceding finite root theorem gives

    q = .591417170853138481798834101735923177964270443...,
    lambda_* = .844850251961768087272984426486720187657...,
    a_0 = -1.215224209362996573970582674...,
    a_2 =  1.249128400110595642565560891...,
    h'  =  2.464352609473592216536143565....

These are deterministic high-precision diagnostics of the supplied finite
polynomials. The q value is already in Jacobsen 2015 Table 2. None is a new
infinite-square critical probability.

## 4. Source balance: exact finite formula and its limit

Use the intrinsic rank source from #337, X=r-1, not a Potts-Q continuation:

    Z(p,s)=P_0 e^-s + P_1 + P_2 e^s,
    M_s(p)=partial_s log Z
           = [P_2 e^s-P_0 e^-s]/Z.                    (6)

For every finite torus and real s, the unique source-balanced probability p_m(s)
is determined by

    log[P_2(p_m(s))/P_0(p_m(s))]+2s=0.                (7)

The normalizer is positive and cancels from the zero condition; it does NOT
cancel from arbitrary responses. P_2 is strictly increasing and P_0 strictly
decreasing on the interior, so the ratio runs continuously from zero to infinity.

Equations (1) and (7) imply for each bounded fixed s,

    p_m(s)=q-2s/[m*h']+O(m^-2+eta^m/m).               (8)

The exact finite differential identity at s=0 is

    d p_m(s)/ds |_(0) = -E(p_m)/M'(p_m)
                        ~ -2/(m*h').                 (9)

For width four, `m*dp_m/ds -> -.8115721720631601937...`.
The numerical reader solves (7) using the full finite probabilities at s=+/-1/2
and m=16,64,256, without substituting the limiting formula into the finite answer.
For example m=256 gives m(p_m(.5)-q)=-.405935683764... versus the limit
-.405786086031..., and the opposite source gives +.405630373026....

An important scope point: at s=0 the common B16 trace cancels from M, but at
nonzero s the unnormalized numerator is

    e^s tr(B15^m)-e^-s tr(B5^m)
      -2 sinh(s) tr(B16^m)+4 sinh(s) t^(2m),           (10)

before division by (1+t)^(4m). The slower shared mode is restored. The original
s=0 remainder rate therefore cannot silently be used for all source responses.
No locality or continuum-field interpretation of s is claimed.

## 5. A large-deviation source limit, without calling it criticality

For fixed p in the neighborhood where (1) holds, set s=m*sigma. Since P_1->1,
log-sum-exp gives

    lim_m (1/m) log Z(p,m*sigma)
       = max{0, log lambda_0(p)-sigma,
                log lambda_2(p)+sigma}.              (11)

At q, the two activation points are +/-[-log lambda_*]. For bounded unscaled s,
Z tends simply to 1, as the geometric theorem already proves. The piecewise
linear limit (11) describes reweighting exponentially rare rank sectors, not
the original iid percolation measure's two-dimensional phase transition.

In particular, the equality of the two rare-event rates at q has real content,
while the unscaled source symmetry Z->1 holds throughout the thin-cylinder
interior and has no threshold-identifying content on its own.

## 6. Statistical consequence — deliberately restricted to one estimator

Suppose n independent Bernoulli configurations are sampled at ONE p, and only
X=r-1 is measured per configuration. Then exactly

    Var(X)=E-M^2,
    P(all n observations equal zero)=(1-E)^n.         (12)

To see at least one nonzero X with probability at least .95 requires

    n >= ceil[log(.05)/log(1-E)].                      (13)

This is only a one-event criterion, NOT enough to estimate a root accurately.
In the local Gaussian inverse-response regime, with nE large,

    se(p_est) approximately sqrt(E-M^2)/(sqrt(n)*M'),
    se(p_est) at balance ~ sqrt(2)/[sqrt(n)*m*h'*lambda_*^(m/2)]. (14)

For nE not large, a normal root-error reference is not justified. At q4 the
full finite formulas give the following numerical diagnostics:

| m | E(q4) | approximate n for 95% chance of one nonzero X |
|---:|---:|---:|
| 8 | .45227368437 | 5 |
| 32 | .009075329013 | 329 |
| 128 | 8.488607461e-10 | 3.529120985e9 |
| 512 | 6.490166177e-38 | 4.615802110e37 |

These rows use q4, NOT the exact finite p_m; the difference is immaterial to
the asymptotic but is recorded rather than hidden. Numerical sample values are
not interval-certified integer ceilings. The formula (13) is exact, while the
listed inputs and rounded evaluations are high-precision diagnostics.

This is not an algorithm-independent lower bound. It does NOT apply unchanged
to Newman-Ziff permutation histograms, conditional/binomial reconstruction,
Rao-Blackwell estimators, rare-event importance sampling, or exact transfer
algebra. In particular the existing N145/N290/N725 productions are different
geometries and different estimators. No new acquisition is priced or requested.

## 7. Numerical precision and executed result

`thin_torus_spectral_controls.py` evaluates only small companion-matrix powers,
so m=1,000,000 is NOT a million-row simulation. All outputs are numerical
controls, not replacements for the proofs. At 80 decimal digits, three probes
per width are re-evaluated at 110 digits; the maximum absolute difference in M
is below 1e-70. This checks numerical precision, not independent physics.

At q and large m the true residual is below the working precision relative to
the two leading terms. Such M entries are explicitly null/unresolved, rather
than interpreting subtraction noise as a physical finite-length correction.
Positive P_0/P_2 and their activity remain separately evaluated in this regime.

The output retains endpoint expansions, central profiles, rare-sector logistic
values, full source-balance roots, source large-deviation checks, quartiles,
and the narrowly defined snapshot-cost diagnostics. No fitted exponent or
sampling budget enters any of these calculations.
