# Common symmetrizing coordinate: an exact necessary condition and an unresolved data test

2026-09-12. Follow-up to #706. The real-data calculation IS NOW COMPLETED;
it was not delegated as another production task. Its result is negative in
information: the tested commutator is unresolved, not evidence for a chart.
No novelty claim is made for the algebra.

## 1. Object and conditional lemma

Let Q_i be increasing quantile functions, F_i their inverses, and

    R_i(p)=Q_i(1-F_i(p)),       R_i o R_i=id.

For this application Q_i is the ALREADY COMBINED spin0 quantile function.
F_i is its inverse, not an orientation-weighted mixture of physical CDFs.
A nonlinear map need not commute with orientation weighting; do not transfer
the outcome to a different observable without a separate map.

Suppose one increasing coordinate psi symmetrizes all three laws, possibly
about distinct centers c_i: psi(R_i(p))=2*c_i-psi(p). Then R_i R_j is
conjugate to translation by 2(c_i-c_j), and translations commute. Therefore

    (R1 R2)(R2 R3)=(R2 R3)(R1 R2),
    R1 R3=R2 R3 R1 R2.                                (1)

Composition acts rightmost first. This proves the necessary identity without
fitting a polynomial or exponent. The evaluated residual is

    C(p)=R1(R3(p))-R2(R3(R1(R2(p)))).                 (2)

Conjugating all maps preserves equality (1), not the numerical magnitude of
(2). Zero is necessary, NOT sufficient. If R1=R2, (1) holds for ANY R3, so
a null can be structurally uninformative. Three zero evaluations also cannot
prove a functional identity on an interval.

## 2. Mathematical and numerical controls

Shared logit coordinates give R_k(p)=k(1-p)/(p+k(1-p)). For k=1,2,3 and rational
p=.3,.5,.7, (1) holds exactly with Fraction arithmetic. Replacing R3 by
(1-p^3)^(1/3) preserves its involution property but gives C(.5)=
0.01884049384875336739614... . Individual symmetrizability is not common
symmetrizability. Three mathematical controls cover these cases.

The data inverse is cheaper than a nested quantile solver. If
Q(u)=w1*F1^-1(u)+w2*F2^-1(u), with positive weights summing to one, solve

    F1(x)=F2((p-w1*x)/w2),
    max(0,(p-w2)/w1) <= x <= min(1,p/w1),

and set F(p)=F1(x). The objective is increasing. This is exactly the inverse
of the quantile combination; it is NOT the CDF mixture. Positivity is checked.
Component CDFs use the beta-mixture formula of the birth-rank histograms;
an independent 55-dps mode-anchored Bernstein sum with rational histogram
weights checks the pooled compositions and selected deleted compositions.
Three more tests check the inverse and beta/Bernstein identity.

## 3. Executed existing-data result

Run 34685079362, job 103530556772, head 2f900dd2066bd38a4d53f02e76ab32cb7b32cedf;
tested merge checkout 55d2b9b5d66918493593f3b20401707ba2bb50a7. The bounded job
passed all 13 new mathematical tests, the earlier nonlinear shape recheck,
and `python scripts/shape_reflection_commutator.py` (78.51 seconds). Raw
histograms were already committed; there was no Monte Carlo or new evidence.

For N145,N290,N725, use the three reference points p_j=Q290(u_j), u=.3,.5,.7.
When deleting N290, these reference points are recomputed, not held fixed.
Each independent size is deleted separately, with 100 batches; batch numbers
across independent seeds are not artificially paired. Both composed maps are
re-evaluated and their full 3x3 covariance is retained.

| u | C(p_j), high-precision pooled | nonlinear jackknife SE |
|---|---:|---:|
|.3|-5.642573955e-10|1.493472756e-9|
|.5|-2.733621556e-13|4.855386934e-13|
|.7|+5.482200467e-10|1.450287413e-9|

All three are below 0.6 marginal standard errors. The descriptive full
Gaussian quadratic form is 0.5746485/3. There is no resolved contradiction
of (1). The result is **UNRESOLVED_COMMON_CHART**, not a positive model choice.
The first and third coordinates are strongly anticorrelated; report the matrix,
not three independent tests. Its condition number is about 4.81e8.

The beta-path Q values agree with the earlier pooled reconstruction within
2.14e-14. Pooled commutators agree with the 55-dps Bernstein path within
1.20e-15; selected deleted controls within 8.11e-16. The largest involution
control is 1.33e-15, or 0.274% of the smallest MARGINAL standard error. These
are numerical checks, not all-direction error bounds or exact confidence
certificates. Pooled composition levels remain in [.23047,.76946]; all deletion
paths remain inside the checked [.05,.95] range.

Result with full covariance/provenance:
`results/research-control-20260912/shape-reflection-commutator.json`.

## 4. Why the null needs particular caution

Write T12=R1 R2=id+epsilon*u+O(epsilon^2) and
T23=R2 R3=id+epsilon*v+O(epsilon^2), in C2 on the tested interval. Expanding
composition gives

    T12 T23-T23 T12=epsilon^2*(u'*v-v'*u)+O(epsilon^3).

Near coincident involutions, the first variation vanishes and the commutator
is second-order small. A noisy nonzero derivative estimated away from this
singular null does not by itself justify ordinary Gaussian/chi-square
calibration. The printed reference tail is therefore DESCRIPTIVE, not a
calibrated acceptance or confidence statement. Nonlinear jackknife propagation
and higher precision address different questions and do not cure nonregularity.

Consequently a null cannot be used to pick a higher-degree chart, price a
sample top-up from a weak point signal, or claim no dynamic finite-size
correction. A common symmetrizing coordinate is also weaker than full location-
scale collapse: it imposes parity but does not pin the symmetric profile.

## 5. Scientific consequence

The N145-fixed quadratic chart removes most measured asymmetry in norm but
fails as an exact full-vector law. The present parameter-free necessary
condition does NOT rule out a more general common chart. Both statements can
be true. The right surviving question is a source/geometry-resolved separation
between analytic-coordinate contributions and genuine finite-size shape
corrections, with explicit predictions and an informative statistic. No
successive polynomial-degree search, larger-N purchase, or GPU run follows
from this delivery. #275's candidate-specific original-U map remains distinct.
