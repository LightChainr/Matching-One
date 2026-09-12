# #622 reviewed: finite shape movement and its dominant smooth-coordinate component

2026-09-12. The original #706 source analysis is retained in Git history and
`results/probe-invariant-shape/quantile-shape-lineage-622.json` (v1). Its pooled
Q/Z/A and per-size covariance remain inputs. Its W standard errors and the
independent-error comparison of two adjacent interval norms are superseded.
Do not run the historical `quantile_shape_lineage_622.py` and interpret those
v1 fields as the current verdict. Current analysis entry points are
`scripts/shape_lineage_review.py` and `scripts/shape_lineage_nonlinear_jackknife.py`.

## Executed validation

Run 34684353419, job 103528608026, head eff8c5e158fb9be93ffde27ba80bc958eff34573,
merge checkout 51f9817c59fcf757858c0338ec724aa4992e0d94. Seven new mathematical
checks passed. The first command read existing Q/covariance; the second actually
reconstructed pooled and all 100 delete-one Q vectors from each of the three
committed histogram sets. It took 94.70 s; all three pooled Q vectors reproduced
EXACTLY in the executed float path. No simulation or new independent evidence.

Summary: `results/research-control-20260912/shape-lineage-reviewed-summary.json`.
The scripts produce full covariance/diagnostics; summary fields were extracted
from the successful job stdout, not invented or inferred from CI colour.
This bounded job is not a claim that the complete repository suite passed.

## Corrected original conclusion

A=Z(u)+Z(1-u)-1 is resolved nonzero and decreases in norm across N145/290/725.
This finite-lineage result survives. Individual adjacent Delta Z/Delta A also
remain nonzero. The width SEs must be multiplied by 99: the original script
mistook delete-one estimates for independent observations. The corrected
primary widths and SEs are

| N | W | corrected jackknife SE |
|---|---:|---:|
|145|0.1192132991|3.673003868e-6|
|290|0.0922037501|2.737674335e-6|
|725|0.0655307632|2.322561571e-6|

Delta1 and Delta2 share N290, so their cross-covariance is -S290 even though
size blocks have independent random streams. Correctly propagating it gives
|Delta Z_2|-|Delta Z_1|=-0.0003775392 with SE 0.00004177061 (about 9.04 nominal
sigma, not 10.9). The analogous Delta A norm change has SE 0.00007209953 and
is only 0.897 nominal sigma from zero. Full nonlinear deletion independently
reproduces these corrected uncertainties. No convergence or exponent follows.

## New finding: most asymmetry is consistent in magnitude with a smooth chart

The ratios ||A||/W are 0.2559635, 0.2563966, 0.2568465. The median coefficient
K_mid=-4*A(.5)/W is approximately -0.511 at all three sizes. This motivates
an analytic-coordinate control, not another freely fitted exponent.

The conditional normal-form lemma in `notes/shape-lineage-review-20260912.md`
states: if Q=h(t+s z) with a common increasing C4 h and reflected z, then
A/W=K*(X^2-1/4)+O(W^2), K=h''/h'^2 and X=(Z-Z_reflected)/2.
It is a lemma about coordinate transformations, not a percolation theorem.

One exact quadratic chart was fixed from the N145 MEDIAN ONLY:

    phi(p)=p+beta*(p-.5)^2,
    beta=0.2684067158 +/- 0.0002742089 (nonlinear jackknife).

It is increasing throughout [0,1]. Applied unchanged to N290 and N725, it
reduces the independent-coordinate asymmetry norms to respectively 0.267081%
and 0.068014% of their original values. This is norm reduction, NOT explained
variance or proof that the physical scaling field is quadratic.

The attractive scalar reading is not the full-vector verdict. Including the
shared source-beta uncertainty, the two target residual vectors have
D=105.9080 on 8 nominal degrees of freedom (log10 p=-18.5793). The source's
other three coordinates already fail exact quadratic symmetry. Likewise,
exact equality of the entire A/W vector across sizes fails (D=309.74895/8).
Thus neither exact scalar collapse nor a universally exact quadratic chart
is established. Do not raise the polynomial degree until something passes.

The full nonlinear jackknife confirms the target covariance against the
linearized calculation in ALL stochastic directions: generalized eigenvalues
lie between 0.999998211 and 1.000002384. This matters because the target
covariance condition number is about 1.3e8; agreement of diagonal errors alone
would not have sufficed. It checks propagation, not exact tail coverage for
an estimated covariance. All references remain nominal/asymptotic.

Equal weighting has the same qualitative reading but is a sensitivity of the
SAME blocks. Its quadratic residual ratios are 0.323681% and 0.052247%; its
full-vector null also fails. No second evidence vote is counted.

## Decision

The raw A signal should not be named a new irrelevant field merely because it
shrinks. Separate the dominant smooth-coordinate-like contribution from the
much smaller resolved residual. Rejecting one exact quadratic chart does not
reject every smooth chart, while arbitrary higher-degree fitting is not an
identifying experiment.

`notes/common-chart-commutator-20260912.md` gives the next parameter-free
necessary condition: three reflection involutions from the effective quantile
laws must have commuting pair compositions if one common symmetrizing chart
exists. Its real-data evaluation is not done in this delivery. It can use the
same histograms; no larger-N or GPU purchase is licensed. #275's original-U
candidate-map question remains separate.
