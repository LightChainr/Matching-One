# P439: direct/plateau decomposition of the same-stream matching response

Retrospective real-archive analysis; **zero new samples**. This is a partition of one
microscopic source, not evidence for two independent sources.

**Result:** neither family reveals a resolved hidden loading. The current data do
not support the story that two individually strong terms cancel to hide M. They
leave both loadings unresolved: direct |z|≤1.11, plateau |z|≤1.28. This is not a
proof of exact zero or a rejection of all possible cancellation mechanisms.

## Fixed-p result

All entries are exact H4 direction contrasts at p=0.592746050790; errors are one
paired-batch jackknife SE. K_A and its delete-one vectors are reused unchanged.

| N | M_direct | M_plateau | M_total | corr(direct, plateau) |
|---|---|---|---|---|
| 85 | 1.8014743e-05 ± 0.000118933 | -0.0001673273 ± 0.00101358 | -0.00014931255 ± 0.00100547 | -0.1266 |
| 170 | 1.240264e-05 ± 1.74813e-05 | 0.0002068762 ± 0.000161921 | 0.00021927884 ± 0.000162068 | -0.0455 |
| 340 | -1.2696525e-05 ± 1.14465e-05 | 7.0724003e-05 ± 0.000136082 | 5.8027477e-05 ± 0.000136836 | 0.0240 |
| 680 | -1.0391719e-06 ± 2.83588e-06 | -3.8722828e-05 ± 4.34603e-05 | -3.9762e-05 ± 4.37148e-05 | 0.0574 |

| Component | zero across four N: chi2/4, p | adjacent wedge: chi2/3, p | loading/K_A, 95% profile | ray p |
|---|---|---|---|---|
| M_direct | 1.8909, 0.75581 | 1.8154, 0.61158 | 0.000254399, [-0.00155397, 0.00213763] | 0.61223 |
| M_plateau | 2.7236, 0.60509 | 2.3435, 0.50424 | -0.00459362, [-0.0221392, 0.0174465] | 0.47192 |
| M_total | 2.8598, 0.58155 | 2.4458, 0.48516 | -0.00466903, [-0.022222, 0.0173954] | 0.44901 |

The joint eight-coordinate (direct, plateau) zero diagnostic is 4.69005/8 df, p=0.790132.
The six independent component wedges give 4.16683/6 df, p=0.654111.
The total wedges are algebraic sums and are not additional degrees of freedom.

## Exact partition and dependence

For H_k(p)=Pr[Binomial(N,p)>=k], use the original total sample denominator:

- M_direct = E[1_DIRECT_RANK2 (2H_tau1-1)].
- M_plateau = E[1_LINE (H_tau1+H_tau2-1)].
- M_total = M_direct + M_plateau = F1+F2-1.

Both raw family contributions are reconstructed independently; neither is rescaled
by its own family count. Raw rank equality/inequality, exhaustive counts and hashes
are checked in the one analysis pass. Every paired delete-one M sum, its covariance
with unchanged K_A and the full total-wedge covariance reproduce immutable P439.

Within each generation the directions and all four coordinates share a deletion.
The four generations have distinct dependency groups (20/80/80/80 batches). Their
nonlinear wedge covariance contributions are summed, with no fictitious cross-N batch alignment.

## Scope and next scientific decision

These are measurement-only asymptotic chi-square/profile summaries, conditional on the
saved covariance estimates and retrospective choice of partition. They do not include
source/model uncertainty and are not prospective model-selection certification.
Opposite point-estimate signs are not by themselves resolved cancellation. Nor does a
surviving common ray establish nonzero coupling to K_A.

The root sensitivity repeats only this partition at the parent's already saved full/delete-one
pooled roots; roots, K_A, transfer order and source definitions are not reselected.

Direct and plateau point signs oppose at N85/N340 but agree at N170/N680. Their
correlations are only -0.127 to +0.057; there is no resolved large anticorrelated
pair at this precision. The root sensitivity changes the joint zero p-value
only from .790132 to .790122. The total-M result remains a missing-loading
observation, not positive evidence of shared marked/unmarked radial dynamics.

**Next output:** on an existing geometry with a separately resolved canonical M
contrast, produce M and the same natural K_A from one paired batch stream and
report their joint loading interval. First consume any compatible saved birth
archive; acquire only the missing same-stream rows if needed. Another K_A-only
N1360 point would not answer this question. Do not select a new M definition or
normalization just because the present loading is unresolved.

Move parent same-stream scoring and this direct/plateau split to **completed
analysis**. Reimplementing those scorers or repeating their synthetic checks is
not the next scientific output. Both derived families remain in the same four
dependency groups; no issue needs to be closed, locked or blocked.

## Reproduce

A clone must contain the immutable source objects listed in the manifest and score.
No server or external raw download is needed when those Git objects are present.

```sh
python3 scripts/analyze_p439_direct_plateau_transport.py
```

Output retains full and delete-one (K_A,M_direct,M_plateau,M_total) vectors, joint
covariances, all nine wedges, per-generation covariance contributions, raw provenance,
environment and one-pass reconstruction diagnostics. No parent analysis or test suite is rerun.

The numerical pass took 11.87 local seconds. The largest fixed-p point/delete-one
reconstruction discrepancy was 9.17e-16. The report interpretation was added after
that one pass; no raw data analysis was repeated to edit the narrative.
