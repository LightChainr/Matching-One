# Covariance-nullspace scorer audit

Issue #543 identified a real fail-open edge case in three generalized
chi-square implementations. For an exact covariance

```text
Sigma = [[1, 1], [1, 1]],   r = [1, -1],
```

the old scorers discarded the only direction containing the residual and
reported `chi-square=0`, rank one and survival probability one. The shared
kernel in `scripts/covariance_nullspace.py` now reports every discarded
projection, freezes five cutoff-sensitivity rows, rejects materially indefinite
covariance, and fails closed when a caller declares a structural null
constraint. Estimated covariances retain the active-subspace statistic but mark
an incompatible discarded projection as requiring sensitivity review.

## Historical rescore

The audit uses only the residual vectors and covariance matrices already stored
in six archived score files. It covers all 16 stochastic vector scores reached
by the three implementations, including the PR #273 norm-4 production result
and PR #277 generation-4 pilot that are not present in this branch. No raw
histogram was reread, no random sample was generated, and no evidence count was
added.

| archived group | vector scores | frozen-cutoff result |
|---|---:|---|
| P50 N145 to N290 full curve | 1 | displayed statistic unchanged; interpretation changed |
| P57 norm-5 thermal jet | 3 | unchanged |
| P57 conjugation parity | 2 | unchanged |
| P57 intrinsic functional cocycle | 2 | unchanged at the frozen cutoff; broad-cutoff sensitivity recorded |
| P154 norm-4 production | 2 | unchanged |
| P154 generation-4 pilot | 6 | unchanged |

Fifteen scores have no discarded mode at their historical `1e-10` cutoff. The
only default-rank reduction is the P50 three-level `DeltaM` transfer:

```text
discarded correlation eigenvalue       1.2378454349e-13
relative eigenvalue                     4.1261538746e-14
standardized residual projection       -1.7749253447e-6
discarded chi-square contribution       25.4503501849
```

Its frozen result remains `chi-square=9.3520036848/2`, survival
`0.0093161871`; the original single-multiplier rejection does not reverse.
The sensitivity rows are decisive for wording:

| relative cutoff | rank | chi-square / df | survival |
|---:|---:|---:|---:|
| `1e-14` | 3 | 34.80235387 / 3 | 1.34119e-7 |
| `1e-10` | 2 | 9.352003685 / 2 | 0.00931619 |
| `1e-6` | 1 | 1.924293999 / 1 | 0.165384 |

The mode therefore cannot be called a harmless numerical null. The archived
default statistic is preserved, but the claim is now explicitly
estimated-near-null and cutoff-sensitive. The two functional-cocycle scores
remain full rank at `1e-10`; their original raw-covariance eigenbasis is
preserved rather than silently replaced with a correlation eigenbasis.

## Decision

This completes the bounded QA task. Existing displayed statistics do not need
replacement, one interpretation needs correction, and no production rescore or
new Monte Carlo is justified. Future structural covariance contracts must opt
into fail-closed mode; estimated covariance reports must carry the nullspace
status and cutoff table emitted by the shared kernel.

Run `python3 verify.py` in this directory to reproduce `RESULT.json` and check
the exact regression conditions.
