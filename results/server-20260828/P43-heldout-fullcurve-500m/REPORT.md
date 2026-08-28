# Issue #43 prospective N=185/265 full-curve score

## Outcome

The preregistered two-sector model does **not** pass as a conjunction.  Its
matching-odd part survives the prospective test and decisively improves on
zero, but the matching-even effect is nonzero with the opposite sign from the
frozen prediction.

All scores below use the frozen model order and zero target-fit parameters.

| sector/model | N=185 observed | N=265 observed | frozen two-size chi-square |
|:---|---:|---:|---:|
| DeltaM, x=21/4 H4 | 1.36830e-4 +/- 3.42352e-5 | 1.27110e-4 +/- 3.46783e-5 | 3.04598 / 2 |
| DeltaM, zero | same | same | 29.40938 / 2 |
| DeltaM, x=17/4 adversary | same | same | 30.24613 / 2 |
| DeltaS, frozen positive N^-1 | -6.08154e-5 +/- 8.08957e-6 | -7.02495e-5 +/- 9.38562e-6 | 240.24721 / 2 |
| DeltaS, zero | same | same | 112.53891 / 2 |

The marginal DeltaM residuals from the x=21/4 prediction are -1.57 and
-0.84 combined standard errors.  The zero-effect marginal statistics are
+4.00 and +3.67 sampling standard errors.  In contrast, the DeltaS residuals
from the frozen positive prediction are -12.77 and -12.44 combined standard
errors; the observed negative effects themselves are -7.52 and -7.48
sampling standard errors from zero.

Thus the prospective result supports a nonzero odd square-harmonic sector
compatible with the frozen x=21/4 radial law, rejects the larger x=17/4
adversarial prediction, and falsifies the frozen sign/amplitude assignment for
the even sector.  It does not support describing the full result as a single
successful two-spin4 model.

The predeclared shared H4+H12 comparison is retained as `NOT_SCORABLE`: only
the exact harmonic columns were frozen, not a complete pre-target H4/H12
amplitude vector and source covariance.  No H12 amplitude is fitted to these
targets.  The invalidated wrong-Kac-branch V_<1,3> artifact is excluded.

## Issue #72 P4[S-prime] score from the same full curves

The derivative target is evaluated at each pair's intrinsic center and scored
only after the Issue #43 primary and secondary stages.

| frozen model | chi-square / df |
|:---|---:|
| original pure N^-5/4 law | 52.71634 / 2 |
| rank-2/Jordan log correction | 1.20360 / 2 |
| analytic inverse-N correction | 0.86221 / 2 |
| zero effect | 1278.55524 / 2 |

The new geometries prospectively confirm that the P4[S-prime] signal is
nonzero and that the original pure law fails.  Both frozen correction models
survive.  The analytic statistic is descriptively smaller, but the frozen
chronology is preserved: the rank-2/log correction is the primary correction
model and is not displaced by post-reveal selection.

## Production provenance

- Huawei environment: DevEnvC_ZyTrST, 16 vCPU / 32 GiB.
- Engine commit: `302464c3a08bdf74a8cea079a50cfebd7fc8843f`.
- Source SHA-256: `0273b1df8d0f91a4e76a7287a897bc82d4e52d003302dcfa33f797c95931ce2b`.
- Executable SHA-256: `f288da0db697efe24766e37f20d07af7a04b2e781e2324fc11731860011c4a17`.
- Compiler reported by production metadata: GCC 10.3.1 with OpenMP.
- N=185: pair (13,4)/(11,8), 500,000,000 replicas, counters [7000000000,7500000000), 100 batches, 8 threads, 1834.252 s.
- N=265: pair (16,3)/(12,11), 500,000,000 replicas, counters [7500000000,8000000000), 100 batches, 8 threads, 2531.543 s.
- Common seed: 2026104301; counter domains are disjoint.
- Both stderr files are empty.  Local downloads match the remote SHA-256 values.

The primary score reconstructs the microcanonical full curves at the frozen
`p_ref=0.592746050790`.  The derivative analysis separately solves each
orientation pair's intrinsic center.  Raw histograms, batch moments, metadata,
stdout, stderr, scorer outputs, derivative covariance, commands, and checksums
are retained.

## Evidence boundary

These are the first revealed N=185/265 target values and were scored against
immutable pre-target artifacts.  No exponent, amplitude, sign, harmonic
coefficient, correction coefficient, or thermal coordinate was fitted to the
target values.  Flexible matrix/Jordan interpretations remain downstream
models; they do not overwrite the failed frozen even-sector prediction.
