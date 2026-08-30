# P334: aggregate active-boundary organization is a no-op

## Frozen result

The fresh 20k/size pilot recorded five one-pass organization summaries:
cut-edge axis anisotropy, corner-versus-opposite contact balance, number of
vacant-frontier components, largest frontier arc, and component concentration.
No configuration or component labels were saved.

Relative to the opposite-size rank-one morphology state from `8fd596d`, the
production-anchored residual age slopes change as follows:

| target | rank-one residual | plus organization | retention |
|---|---:|---:|---:|
| N325 first | `-0.02907` | `-0.03152` | `108.4%` |
| N325 second | `-0.04318` | `-0.04001` | `92.7%` |
| N425 first | `-0.03129` | `-0.03132` | `100.1%` |
| N425 second | `-0.03110` | `-0.03155` | `101.5%` |

The frozen absorption rule fails decisively on effect size: it required every
retention to be at most `25%`.  These coordinates leave `92.7--108.4%` of the
remaining signal.  Larger samples of the same global aggregates have low
expected information value.

## Temporal bridge and the actual positive result

The five new fields add only `0.025--0.129` percentage points of line-centered
variance to temporal modes 2 and 3.  The largest increment is `1.79` delete-one
standard errors.  Thus they are also nearly a no-op for the temporal kernel.

The baseline used for that test reveals a more useful structural bridge.  The
rank-one morphology direction was learned from the *opposite size* next-exit
hazard and applied to fresh target paths.  After fitting only its target
amplitude, it explains:

| target | mode 2 R2 | mode 3 R2 |
|---|---:|---:|
| N325 first | `54.40% +/- 0.68%` | `11.66% +/- 0.54%` |
| N325 second | `55.19% +/- 0.70%` | `11.45% +/- 0.70%` |
| N425 first | `43.52% +/- 0.89%` | `23.29% +/- 0.90%` |
| N425 second | `43.51% +/- 0.63%` | `23.72% +/- 0.58%` |

So the transferable morphology direction is not merely a weak hazard
regressor: it is strongly aligned with temporal mode 2 and moderately with
mode 3.  The missing age information is not recovered by coarse organization
totals around that state.

This crosswalk is post-reveal.  The subspace is held out across size, but its
target amplitude is fitted; it is a state-space alignment, not a point
forecast or an independent evidence block.

## Scientific card

- Mechanism removed: global active-boundary anisotropy, component count,
  largest arc, and concentration do not complete the cross-size rank-one state.
- Positive bridge: the existing opposite-size morphology state carries about
  `44--55%` of temporal mode 2 and `11--24%` of temporal mode 3 variance.
- Not proved: intrinsic memory, universality of the latent direction, or
  absence of a local marked-boundary sufficient statistic.
- Next lift: if current geometry is pursued again, record a localized boundary
  pattern relative to `ell` or the completing sites, not another global total.

## Provenance

Freeze/runner/scorer commit: `3f190ed`.  XP ran both fresh seed/counter blocks.
The first N325 shell command completed before its wrapper could persist a PID
or exit marker; its completion stdout, exact metadata, and empty stderr are
preserved.  N425 has exit code zero and a completion timestamp.  Raw was locked
at `9315f3e` before any organization distribution or score was inspected, and
the frozen score was executed once.
