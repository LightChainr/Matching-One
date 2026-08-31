# A coarse trigger mode with an Euler-invisible local response

The exact-census result `22952d75` and the fixed-degree contact result
`b9f79bfb` now share one saved twenty-batch covariance factor. This final join
reads their committed summaries only; it performs no new fork sampling,
contact enumeration, DP or refitting.

## What the refinement resolves

| N | binary between Gamma share, canonical | four-trigger between share, canonical | four-trigger between share, integrated |
|---|---:|---:|---:|
| 325 | 12.555 ± 5.873% | 93.943 ± 2.781% | 89.207 ± 6.517% |
| 425 | 5.583 ± 2.869% | 87.540 ± 3.360% | 76.992 ± 6.335% |

The increase in canonical between Gamma when retaining **which orientation
triggers** is .001124238 ± .000078218 and .000817410 ± .000043086.
This is the shared-covariance difference between two nested partitions of the
same next-label response. It shows where the large covariance previously
hidden inside the combined birth class goes.

For checkpoint group 01+10, four-trigger canonical between shares are
96.154 ± 5.657% and 86.322 ± 4.350%; the integrated shares are
90.556 ± 9.074% and 77.014 ± 8.861%. The binary between Gamma shares of that
group are negative, approximately −44% and −46%. Cross-covariance shares are
signed quantities; these signs do not contradict positive covariance matrices.

## What remains inside the coarse state

The contact readout fixes own-orientation R0, preservation of that rank, and
the next site's occupied-neighbor degree e. Its feature is loop=e−c, with c
the number of distinct occupied components touched. At fixed e, the change in
loop is minus the change in component mergers. The following are the original
contact result's pooled response slopes, carried into the common covariance:

| N | K1 slope | K2 slope | C slope | W slope |
|---|---:|---:|---:|---:|
| 325 | .59009 ± .04619 | .88512 ± .05995 | .73760 ± .04509 | .29502 ± .05763 |
| 425 | .71019 ± .06240 | .96502 ± .09445 | .83760 ± .06917 | .25483 ± .08056 |

This named contrast changes the expected future birth clock even though both
immediate rank and Euler increment Δχ=1−e stay fixed. In the source-defined
conditional domain it excludes a next-response closure depending only on e.
The equal-degree covariance also specifies the linear response to the
rank/Euler-distribution-preserving label tilt discussed in the parent result;
it is already in the saved raw contact span and needs no new covariance block.

These two observations fit a concrete finite-source picture: the coarse
paired response is largely organized by orientation-specific trigger type,
while contact topology resolves additional future-clock information on a
fixed-rank/fixed-Euler level set.

## Estimand and dependence boundary

The four-trigger partition concerns a **paired H4-normalized** observable with
both orientations' trigger labels. The contact result pools **single-orientation**
R0-safe labels, including an equal-degree restriction. Its signal does not
numerically equal the partition's remaining within Gamma. They are measured
on the same `e32a8593` source block and are not independent confirmations.
No field count, path-memory identity or complete state sufficiency follows.

All old covariance coordinates, the contact 87 raw batch columns, its 128
existing derived LOO columns, and the named partition ratios/differences are
retained as a rank-at-most-19 factor. This intentionally preserves redundant
linear relations, with no matrix inversion or clipping. The compact score
provides the focused joint covariance. No higher-dimensional omnibus is added.

Reproduce the thin join:

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_trigger_contact_joint.py
```

Outputs are under `results/p334-trigger-contact-joint/`. The underlying partition
theory, source hashes and one-pass scorer are documented in
`notes/p334-exact-birth-partition-response.md`.

Scientific card: the main mode is now assigned to a measured trigger-type
coordinate; a distinct local loop-versus-merger response survives equal rank
and equal Euler degree. Observer/source/geometry remain the original paired
N325/N425 fork experiment and its own-orientation contact projection. The
dependence group is unchanged; all covariance is coordinated once here.
