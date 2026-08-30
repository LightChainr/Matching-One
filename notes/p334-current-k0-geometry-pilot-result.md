# P334 result: cheap current geometry does not explain birth age

## Direct result

The fresh pilot observed the actual rank-one configuration at the frozen
intrinsic layers, rather than another birth-local or future mark.  It ran
50,000 paths per size, with two common-randomness orientations per size.

The cleanest exact result is the one-step carrier classification:

```text
all four size/orientation rows:
H2_figure8 = 0
H2_separate = 0
H2 = H2_theta
```

Every vacant site that can complete rank two in one insertion does so by
closing a nonparallel cycle inside an already essential component.  No joined
figure-eight trigger and no disconnected completing carrier occurred.  This is
a sharper finite-geometry statement than the age regression.

As required, H2 is treated as a calibration ceiling.  Its fitted hazard
coefficient is about one, and after H2 the N425 age slope retains only `9.2%`
and `7.1%` with joint `p=0.984`.  That is the identity
`P(next exit | current configuration)=H2/(N-k0)`, not evidence that H2 is a
new causal mechanism.

## Cheap geometry versus age

The frozen 50k primary is intentionally a pilot and its raw line-only slopes
are noisy.  The cheap vector (essential size/count plus occupied/vacant
frontiers) gives heterogeneous point attenuation:

| size | first retention | second retention |
|---|---:|---:|
| N325 | `22.6%` | `84.0%` |
| N425 | `90.7%` | `94.5%` |

More importantly, the frozen N325 geometry map transferred to held-out N425
leaves age slopes `-0.0720` and `-0.0650`.  The 50k joint test is underpowered
(`p=0.235`), but the effect sizes are essentially the production-scale age
signal rather than attenuation toward zero.

For a post-reveal precision description, the independent paired pilot shift
`beta_Mcheap-beta_M0` was added to the disjoint 2M production `M0` slope.  This
uses the large archive only as a precision anchor; it is not the frozen primary
score.  The resulting retained magnitudes are

```text
N325: 62.3%, 89.7%   joint p=1.33e-10
N425: 88.1%, 92.7%   joint p=1.25e-13.
```

Thus the compact current size/frontier/carrier vector does not absorb the
production birth-age association.  It removes a real part in N325-first, but
not a common mechanism across size and orientation.

## Crosswalk to the temporal modes

The retained K1/K2 rows permit a zero-new-simulation crosswalk to the external
two-time kernel result at `5a7f2d9`.  Applying its production eigenvectors to
the fresh paths, then regressing mode scores on the same line-centered cheap
geometry, gives:

```text
mode 2 R2: 1.61% -- 2.32%
mode 3 R2: 0.93% -- 1.39%.
```

So the stable nonleading temporal modes are not merely size/frontier modes.
This remains descriptive: the eigenvectors are convenient kernel coordinates,
not exact state identities, and the pilot shares the same observable family.

## Scientific card

- Mechanism space changed: cheap current mass/frontier/carrier summaries are
  removed as a common explanation of the production age slope; one-step
  completions are empirically all theta-like in these four rows.
- Not proved: intrinsic temporal memory, complete-state non-Markovianity, a
  scaling-limit memory field, or an exact three-state temporal theory.
- Observer/source/geometry: rank-one state at N325 `k0=193` and N425 `k0=252`,
  two paired norm-five quotient orientations, fresh counter streams.
- Dependency groups: one fresh block per size; the production anchor and mode
  coordinates come from existing disjoint 2M archives and are post-reveal
  descriptions, not extra independent votes.
- Next lift: record one compact shape/bottleneck coordinate that approximates
  H2 without enumerating every vacant site; repeating size/frontier counts has
  low expected information value.

## Provenance

The first launch on each Huawei node stopped before runner invocation because
`/usr/bin/time` was absent.  Both exit-127 logs are preserved.  K2 retention
was added transparently before any CSV existed.  The successful runs used the
same frozen seeds/counters, runner `12c79dc`, source SHA256 `d74b6dce...`, and
the frozen scorer was completed once after raw lock `424afdc`.
