# P250 result: the radius-five R2 line does not extend to the rank-eight plane

## Result

The separately frozen existing-data test rejects the full rank-eight
Alexander-R2-plus-conjugation kernel-plane bridge:

| quantity | value |
|---|---:|
| frozen projector coordinates | 200 |
| algebraically zero Hermitian coordinates removed | 10 |
| identifiable input coordinates | 190 |
| resolved covariance modes | 99 |
| asymptotic chi-square | 7373.6811 |
| finite-batch Hotelling F | 56.1879 |
| finite-batch p | `4.2173e-149` |
| frozen alpha | 0.01 |

The decision is

```text
rank8_R2_kernel_plane_bridge_rejected.
```

All old4, radius-five and radius-six delete-one covariance contributions were
retained as one declared dependency group.  Identity and Alexander R0/R1/R3
were not reopened.

## What changed structurally

The two principal cosines between `conj(ker H3_plus)` and
`ker H3_minus_R2` are

```text
0.89870, 0.14513,
```

or principal angles approximately

```text
26.0 degrees, 81.7 degrees.
```

Thus one direction remains moderately aligned while the other is almost
orthogonal.  This gives a direct explanation of the chronology:

- radius five saw one degree-two null line and selected R2 at that truncation;
- radius six exposed a two-dimensional relation plane at the first compatible
  rank-eight class;
- the second relation does not obey the same R2 bridge.

The radius-five result was therefore not a false calculation.  It was an
accurate lower-order relation that does not extend to the full observed plane.

## Protocol chronology

The protocol was frozen and pushed at `aaf350e`.  The scorer and tests were
pushed before reveal.  The first invocation at `1923aa1` stopped before a
statistic or output because the inherited whitening routine rejected ten
algebraically zero projector coordinates, chiefly imaginary diagonal entries
of Hermitian projectors.  Commit `c0dc4f3` removes only coordinates with
exactly zero delete-one variance; the frozen projector difference, covariance,
eigenvalue cutoff, Hotelling rule and alpha are unchanged.  One completed
score was then produced from `c0dc4f3`.

## Boundary

This is a conditional, truncated relation-space result.

- Rank eight is the first class not eliminated by the upstream rank ladder;
  it is not proved exact or flat at the next order.
- The failure does not imply `Tx Ty != Ty Tx`, path memory or a physical
  eight-field spectrum.  Only endpoint displacements were stored.
- It rejects the selected R2 equality only.  It does not authorize a
  post-reveal vote among the four discarded Alexander maps.
- The radius-five single-line result remains valid in its original domain.

## Scientific card

- Mechanism space changed: a common plus/minus rank-eight quotient in the
  selected R2-conjugate gauge is removed.
- Not proved: exact state dimension, flatness, noncommutation or context
  memory.
- Observer/geometry: N505 C4-gauged projective-leg Z5 charged two-point
  endpoint moments, charges 1/2, degree-three H3.
- Dependency group: old4 80k plus radius5 1.2M plus radius6 1.2M.
- Next lift: record a signed, ordered two-morphism rectangle.  Another endpoint
  shell or Alexander-map revote has lower information value.
