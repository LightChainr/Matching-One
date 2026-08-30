# P250 rank-eight R2 kernel-plane bridge freeze

## Question

The radius-five reveal selected one Alexander-R2-plus-conjugation map for a
single degree-two null line.  The radius-six reveal then eliminated every
endpoint-Hankel class through rank seven in both hands.  It did **not** reject
the full R2 bridge: the old bridge was locked behind a rank-five prerequisite
that failed.

This existing-data score asks the remaining minimal question:

```text
does the selected R2 map carry the complete two-dimensional rank-eight
right-kernel plane of H3_plus to that of H3_minus?
```

No alternate Alexander map or identity alignment is reopened.

## Frozen object

For each hand, reconstruct the same complex `20 x 10` matrix used by the
radius-six rank ladder,

```text
H3[(charge,u),v] = G_charge(u+v),
charge in {1,2}, |u|<=3, |v|<=3.
```

The minus entries are placed in the already-selected Alexander-R2 spatial
chart before `H3` is built.  Let `P_plus` and `P_minus` be the orthogonal
projectors onto the two smallest right-singular directions.  The frozen null
is

```text
conj(P_plus) - P_minus = 0.
```

This formulation is invariant under a basis change within either two-plane.
The decision statistic is the real/imaginary vectorization of the projector
difference, with delete-one covariance computed separately from:

1. the 80k radius-four block;
2. the 1.2M radius-five block;
3. the independent 1.2M radius-six block.

The three covariance contributions are summed.  The finite-batch Hotelling
calibration, relative covariance eigenvalue cutoff `1e-10`, and alpha `0.01`
are inherited unchanged from the locked radius-six scorer.

Principal angles and singular-value gaps are descriptive only.  They cannot
change the decision.

## Decision and boundary

```text
p < 0.01:
  rank8_R2_kernel_plane_bridge_rejected

p >= 0.01:
  rank8_R2_kernel_plane_bridge_compatible_at_this_truncation
```

A compatible result would extend R2 from one low-degree relation to the full
observed rank-eight relation plane.  It would not prove that rank eight is the
exact state dimension, that the moment table is flat at the next order, or
that eight physical fields exist.

A rejection would localize the radius-five R2 selection to its truncated null
line.  It would not establish noncommutation or context memory, because the
archive contains endpoint displacements rather than ordered `TxTy/TyTx`
histories.

This is one post-reveal structural diagnostic within the already declared
P250 dependency group, not a new independent evidence row.
