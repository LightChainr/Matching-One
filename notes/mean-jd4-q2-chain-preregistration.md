# Preregistered three-generation mean-`J_D4` q2 chain

This experiment follows the exact contact no-go: the primary observable is
only the **mean** lifted-line odd source, not `Cov(A_top,J_D4)`.

The frozen chain is

```text
N65  : (8+i, 7+4i),       Delta cos4 = +1152/845
N130 : (9+7i, 11+3i),     Delta cos4 = -1152/845
N260 : (16+2i, 14+8i),    Delta cos4 = +1152/845.
```

Each size uses 5,000,000 fresh counter replicas, 100 batches, seed
`202608290315`, offset `9300000000`. N260 is unobserved when this artifact is
committed.

The primary complex targets are

```text
J130/J65  = -2^(-13/8),
J260/J130 = -2^(-13/8),
J260/J65  = +4^(-13/8).
```

Phase is scored before magnitude. A real exponent is fit jointly to both
edges, but the two complex transfers are also reported freely. Intrinsic
centers are recomputed per size and inside every delete-one score.

Existing timing implies roughly 82, 157, and 315 seconds on the three 16-core
Huawei hosts under a conservative 75% parallel-efficiency assumption. The
largest raw stream should remain below 1.6 GiB. Extrapolated N260 combined SNR
is about four, so this is a discovery experiment rather than a precision
exponent measurement.

The N65/N130/N260 counter IDs are aligned for provenance. Only the two
orientation geometries within each size share an exact permutation field;
cross-size samples are not mislabelled as configuration-identical.
