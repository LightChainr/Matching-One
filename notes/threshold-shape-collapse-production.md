# Production threshold-shape flow from existing rank archives

Status: retrospective production analysis under the exact seven-quantile contract merged in `8e5a236`.

The finite threshold distribution is reconstructed as the equal mixture of the
`Kminus/Kplus` beta order-statistic laws and both Gaussian orientations. Median
and IQR are removed inside every delete-one batch. Since the standardized q25,
q50 and q75 entries are structural, the covariance test uses q05, q10, q90 and
q95 only.

The four parent-to-cover comparisons reject exact finite-size shape equality at
the available precision. That is not the main discovery. The norm-2 and norm-5
shape changes from the same parent are almost collinear: the Euclidean cosines
are approximately `0.999995` for N65 and `0.994739` for N85. Thus cover size
mainly advances one common tail-deformation coordinate rather than producing a
new arbitrary shape for each multiplier.

The deformation is strongly one-sided after location/IQR removal: the largest
component is the upper q95 tail. The norm-5/norm-2 projection scales are about
`1.8000` and `1.9386`. Interpreted as

```text
Delta_Q shape = amplitude * (1-Q^-q) * v_shape,
```

they give effective powers `0.632` and `0.429`. The first is strikingly close
to `5/8` (which predicts a cover ratio `1.8041`), while the second shows that a
single universal amplitude/power is not yet identified. This is a nominated
shape-flow direction, not a fitted CFT field.

The useful next test is no longer “does the full curve collapse exactly at
finite N?” It is to freeze this tail vector on the N65/N85 families and test
whether a fresh lineage changes only its amplitude, then connect that amplitude
to the existing Hermite--Krawtchouk thermal jet.

## Independent archive replication

After commit `2b9b91b` had recorded the cover-family direction, the independent
P43 N185/N265 500M blocks were inspected with the same scorer. Their residual
vector is

```text
[0.000129, 0.000880, 0.000676, 0, 0.000676, 0.003951, 0.007753].
```

Its cosine with the four cover-flow vectors ranges from `0.99155` to
`0.99962`. This is not a prospectively generated target, but it is an
independent counter block and a different size ratio. The rank-one direction
therefore replicates beyond the two cover pairs and can now be frozen as a
genuine next-lineage prediction.
