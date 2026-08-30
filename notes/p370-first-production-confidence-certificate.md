# Issue #370: first production confidence certificate

The first real archive connected to the proof-carrying envelope is the
eight-row `(A_top,E_top)` rank-plane crosswalk.  It contains P49 and P50 100M
production pairs and P43 and P57 500M production pairs, each with 100 aligned
batches and its full intrinsic-center `2x2` covariance.  It is not a synthetic
fixture, smoke stream or compiler control.

For a line model `E_top=r A_top`, distance to that line inside one Gaussian
outer confidence ellipse reduces to

```text
(E-rA)^2 <= K^2 (var(E)+r^2 var(A)-2r cov(A,E)).
```

The adapter freezes `K=13/4`, conservatively above the 99% familywise
Gaussian-Bonferroni critical value for eight rows.  Because every archived
decimal is parsed as a rational number, each feasibility margin and its sign
is then verified exactly; no optimizer, square root or fitted model parameter
enters the three exclusions.

The production confidence set eliminates all three frozen fixed lines:

```text
E_top=0, E_top=-A_top, E_top=A_top.
```

The free common-line class is not eliminated: the explicit rational witness
`r=-2/3` lies inside all eight simultaneous bands.  Thus model ordering changes
from three named cancellations to one surviving nuisance ray, without claiming
that this ray is the true or unique physical model.

The exactness is conditional on the declared Gaussian outer-confidence
contract.  It does not turn that asymptotic coverage statement into a
finite-sample theorem and does not use cross-archive independence or pooled
p-values.
