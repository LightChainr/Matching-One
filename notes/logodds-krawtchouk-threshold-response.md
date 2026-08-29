# Log-odds / Krawtchouk threshold response

Issue #182 proposes reorganizing the threshold-rank sufficient statistics in
the Bernoulli natural coordinate

```text
eta = log(p/(1-p)).
```

The first executable step is now `scripts/threshold_score_modes.py`.  It
constructs the orthonormal binomial Krawtchouk modes at the intrinsic center,
projects the microcanonical primal and matching-partner curves mode by mode,
forms the existing S/D and P4 combinations, and repeats the complete operation
inside every delete-one replicate.

## A useful no-go result for the S-prime anomaly

With the sign convention

```text
H1(K) = (K-Np)/sqrt(Np(1-p)),
c1 = E[O H1],
```

the exponential-family score identity is exact:

```text
dR/deta = sqrt(Np(1-p)) c1,
dR/dp   = sqrt(N/[p(1-p)]) c1.
```

Therefore the P4[S-prime] observable and the first Krawtchouk S mode differ
only by the known scalar `sqrt(N/[p0(1-p0)])`, evaluated at the same intrinsic
center.  Changing from bare `p` to `eta` cannot, by itself, cure the observed
first-derivative scaling failure.  Any simplification must come from one of:

1. center transport/nonlinear intrinsic-coordinate effects;
2. coupling to score modes of order two or higher;
3. an ordinary additional eigenfield;
4. a genuine triangular/Jordan transfer block.

This narrows #182 before fitting anything.  Mode 0 and mode 1 must be treated
as exact alternate views of the already-scored value and derivative, not as
new evidence.

## Two interlaced score-mode towers

The natural basis also exposes a parameter-free exponent ladder.  A thermal
response of order `r` supplies `N^(3r/8)` through the near-critical variable,
whereas the orthonormal binomial score removes `N^(r/2)`.  Relative to the
center amplitude, score order `r` therefore advances the decay exponent by
exactly `r/8`.

Matching parity selects two interlaced towers:

```text
matching-even I family, base alpha=1:
  P4[S_0] N^1, P4[D_1] N^(9/8), P4[S_2] N^(5/4), ...

matching-odd T family, base alpha=13/8:
  P4[D_0] N^(13/8), P4[S_1] N^(7/4), P4[D_2] N^(15/8), ...
```

This turns the four familiar value/derivative laws into the first two rungs of
two full orthogonal towers.  Modes `r>=2` are the new predictions.  Coherent
collapse across independent sizes would support a common scaling function;
one drifting rung, especially `T:S_1`, would localize the correction inside a
family rather than motivate unrelated scalar fits.

## Source-data result: a coherent thermal tower through order six

The basis and exponent ladder above were fixed without using the norm-5
target.  Applied to the existing doubled-family source blocks, they reveal a
high-signal matching-odd tower rather than isolated derivative behavior.  The
three even `D` rungs have the following scaled amplitudes:

| score order | scaling | N=130 | N=170 | N=185 | N=265 |
| ---: | :--- | ---: | ---: | ---: | ---: |
| 2 | `N^(15/8) P4[D_2]` | 1.327(102) | 1.392(115) | 1.666(82) | 2.026(109) |
| 4 | `N^(17/8) P4[D_4]` | -2.268(194) | -2.459(249) | -2.700(179) | -3.544(246) |
| 6 | `N^(19/8) P4[D_6]` | 2.833(437) | 3.269(669) | 3.338(473) | 5.497(721) |

The interlaced odd `S` rungs are negative at orders three and five.  Across
N=130 and N=170, all five new thermal rungs `D_2,S_3,D_4,S_5,D_6` keep their
sign and differ by less than `0.7` combined standard errors after the frozen
scaling.  Across N=185 and N=265, the even-`D` subvector has cosine similarity
`0.9924`; its magnitude drifts, but its alternating `+,-,+` direction does
not.  By contrast, the matching-even high-order tower is weak and unstable in
the same data.  The present source evidence therefore supports one coherent
thermal/interlaced response spectrum through order six, not two equally
resolved towers.

The order-one `S` mode is positive at about 25 standard errors in both N=185
and N=265.  Its reconstruction of the directly computed `P4[S-prime]` closes
at roughly `1e-40`, while the naively scaled amplitude still drifts from
`1.267(49)` to `1.404(57)`.  Thus the S-prime anomaly is a real matching-odd
channel, not a finite-difference artifact.  The simultaneous drift of
`D_2,D_4,D_6` suggests family-wide correction or triangular/logarithmic
mixing; it is not well described as an unrelated scalar slope failure.

These are descriptive source results.  Modes zero and one remain exact
coordinate views of already scored observables, and the norm-5 target remains
unseen when defining the basis, order, sign, and exponent ladder.

## Next score

Use the now-frozen modes `r=0..6` to predict the held-out norm-5 lineage.  The
primary target is the signed thermal vector `D_2,S_3,D_4,S_5,D_6`, with the
even-`D` direction as the high-signal subtarget.  Do not rotate it using the
target.  Independent Gaussian-multiplier stability can distinguish a shared
response direction from one noisy high-order coefficient.
