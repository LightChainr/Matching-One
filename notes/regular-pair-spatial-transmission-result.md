# A regular local interaction has noncontact Q-activated transmission

**The frozen finite-distance zero-transmission null is rejected.** The
canonical regular interaction has a positive occupation-averaged connected
Q response at separation 16 on the L64 square torus. This comes from two
fresh random blocks, not another projection of the old N25 population.
The fixed production question is complete; no samples or completion
coefficients are added after this readout.

The same work gives a complementary exact result: an entire class of
regular homogeneous local deformations cannot change the original
moving-root U directly at Q=1. Thus direct Q1 transmission and
Q-activated spatial transmission are now separated constructively.

## 1. The measured interaction and its exact spatial mechanism

Keep the canonical completion

```
Kreg = average_C4 i(I-P1)i^dagger,     c(Q) = 1 identically.
C(x,y) = partial_logQ partial_lambda_x partial_lambda_y log Z | Q1,lambda=0.
```

The lambda parameters multiply the vacant-site vertex at x and y and are
not divided by the number of sites. Both points must be vacant for a
nonzero contribution; all configurations remain in the denominator.
Original NN site occupations are sampled at the fixed reference
`p_ref=0.592746050790`, implemented as
`10934234699625173385 / 18446744073709551616`. This is a prescribed
sampling parameter, not a new certification of the critical threshold.

The [observer derivation](regular-pair-spatial-observer.md) proves

```
C(x,y) = E_iid[g_pi(A)],
g_pi = partial_Q [Q^(-|pi|) sum_colours Kreg_x Kreg_y D_pi] | Q1.
```

Here pi is the exterior-component partition of the eight incident edge
ports. The background-measure derivative and disconnected-product
derivative vanish because every nonempty regular insertion is zero at
Q1. Consequently this is the actual connected colour contraction, not
`Cov(a_x,a_y)`, not a covariance of two separately closed marks, and not
the original U.

The [exact kernel](regular-pair-spatial-kernel.md) gives a sharp microscopic
constraint: **zero or one shared exterior component implies g_pi=0**.
With exactly two shared components,

```
g_pi = d_x(1) d_y(1),
d_x = F_x(equal shared colours) - F_x(unequal shared colours).
```

Thus this observer couples colour contrasts carried by two components
reaching both marked neighbourhoods. It is not an unsigned bridge count.
For example the physically realizable port patterns `00110011` and
`00110012` have weights `+1/16` and `-1/8`. The Bell8 calculation contains
1874 nonzero partitions, including 480 negative ones; these are algebraic
partition counts, not probabilities or claims of universal realizability.

## 2. Two fresh blocks: the completed finite-distance decision

| Square torus L | Distance r | C | Monte Carlo SE | 99% interval |
|---:|---:|---:|---:|---:|
| 32 | 8 | 3.6591796875e-5 | 1.3775571881e-6 | [3.3009100945e-5, 4.0174492805e-5] |
| 64 | 16 | 6.8554687500e-6 | 6.3522637112e-7 | [5.2033972758e-6, 8.5075402242e-6] |

Each size used 200000 fresh configurations in 200 batches. The 16 fixed
anchors and two fixed directions give 32 correlated pair contributions
per configuration, averaged before inference. The sample size is not
6.4 million independent pairs per geometry. The two geometries use
separate frozen RNG streams and share no production configurations.

The primary two-sided Student-t readout on the 200 L64 batch means gives
`p=1.1328120879e-21`, rejecting `C64=0` under the declared 99% Monte Carlo
decision rule. This is a Monte Carlo inference, not an exact probability
bound. The secondary simultaneous size/distance ratio is

```
C64(r16)/C32(r8) = 0.1873498799,
99% Fieller interval = [0.1401616411, 0.2381648820].
```

The full shared-component vector and its joint batch covariance are in
[score.json](../results/p337-regular-pair-spatial/score/score.json).
The interval describes this fixed ratio; no exponent has been fitted.

## 3. What carries the observed signal

All observed nonzero contributions in both blocks have exactly two
shared exterior components. The predeclared decomposition is:

| Count | L32/r8 | L64/r16 |
|---|---:|---:|
| Total prescribed pair readouts | 6400000 | 6400000 |
| Both endpoints vacant | 1060448 | 1062314 |
| Exactly two shared components | 2436 | 445 |
| Nonzero signed kernel | 1060 | 209 |
| Sum of g16 = 16g | 3747 | 702 |
| Observed pairs with three or four shared components | 0 | 0 |

No observed three-/four-component event does **not** prove such events
impossible or their population contribution zero. The zero/one-component
zeros, by contrast, follow from the exact factorization theorem.
Negative kernel weights were retained. The collected sums do not resolve
positive and negative masses separately, so they do not establish absence
of cancellations.

One descriptive accounting identity makes the size change transparent.
On the observed blocks, since no higher-shared-count contribution occurs,

```
C_observed = frequency(two shared components) * mean(g | two shared).
frequency: 0.000380625 -> 0.00006953125       (ratio 0.1826765189)
conditional mean: 0.0961360837 -> 0.0985955056 (ratio 1.0255827134)
product ratio: 0.1873498799.
```

This is post-readout arithmetic on the predeclared counts, not a fitted
descriptor, another independent result, or a test of constant conditional
amplitude. Descriptively, the observed decrease is chiefly a decrease in
the frequency of two-shared-component configurations.

## 4. A broader exact exclusion for direct Q1 transmission

The [local scalar-reduction theorem](q1-regular-local-scalar-reduction.md)
does not assume the tensor vanishes at the all-one colour assignment.
For any finite homogeneous site-local deformation using only the original
edge colours, regular at Q1 and retaining the original vacant/occupied
summands, singleton-colour specialization leaves two scalar weights w0,w1:

```
P(A) = p_eff^K (1-p_eff)^(N-K),     p_eff=w1/(w0+w1).
```

If the weights are the same on the paired geometries and the moving
root and thermal Jacobian are regular, then

```
U(epsilon) = A_N Y0'(p0)/M0'(p0),
```

independent of epsilon to every order. A regular completion within this
class cannot preserve the old nonzero direct pure-pair Q1 tangent.
This rules out a microscopic model class rather than merely finding a
small coefficient. Surviving alternatives must change an assumption,
such as retaining additional states, introducing genuine multisite
occupation interactions, using nonuniform/geometry-dependent weights,
or changing the endpoint/limit prescription. It does not rule out all
Q1 models, nor the Q derivative measured in this spatial experiment.

## 5. Lifecycle, dependency and stop boundary

The [contract](../analysis/p337_regular_pair_spatial_contract.json) was
frozen at `3210aeb338ca7bb52c799d1de9048232f50ab921`, and the scorer at
`0096e79469fc9b8de00ebdafe52226345c65c364`, before data collection.
The exact kernel `32ff99fa5361ba0fe435fac835be2dbb206e0a6c` is public as
`c29d8bce`; the scalar theorem `7f60e92d5cdb58e7542db06cd49547a4451ba022`
as `5d24f99b`; the final producer `9f6ff44dd41764fc34a251b202494172e62228b6`
as `3477dfe0`. Raw production `00bfeda6d7afccbde81106decf6f693b741bed81`
is public as `6425bffbe56c4eca42d74410dcf45dba3899e74b`.
The deterministic accounting-bound addition `49f20d0c` preceded readout
and changed neither estimator nor primary test.

The [raw receipt](../results/p337-regular-pair-spatial/run.json) records
commands and file hashes. L32 took 2.134110625 seconds and L64
7.974634458 seconds on two local workers; the single numerical score
took 0.04900025 seconds. There was one compilation, no production
restart, no top-up, no old-source rescore and no cloud startup.

**Science card:** observer = connected two-insertion Q susceptibility;
sector = unprojected canonical C4-averaged colour interaction;
source = two vacant-site Kreg coefficients; geometry = square L32/L64
periodic tori at r=L/4. The two fresh size blocks are separate dependency
groups; every within-size decomposition reuses its same configurations.

This excludes the specified finite contact-only explanation. It neither
identifies a continuum H4/Jordan field nor establishes that this spatial
signal generates the original Matching-One anomaly. The finite question
is closed. A next field claim needs a named projection and a prediction
that separates candidate fields; it must not recycle these two points
into a fitted exponent or a new completion choice. P154/P334/F4 production
stops and other teams' task priorities are unchanged.
