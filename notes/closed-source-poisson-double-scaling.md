# Poisson rank coexistence closes the oblique pooled-U double-scaling limit

**New result.** In the joint limit `N/m² -> zeta < infinity`, the same
closed source has a finite Poisson particle cloud in its rank-zero phase,
and a pure full configuration in its rank-two phase. This holds on growing
honest square-cell quotients, including the oblique Gaussian companion.
At the pooled matching root the two phases have equal limiting weight,
the within-geometry denominator remains macroscopic, and original U
vanishes faster than any inverse power of N.

The statement covers both already fixed laws, `Sstar` and `Sdrop=Sstar+r`.
Their leading phase distributions agree at their respective roots, but
their logit roots differ by `2 log(m)/N+o(1/N)`. This is a controlled joint
limit, not an interchange of the old fixed-N and fixed-t limits.

No new configurations, couplings, fitted exponent or descriptor enter the
proof. It starts from `e17b286b` and the same occupied-NN / vacant-matching
corner convention as [the winding barrier](closed-source-winding-barrier.md).

## 1. Two fixed laws and the theorem's scope

On `G=Z²/Lambda`, let N be its number of sites, K the occupied count,
Bmix the mixed NN edge count, C_B the number of occupied NN components,
and r their ambient image rank. Put `g=Bmix-2C_B+r`, `m=exp(t)`, and
`h=p/[(1-p)m]`. After removing the common constant, the two weights are

```text
w_delta(A)=h^K m^(-g+delta r),
delta=0: Sstar;       delta=1: Sdrop.                              (1)
```

Delta is not a tunable source in this note; it indexes these two specified
laws only. Below we write delta in lower case to distinguish it from the
angular normalization `Delta_angle`.

Consider a fixed finite number of equal-area quotient sequences with
`N -> infinity` and shortest nonzero Manhattan period `ell_N -> infinity`.
The local square stencil is injective; eventually every fixed-size patch
embeds in the covering square lattice. Assume the NN dual graph has no
cycle shorter than4, as holds eventually for these sequences. In particular
the usual `(5k,0)` / `(4k,3k)` pair is included. No reflection symmetry,
checkerboard parity, or integer restriction on m is needed.

Let `m_N>3`, `N/m_N² -> zeta in [0,infinity)`, and use the common chart

```text
h_delta(s,N)=exp[(s-2delta log m_N)/N].                           (2)
```

Uniformly for s in compact real intervals, in each geometry separately,

```text
Z_delta -> exp(zeta)+exp(s),
P0 -> exp(zeta)/(exp(zeta)+exp(s)),
P2 -> exp(s)/(exp(zeta)+exp(s)),       P1 -> 0.                   (3)
```

Here Z_delta is the partition function of (1), with the empty configuration
having weight1. Conditional on rank0, K converges in total variation to
`Poisson(zeta)`. Conditional on rank2, the probability of the single full
configuration tends to1. The rest of the note proves this statement and
the original pooled-root consequence, rather than replacing U by a proxy.

## 2. A contour bound that works on oblique tori

Use the resolved mixed-edge dual curves from the winding-barrier proof.
At alternating corners the arcs round occupied corners, preserving NN
black connectivity and matching white connectivity. The curves are embedded
after smoothing, use every cut edge exactly once, and have total integer
length Bmix. Before smoothing a walk may revisit a dual vertex but does
not immediately backtrack and never repeats a cut edge.

For a proper nonempty occupation set every black component has a nonempty
boundary; each boundary curve belongs to one black component. Thus, if
Gamma is its family of resolved curves,

```text
C_B <= |Gamma|,
2C_B-(1-delta)r <= 2|Gamma|,           delta=0 or1.                (4)
```

The full configuration is the sole exception to the first inequality;
it has no contours and will never belong to a contour event below.
Given a curve family, its cut edges determine the occupation on the
connected torus graph up to a global complement, hence at most two states.
Not every family or complement is admissible; retaining them only enlarges
an upper bound.

There are at most `4N 3^(n-1)` rooted, oriented, nonbacktracking dual walks
of length n. This also bounds the number of possible curves of that length.
Assign each curve the positive activity `u_gamma=m^(2-|gamma|)`.
Relaxing disjointness, nesting and cut consistency to arbitrary subsets
of the finite set of possible curves gives the product bound

```text
sum_gamma u_gamma <= A_N(m)=108N/[m²(1-3/m)],
sum_(|gamma|>=R) u_gamma <= B_R(N,m)
                         =4Nm²(3/m)^R/[3(1-3/m)],     R>=4.      (5)
```

For either law, `h^K<=M_h=max(1,h^N)`. Its partition function is at least
`1+m^(2delta)h^N`, from the empty and full states, and so is at least M_h.
Equations (4)-(5) therefore give the completely finite, activity-uniform
bound

```text
P_delta(any resolved contour of length >=R)
 <= min{1, 2 exp[A_N(m)] B_R(N,m)}.                              (6)
```

Indeed the unnormalized bad sum is at most
`2 M_h product_gamma(1+u_gamma) sum_(|gamma|>=R)u_gamma`; divide by
Z_delta and use `product(1+u)<=exp(sum u)`. The extra full-phase factor
of Sdrop causes no loss: the bad event has nonempty contours and (4)
already bounds its cluster reward. The denominator retains that full
factor. In particular no rank fugacity has been silently discarded.

Rank1 requires an essential boundary curve of length at least ell_N.
Consequently both laws obey, at every h>0,

```text
P1 <= epsilon_N(m)
    = min{1, 8Nm² exp[A_N(m)] (3/m)^ell_N/[3(1-3/m)]}.           (7)
```

For fixed m, the volume factor exp[A_N(m)] prevents (7) from proving
thermodynamic suppression. It is not a replacement for the previous
fixed-m axis chessboard theorem. But if N/m² is bounded, A_N(m) is bounded
and (7) decays faster than every inverse power of N as ell_N tends to
infinity. This is precisely the regime needed here. Unlike that earlier
theorem, (6)-(7) apply directly to oblique quotients and real m>3.

## 3. Short contours leave only two backgrounds with bounded defects

Fix R=16 and take ell_N>R. Equation (6), with N/m² bounded, gives

```text
P_delta(bad_R)=O(N m^-14)=O(N^-6),                               (8)
```

uniformly in h; constants depend on the bound on N/m², not on the
geometry. All remaining contours are contractible and have length<R.
Their lifted disks contain a bounded number `C_R` of lattice sites,
independent of N. For example `(R+2)²` is a safe bound from the lifted
bounding box. The disjoint smoothed curves have nested or disjoint
contractible interiors. Outside the outermost disks there is one connected
genus-one region of a single occupation colour.

If this exterior is vacant, every black NN component is contained in
one bounded disk and r=0. If it is occupied, the exterior supplies the
two torus directions and r=2. Thus a good configuration is a dilute
black-component gas or a full background with bounded defect patches.
There is no third macroscopic background hidden in the short-contour
restriction. Nested islands are retained in the patch weights below.

### Rank-zero gas

A black NN component C of size k and mixed perimeter b has exact activity
`h^k m^(2-b)`. Distinct components must be nonadjacent. The single-site
activity is `x=h/m²`. Every other bounded component has perimeter at
least6 and activity at most a constant times `m^-4`: h is bounded above
in (2) on compact s intervals, including the projection-deleted chart.
There are at most a fixed C_R-dependent constant times N such embedded
component shapes.

The subensemble of only single-site components is the NN independent-set
gas. Its partition function Z_ind satisfies

```text
(1+x)^N [1-2N(x/(1+x))²] <= Z_ind <= (1+x)^N.                   (9)
```

This is simply the iid Bernoulli measure of no adjacent occupied pair,
bounded by the union bound over 2N edges. Every such configuration is good
under the occupied-corner resolution. Adding all nonsingleton polymers
and dropping their compatibility constraints supplies the upper bound

```text
Z_ind <= Z0_good <= (1+x)^N exp[O(Nm^-4)].                       (10)
```

In chart (2), `Nx -> zeta`. If zeta>0, log(m)/N tends to0, so h tends
to1. If zeta=0, `Nx<=exp(s/N)N/m² ->0`, even when m grows very fast.
Since `Nx² ->0` and `Nm^-4 ->0`, equations (9)-(10) give

```text
Z0_good -> exp(zeta).                                          (11)
```

The nonsingleton mass vanishes. The independent-set gas is within
O(Nx²) total variation of iid Bernoulli occupation; its occupied count
therefore converges to Poisson(zeta). Bounded component sizes also control
the differentiated polymer bound, so the conditional first moments
converge. This is an actual dilute-particle limit, not a fitted Poisson
approximation to an old histogram.

### Rank-two background

Divide a good rank-two configuration by the full-state weight
`m^(2delta)h^N`. Its outermost vacant disks are disjoint bounded patches,
possibly containing finite black islands. For a patch with v vacant sites,
i internal black NN components and total mixed perimeter b, the exact
relative activity is

```text
h^-v m^(-b+2i).
```

The outer boundary costs at least4 edges and each internal black component
adds at least4 of its own boundary edges. Therefore `b>=4+4i` and
`b-2i>=4`. There are only finitely many local patch patterns per site.
This explicitly retains the cluster reward of nested islands.

In the Sdrop chart, the potentially small h is harmless for a bounded
patch: `h^-v<=exp(|s|C_R/N)m^(2C_R/N)`. For large N the total activity
of all patches is bounded by `O(Nm^-3)=o(1)`; for Sstar the sharper
`O(Nm^-4)` suffices. Dropping compatibility again bounds the relative
partition sum by an exponential of this total activity. Hence

```text
Z2_good/[m^(2delta)h^N] ->1,
Z2_good ->exp(s).                                               (12)
```

The conditional probability of any defect patch tends to zero, and its
mean vacant count does too. Equations (8), (11) and (12) prove (3), since
`Z=Z_good/(1-P(bad_R))`. The normalized bad mass also remains negligible
for the first two occupation moments by (8).

## 4. The pooled root and the denominator now close together

At chart value s, every geometry has the same limiting mean

```text
q_g -> tanh[(s-zeta)/2].                                       (13)
```

Both fixed laws have unique simple homogeneous roots, as proved earlier.
Evaluate (13) at s=zeta-epsilon and s=zeta+epsilon; their signs bracket
the equal-weight pooled root. Uniform convergence on this interval and
monotonicity imply, for each of the two laws,

```text
N log h_delta,root +2delta log m -> zeta.                        (14)
```

Thus individual geometric q means tend to0 at the **same pooled root**,
not merely to equal and opposite nonzero values. At this root each has

```text
P0 ->1/2, P2 ->1/2, P1 ->0,
K | r=0 -> Poisson(zeta),       P(K=N | r=2) ->1,
Cov_g(q,K)/N ->1/2.                                            (15)
```

Consequently the average within-geometry variance kappa tends to1, and
the actual pooled thermal denominator satisfies

```text
D_h=mean_g Cov_g(q,K)/h ~ N/(2h).                               (16)
```

This resolves the denominator distinction in the controlled joint limit.
It never substitutes the variance of a mixture of geometries for the
within-geometry covariance.

There is also a concrete microscopic correction to the naive two-state
picture. The probability of the single empty configuration tends to
`exp(-zeta)/2`, while that of the single full configuration tends to1/2.
The remaining mass tends to `(1-exp(-zeta))/2`, entirely on the dilute
rank-zero cloud. At a single geometry's own root one further obtains
`E[K]=N/2+zeta/2+o(1)`; that subleading formula is not asserted separately
for each geometry at a pooled root. The earlier necessary pure-state
criterion N/m²->0 is now sufficient in this growing-systole regime.

## 5. Original angular U and the two-law root displacement

Keep the same pair of geometries and the original nonzero constant
`Delta_angle=cos4(theta_f)-cos4(theta_s)`, with
`Y=(E_f-E_s)/Delta_angle` and `A_N=N^(13/8)/2`. Since E=1-1_(r=1),

```text
|Y_h| <= N(P1_f+P1_s)/(h |Delta_angle|).
```

Together with (16), for all sufficiently large N this gives

```text
|U_delta| <= [8A_N/|Delta_angle|] max_g epsilon_(N,g)(m_N).        (17)
```

Equation (7) makes the right side smaller than every inverse power of N.
This proves superpolynomial suppression of **original, separately
normalized, pooled-root U** for both fixed laws along N/m²->zeta<infinity.
It is not only a winding-probability estimate or a single-axis proxy.

Finally use `logit(p)=log(m)+log(h)` in (14). For the two corresponding
pooled roots,

```text
logit p_star - logit p_drop = 2log(m)/N+o(1/N).                  (18)
```

The opposite fixed-N tail signs previously established are not contradicted:
both amplitudes vanish in this joint limit, and (17) does not assign their
subleading sign. Nor does their equal leading phase distribution erase the
root displacement (18).

## Scientific card and the remaining fixed-t question

- **Mechanism changed:** the balance is between a Poisson dilute phase and
  the full phase, not generically between two individual configurations.
  Its common partition limit aligns the two roots and closes the pooled
  denominator, including oblique geometry.
- **Exact finite input:** (6)-(7) are explicit all-activity contour bounds
  for both already named laws and all real m>3 in the stated graph scope.
- **Asymptotic outcome:** (3), (14)-(18) take N and m to infinity together
  with N/m² bounded. They are not simulation results or additional evidence
  votes from the existing N25 population.
- **Not proved:** fixed-m oblique suppression, a thermodynamic critical
  line, a continuum field, or a useful finite sampling window for opposite
  U signs. The exp[O(N/m²)] factor is the precise reason this argument does
  not settle arbitrary fixed positive t.
- **Next discriminant:** a fixed-coupling treatment must cancel that bulk
  small-contour partition factor rather than just extend this count or
  reuse the pure-state limit. P154/P334/F4 stop decisions are untouched.
