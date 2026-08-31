# A thinning early shoulder retains its rank-three moment leverage

The early-center weight falls `0.1806 -> 0.0654 -> 0.0320`, but this alone
does not imply convergence to Gaussian plus two centers. One fixed invariant
was defined in `fa562f47` before evaluating it on the stored three-scale LOO
vectors. No Gaussian boundary or center was refitted, and no new sample used.

## A unique, affine-invariant rank-two residual

Let `C` be the positive center measure of the existing maximal-Gaussian
degree-six realization, `Y=C-E[C]`, `v=E[Y^2]>0`, and `q_j=E[Y^j]`. Define

\[
\Delta_2=\min_{a,b}\frac{E[(Y^2-aY-b)^2]}{v^2}
=\frac{q_4}{v^2}-1-\frac{q_3^2}{v^3}
=\frac{\det H_2}{v^3},\qquad
P_2(Y)=Y^2-\frac{q_3}{v}Y-v.
\]

`H_2` is the 3-by-3 moment Gram matrix in `(1,Y,Y^2)`. For three centers,
its Vandermonde factorization gives the exact expression

\[
\Delta_2=\frac{w_0w_1w_2\prod_{i<j}(c_i-c_j)^2}{v^3}
=\frac{w_0w_1w_2\rho^2(1-\rho)^2}{v_a^3},
\]

where `rho=(c1-c0)/(c2-c0)` and `v_a` is the weighted variance of `(0,rho,1)`.
This is one residual, expressed three ways, not three selected distances.
Translation/scaling cancel exactly. It is also unchanged by common Gaussian
blur of the source because it depends on center weights/relative geometry,
not the Gaussian fraction `alpha`. It is not invariant under arbitrary
nonlinear thermal warps, nor is it a metric such as Wasserstein distance to
the whole two-component model class.

For a positive measure, `Delta2=0` means that the monic quadratic vanishes
almost surely: there are at most two distinct centers. Thus, in the finite,
nonzero-center-variance closure of this realization, it is exactly the
**Gaussian-plus-two-center boundary through degree six**. Conversely, if a
source already has two centers and Gaussian variance `T`, inverse heat at
any slightly larger variance makes
`d E[P2^2]/dt = -E[(P2')^2] < 0`; hence `T` is the maximal feasible variance
and the same boundary is found. This is a truncated-moment equivalence,
not a claim that matching six moments identifies the whole signed profile.

## Three-scale readout

| N | Delta2 ± LOO SE | early weight | leverage Delta2/w0 ± SE |
|---|---:|---:|---:|
| 100 | 0.2096771 ± 0.0104424 | 0.1805891 | 1.16107 ± 0.06774 |
| 400 | 0.1732677 ± 0.0287416 | 0.0653620 | 2.65089 ± 0.51658 |
| 900 | 0.1727232 ± 0.0411192 | 0.0320100 | 5.39592 ± 1.62250 |

The N400-to-N900 gap change is `-0.0005446 ± 0.0501684`; from N100 to N900
it is `-0.0369540 ± 0.0424244`. There is no resolved decrease of the invariant
gap, despite the much smaller early weight. The identity
`Delta2=w_early*leverage` exposes the compensation: lower weight retains
substantial influence on the residual because its relative center leverage
increases. The ratio is explanatory algebra, not a second optimized distance.
LOO covariance includes the strong dependence among gap, weight and ratio.
The rank-two null is singular, so these regular-branch SEs are not turned
into a purported calibrated normal boundary p-value.

With finite noncolliding centers and other weights bounded away from zero,
`w0->0` does force `Delta2->0`. But small weight alone does not provide those
uniformity conditions: a moving low-weight shoulder can retain fourth-moment
leverage. The measured gap, not a 3% component label, distinguishes those
possibilities. Nothing here says that a third physical field is disappearing.

## Why this is not the earlier common-symmetric two-lobe test

`ddf7d564` allowed **any** common positive symmetric kernel `Z` in `X=B+Z`
with two-valued `B`; its odd-moment reconstruction forced impossible kernel
sixth moments on N100/N400. That class is broader than Gaussian-plus-two
centers. A positive rank-three gap alone does not exclude it.

An exact counterexample makes the distinction concrete: take independent
`B=+-a` and `Z=(+-a)+G`, where both signs have equal probability and `G` is
Gaussian. `Z` is symmetric, and `X` is two translates of this one kernel.
Yet its Gaussian-deconvolved centers are `(-2a,0,2a)` with weights
`(1/4,1/2,1/4)`, giving **Delta2=1**, not zero. The broader old no-go and this
Gaussian-center boundary are therefore different questions. Neither counts
visible peaks or physical fields. A Fraction-based focused check verifies the
Hankel/energy/Vandermonde identities, affine invariance, two-center zero and
this exact counterexample.

Scientific card: the new positive interpretation is **persistent polynomial
shape leverage of a thinning early shoulder**, rather than automatic reduction
of state count. The next discriminator is the already named gap together with
its weight/leverage factorization at subsequent scales: compact two-center
approach versus vanishing weight without moment collapse. No new collection
is authorized by this note. Sources are the same N100/N400/N900 blocks in
`54430ea7`; all 1400 stored LOO vectors are reused, not new independent evidence.

Reproduce: `python3 scripts/p267_two_center_boundary_gap.py`.
Outputs: `results/p267-two-center-boundary-gap/{score.json,REPORT.md}`.
