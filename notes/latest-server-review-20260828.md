# Latest server review: from discovery to mechanism closure

Status: review of PR #21 at server head `28fd1d872f97a83c6ddd59bde996dae2343d3c77`.

This note records the current scientific interpretation after the completed C01, C03, C04, C05 and C07 follow-ups. It supplements `post-server-third-wave.md`; raw server outputs remain unchanged.

## 1. What is now resolved

### 1.1 Small-Pell calibration

The `(7,5)` Pell pair was measured in two independent five-million-replica runs. The pooled diamond-minus-axis root gap is

\[
\Delta p^*=(4.3365\pm0.359)\times10^{-4},
\]

about twelve standard errors from zero. Halving the finite-difference step changed the result by less than `1e-6`, far below the sampling error.

Together with the exact `(3,2)` pair, the two-point effective exponent is about `4.23`, and the `(7,5)` scaled amplitude `L^4 Delta p*` is about `1.06`. This is a positive power calibration for an `L^-4` orientation gap. It does not make large-Pell fixed-p scans efficient; those remain statistically dominated by the rapid decay of the signal.

### 1.2 Same-N Gaussian orientation signal

The 30-million-replica production run used primitive Gaussian quotients with identical `N`, physical area, modulus `tau=i`, cyclic vertex set and counter-keyed random field. For

\[
A_4(N)=\frac{N^{13/8}\Delta M_N}{\Delta\cos4\theta},
\]

the matching-function results are

| N | pair | Delta M | SE | z | A4(N) |
|---:|---|---:|---:|---:|---:|
| 65 | `(8,1)/(7,4)` | `+1.00377e-3` | `1.68448e-4` | 5.96 | `0.65016 +/- 0.10911` |
| 85 | `(9,2)/(7,6)` | `+7.60333e-4` | `1.58350e-4` | 4.80 | `0.65117 +/- 0.13562` |
| 145 | `(12,1)/(9,8)` | `+3.28667e-4` | `1.95415e-4` | 1.68 | `0.55743 +/- 0.33143` |

All signs agree with `Delta cos(4 theta)`. The first two amplitudes agree much more closely than their Monte Carlo precision requires. This promotes the orientation hypothesis to confirmation-stage P0, but it is not yet an asymptotic exponent or operator identification.

The primary working law is

\[
\boxed{
\Delta M_N
=A_4\,\Delta\cos(4\theta)\,N^{-13/8}
+\text{subleading terms}.
}
\]

## 2. What has been rejected or demoted

### 2.1 Matching-even dominance is not a valid gate

The earlier hypothesis that a large matching-even `L^-2 cos(4 theta)` sector must be resolved before the matching-odd residual is interpreted is not supported. The even-sector orientation differences vary strongly with wrapping convention and do not show a stable amplitude collapse across `N=65,85,145`.

This is a completed negative result. It does not erase the resolved matching-function difference, which is channel-identical by exact topology.

### 2.2 Wrapping-only GLS is structurally impossible

For the complementary primal/matching construction, the five validated differences

```text
D_cross
D_both
D_either
D_direction_0
D_direction_1
```

are identical configuration by configuration. Their covariance matrix is therefore rank one and no fixed linear combination can reduce variance. The server evaluation measures exactly `1.0x` variance reduction.

The replacement program is the Euler/local-motif hierarchy in issue #34, based on exact zero-mean local controls rather than duplicate topology views.

### 2.3 Padé corrections do not cure the width drift

The leakage-safe Stage-A challenge selected a degree-four polynomial correction with `n_min=9`; it improves withheld RMSE but retains positive, increasing errors at widths 19--21. The best admissible Padé `[2/2]` is worse and has the same signed trend.

Finite-width predictive accuracy remains distinct from a calibrated infinite-width uncertainty.

## 3. The next decisive test is amplitude closure

The current project has separately observed

- a matching-function residual compatible with `N^-13/8`; and
- a root gap compatible with `N^-2`.

These must be connected using one threshold-rank dataset, not merely compared by exponents.

For two same-N orientations define

\[
A_M(N)=\frac{N^{13/8}\Delta M_N}{\Delta c_4},
\qquad
B(N)=N^{-3/8}\overline{M'_N},
\]

and

\[
A_p(N)=-\frac{N^2\Delta p_N^*}{\Delta c_4}.
\]

Linearization predicts

\[
\boxed{A_p(N)=A_M(N)/B(N)}
\]

or

\[
\boxed{C_N=-\Delta p_N^*\overline{M'_N}/\Delta M_N\to1.}
\]

Issue #35 freezes this closure test. Failure is informative: it separates wrong radial scaling, orientation-dependent thermal metric factors and nonlinear root effects.

## 4. Two-angle data need an angular-radial challenge

Every current `N` supplies only two orientations. The sign pattern and the near equality of the first two `A4` values are strong, but they cannot exclude higher harmonics or logarithmic mixing.

Issue #36 therefore fits the cross-size design

\[
N^\alpha\Delta M_N
=A_4\Delta\cos4\theta
+A_8\Delta\cos8\theta+\cdots
\]

with `alpha=13/8` frozen first, training on `N=65,85,130` and holding out `N=145,170`. This is more information-efficient than immediately attempting a very large four-angle `N=1105` run.

The large multi-angle experiment remains gated on held-out success of the smaller cross-size design.

## 5. Threshold-rank production is now central infrastructure

The Python reference has frozen the exact `K_minus/K_plus` conventions, reconstructed `M`, `M'` and roots from integer histograms, and passed all tiny-system checks. Its `L=8`, 100,000-permutation runtime is `83.85 s`, so production must move to C++ or GPU.

The production engine must preserve, per batch and orientation:

- integer `K_minus` and `K_plus` histograms;
- joint moments and sparse joint histogram where affordable;
- RNG counter ranges;
- enough information to reconstruct `M(p)`, derivatives and roots after the run.

The same run should feed issues #35, #36, #25 and later full-profile analysis. A dense probability grid is not an acceptable substitute.

## 6. Kappa3 requires same-modulus controls

The triangular-site control uses a 60-degree rhombic torus while the square-bond sequence uses a square torus. The limiting derivative ratio can depend on the torus modulus. A common-intercept fit across these two shapes is therefore only an estimator diagnostic, not a universality test.

The current data retain `-5/3` as a candidate but do not establish it. The next useful control must match both the continuum modulus and the observable normalization. Increasing square-bond `L` with a fixed underpowered sample count is lower priority than implementing a same-shape exact-threshold realization or a coordinate-free profile comparison.

## 7. Updated computation order

1. independent confirmation at `N=65,85,130,145,170` under the frozen PR #33 protocol;
2. C++ threshold-rank production engine with Python-oracle regression;
3. issue #35 amplitude closure;
4. issue #36 angular-radial model challenge;
5. issue #34 Euler/motif control variates;
6. same-modulus `kappa3` control;
7. only then a GPU production campaign or the `N=1105` multi-angle test.

The GPU target is not a large-Pell three-point scan. It is batched threshold ranks plus same-N orientation statistics, because those sufficient statistics support several independent scientific tests from one computation.

## 8. Current strongest conjecture

The most economical statement consistent with the data is

\[
M_N(p_c;\theta)
=N^{-13/8}
\left[A_0+A_4\cos4\theta+\cdots\right],
\]

with a nonzero matching-odd spin-4 amplitude. A natural but unproved continuum candidate is an `x=21/4`, spin-4 field, possibly a level-four thermal-family descendant or quasiprimary in the `c=0` logarithmic theory.

That interpretation is promoted only after the angular-radial and amplitude-closure tests succeed. The immediate empirical claim is narrower: a same-area, same-shape orientation signal has been resolved with the predicted sign and a first two-point radial collapse.