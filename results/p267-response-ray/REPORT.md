# P267: changing amplitude does not close the four-coordinate response

The best common amplitude is `3.321708`, but the two fixed-p response vectors still disagree at nominal `chi2=19.87177/3, p=0.0001804517`.
This excludes the declared finite four-coordinate amplitude-only model. It does not count microscopic sources or identify a continuum operator.

## Definitions and evidence

x and y are the exact cos4-normalized `(A_top,E_top,C,W)` contrasts at tau=i and 2i, fixed N50 and p=.59274605079. C and W are normalized integrated clocks. The old 12M square block and new 100k rectangular block are independent; each retains its internally paired direction covariance.

| coordinate subset | fitted amplitude | residual chi2 / nominal df | nominal p |
|---|---:|---:|---:|
| A/E/C/W | 3.321708 | 19.87177 / 3 | 0.0001804517 |
| A_top/E_top | 3.006755 | 0.1967398 / 1 | 0.6573654 |
| C/W | 3.410415 | 2.041993 / 1 | 0.1530091 |
| A_top/C/W | 3.485093 | 2.143139 / 2 | 0.3424706 |
| E_top/C/W | 3.37264 | 19.82615 / 2 | 4.952295e-05 |
| A_top/E_top/W | 4.332734 | 17.7898 / 2 | 0.0001370862 |
| A_top/E_top/C | 3.52998 | 0.7805164 / 2 | 0.6768821 |

These subsets localize one shared result, not independent discoveries. Their nominal p-values are retrospective and not multiplicity-adjusted. The source vectors and complete joint covariance are saved in score.json.

## Model and uncertainty

The joint Gaussian model is `x~N(mu,Cs), y~N(lambda*mu,Cr)`. Profiling mu gives `D(lambda)=(y-lambda*x)^T(Cr+lambda^2*Cs)^-1(y-lambda*x)`. The original joint covariance determinant is constant, so no contrast log-determinant is added. A compact projective-angle search includes lambda=infinity and refines the grid minima.

Fixed lambda=1 gives `236.75639` on nominal df4. Allowing one amplitude improves distance by `216.88462` on nominal df1, but leaves `19.871771` on df3. The infinity boundary has distance `16863.787`.

The 95% conditional amplitude profile interval is `[3.0087921019116686, 3.63660159784094]`. This is a parameter range inside a rejected model, not an identified physical transport coefficient. All inference conditions on the saved estimated covariance; no exact finite-sample coverage is claimed.

## Where the four-coordinate mismatch remains

At the full fitted amplitude, A_top/E_top marginal distance is `0.3990675`; C/W conditional on it contributes `19.4727`. Their sum is the full distance.
At the full fitted amplitude, C/W marginal distance is `2.138045`; A_top/E_top conditional on it contributes `17.73373`. Their sum is the full distance.

These order-dependent Schur partitions describe correlated residuals; they are not causal shares or separate tests. Two observed vectors always fit some rank-two plane, so a saturated rank-two fit would add no mechanism information.

## A third modulus needs a realizable geometry

At N50, Smith(5,10) forces L=5B with det B=2. Up to square-lattice symmetry, the index-two sublattices provide only the square and 1x2 rectangle. A third same-N, same-Smith modulus is therefore not available. A new map must change a concrete endpoint/readout, not merely rename a Smith class.

The saved exact future design uses N250 and Smith(1,250)/(5,50), with the same O=(1/5)[[4,-3],[3,4]] at rho=1,2,5. Its integer matrices, determinants, orthogonal columns and exact direction factors are recorded in score.json. It is not acquired, powered or a request to start production.

Two explicit exploratory shape candidates are `f_log(rho)=(log(rho)/log(2))^2` and `f_geom(rho)=2(rho+1/rho-2)`. Both interpolate the two anchors and are even under axis exchange. At rho=5 they predict coefficients `f_log=5.3913501` versus `f_geom=6.4` in `u5=(1-f5)u1+f5*u2`. These are distinguishable phenomenological hypotheses, not laws selected by the old two points.

Any N250 prediction needs N250 anchors or an explicitly justified coordinate-wise scale law. Current N50 vectors are not numerical predictions at N250. Likewise, fixed-p rows must not be concatenated with intrinsic-center norm-4/norm-5 production rows without recomputing the same observable definition and covariance.

## Next scientific output

Preserve the amplitude-only exclusion. Use the coordinate localization to define the next physical response model; compare its explicit predictions on a realizable design. Do not treat an arbitrary rank-two interpolation, a free surface, another compiler or more N50 factorial replicas as identification.

## Reproduce and provenance

Input: `ce01e4d10abb03abf1b278192510937ee96db29d:results/p267-etop-tau-topology-factorial/score.json`; SHA256 `379eedf05559f259a237b10bc0a5b0d3e26540cd78f222218d4d31a56ae4467a`. This is a dependent, post-reveal reanalysis of the existing factorial. No Monte Carlo or test suite was run.

```sh
python3 scripts/analyze_p267_response_ray.py
```

One in-run arithmetic check recovers the parent's fixed-amplitude distance and both Schur sums. No chart is needed: the small exact-value table exposes the coordinate comparison directly. All raw covariance, fitted residuals, profile intervals and prospective integer designs are in the companion JSON.
