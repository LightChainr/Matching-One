# C04 exact-threshold kappa3 controls

## Stage A: triangular-site self-matching control

The fixed geometry is a 60-degree rhombic torus in triangular-lattice basis
coordinates, with periods `(L,0),(0,L)` and undirected steps
`(1,0),(0,1),(1,-1)`. The observable is black-site wrapping minus
complementary white-site wrapping at the exact threshold `p=1/2`.

Tiny exhaustive results are:

| L | exact kappa3 |
|---:|---:|
| 2 | `-8/9` |
| 3 | `-2424832/1975509` |
| 4 | `-186000605184/136352311523` = `-1.3641177264` |

After a 20,000-replica variance pilot, an independent production seed used
300,000 replicas and 60 blocks per size:

| L | kappa3 | jackknife SE |
|---:|---:|---:|
| 8 | `-1.59309` | `0.03350` |
| 12 | `-1.53251` | `0.03941` |
| 16 | `-1.59666` | `0.05350` |
| 24 | `-1.67803` | `0.07379` |
| 32 | `-1.55638` | `0.07371` |

A weighted fixed-`L^-3/2` model trained on `L=8,12,16,24` predicts `L=32`
within `0.47` standard errors. Fitting all five sizes gives
`kappa_inf=-1.581 +/- 0.046` from Monte Carlo errors alone; `-5/3` is about
`1.85` such standard errors away. A free-power training scan selects its lower
bound `q=0.5` and predicts `L=32` within `0.70` standard errors, so this data
does not identify the correction exponent.

## Cross-control interpretation

The earlier square-bond sequence and the triangular-site sequence agree within
about `0.52`, `2.26`, and `1.31` combined standard errors at `L=8,12,16`.
As an exploratory diagnostic only, a common-intercept fixed-`L^-3/2` fit with
separate microscopic amplitudes gives `-1.636 +/- 0.020` (`chi2=10.31` for 7
degrees of freedom), leaving `-5/3` about `1.50` statistical standard errors
away.

The two observables have not been proved to be the identical universal
normalization, and finite-size systematic errors are not included. Therefore
the present result validates the estimator and retains `-5/3` as a candidate;
it neither establishes universality nor a rational limiting value.
