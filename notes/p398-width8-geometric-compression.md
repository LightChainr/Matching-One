# P398: named pair geometry nearly closes the crossing, not the whole tail

Geometry spans and comparison distances were fixed at `588a2fdd` before
these new projection outputs. They operate on the already known width-eight
generator and correlations, so this remains retrospective mechanism analysis,
not a fresh experiment. No decay, amplitude, source or lag is fitted.

**Outcome:** T3/S11/B2 and their Kreweras partners give a 5+5-dimensional
geometric dynamics with crossing error .218% and plus slow-mass error .512%.
Triplet boundary incidence improves the crossing and minus tail, but worsens
the plus tail. The next specifically missing geometry is size-four charge;
the right and left spectral diagnostics do not support calling it a complete
explanation yet.

## Fixed named dynamics, with exact relations

Every seed comes with its Kreweras pullback. The fixed stages are:

| Stage | Seeds before adding partners | Exact dimension per ray |
|---|---|---:|
| Source | A | 1 |
| First force | A,T2 | 2 |
| Pair hierarchy | A,T2,T3,S11,B2 | 5 |
| Triplet incidence | previous seeds plus T3_r0,T3_r1,T3_r2 | 7 |

Here r is the number of nearest neighbors of a triplet member in the same
frontier block. The exact relations are `T3=T3_r0+T3_r1+T3_r2` and its
Kreweras pullback. These are explicitly removed, not additional fit degrees
of freedom. All other displayed columns are independent. A modular nonzero
minor and Gaussian-integer null vectors establish these finite column ranks.

The result saves the full effective mass matrix **in independent named
columns in declaration order**, not only in an opaque spectral basis.
Numerical orthonormalization merely represents the same prescribed span.
The matrix is the stationary-L2 Galerkin projection of M=-G; source
normalizations remain exactly those of the inherited A/L readouts.

The pair hierarchy contains M^2 psi, because

```text
GA=-3A+R; GL=-3L+T2;
GT2=2T3-2T2+S11-B2,
```

and Kreweras symmetry supplies the corresponding GR relation. Therefore
the normalized source moments through M^3 are exact in the 5+5 model;
this includes both initial memory curvature and its first derivative.
No exact closure of the next moments or long tail follows.

## Which propagation is explained?

The complete reference crossing is .2656573200; its lowest masses are
2.819658633 in minus and 1.955750138 in plus.

| Model | Crossing | Crossing error | Minus mass | Plus mass |
|---|---:|---:|---:|---:|
| Source only | none in fixed bracket | — | 3.368820241 | 3.691415268 |
| First force, 2+2 | .2541213924 | -4.3424% | 3.069281688 | 2.005141014 |
| Pair hierarchy, 5+5 | .2650793593 | -.2176% | 2.878004560 | 1.965757683 |
| Triplet incidence, 7+7 | .2654226141 | -.08835% | 2.842890287 | 1.931369756 |

The identical frozen lag grid was .05,.1,.25,.5,1,2,4. Relative errors
of the normalized covariance in its long end are:

| Model | Minus at t=2 | Plus at t=2 | Minus at t=4 | Plus at t=4 |
|---|---:|---:|---:|---:|
| First force | -17.360% | -.2340% | -47.588% | -9.385% |
| Pair hierarchy | -2.850% | -.1818% | -11.830% | -2.094% |
| Triplet incidence | -.6614% | +2.920% | -4.011% | +8.062% |

Thus pair creation, pair loss at the cut boundary and triplet chipping are
enough for a very accurate crossing and a useful plus slow-mode approximation.
Resolving triplet boundary incidence is genuinely helpful for minus: the
slow-mass error falls from +2.069% to +.824%, and the t=4 error from -11.83%
to -4.01%. Its strong signed budgets were not merely algebraic decoration.

But the same refinement pushes the plus mass **below** the true value by
1.247%, producing a long-tail overshoot. Source variance is unchanged, and
all low source moments through M^3 remain correct. This is truncation of a
nonselfadjoint dynamics, not a newly identified physical slow state.
Stationary-L2 Galerkin masses need not converge monotonically or bound the
true mass.

## The residual is specifically size-four charged geometry

Only the pre-named next hierarchy columns T4,Q3,B3 were diagnosed, through
`GT3=3(T4-T3)+Q3-B3`; none was added to the scored models. Q3 joins adjacent
singleton/pair blocks into a triplet, and B3 is triplet charge weighted by
its cut-boundary edge count.

After the final 7+7 projection, the plus-direction diagnostics are:

| Named residual | Remaining fraction of its variance | Alignment squared with omitted slow right mode | With omitted slow left mode |
|---|---:|---:|---:|
| T4 | .59972 | .64800 | .12562 |
| Q3 | .14355 | .28840 | .05771 |
| B3 | .00948 | .00478 | .00087 |

Size-four charge is an explicit independent direction, and the strongest
of these three candidates for the missing plus right mode. The left-mode
alignment is much smaller, so this does **not** predict that adding T4 alone
must repair the tail. In minus, T4's corresponding right/left alignments are
.31278/.00312, again exposing the importance of the adjoint side.

Indeed, moving from 5+5 to 7+7 increases plus left-mode capture from 95.49%
to 97.73% and right-mode capture from 98.89% to 99.15%, yet the mass overshoots.
One-sided overlap, or even separately improved left/right overlaps, does not
fix the small nonselfadjoint coupling error. This distinguishes a promising
new geometric variable from an already established closure.

## Reproduction and boundary

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /Users/lc/python-envs/research-py311/bin/python scripts/p398_width8_geometric_compression.py
```

One single-thread deterministic calculation emits exact column relations,
named matrices, all fixed-lag responses, masses and two-sided diagnostics.
It does not rerun the 93+93 rank campaign, add samples, scan alternative
training bases, or operate a server. The result concerns the existing finite
positive frontier process, not continuum field counts or Jordan identity.
