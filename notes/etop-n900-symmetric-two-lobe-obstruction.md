# N900 still rules out the common-positive-symmetric-kernel two-lobe moment picture

At the completed N900 source moment estimate, **two translated copies of any
common positive symmetric kernel require a negative sixth kernel moment**:
`E[Z^6] = -2.142812 ± 0.481800`. The accompanying moment determinant is
`Var(Z) E[Z^6] - E[Z^4]^2 = -0.260440 ± 0.059149`.
Both must be nonnegative for a positive kernel. This is the same necessary
condition as `ddf7d564`, with no new kernel family or branch choice.

The branch question is resolved on this archive: the point estimate and all
800 aligned delete-one estimates have the unique admissible branch; every
one has both negative quantities. Thus the small early weight in N900's
three-center chart does not rescue this broader two-lobe mechanism.

## Fixed moment condition and actual readout

Standardize the whole candidate profile to zero mean and unit variance.
Write `X=B+Z`, with two-valued B and independent mean-zero symmetric Z,
common to both lobes. Z may be discrete or continuous; finite sixth moment
is the only moment-existence requirement. With between-center variance r,
right-center weight w, `t=(1-2w)/sqrt(w(1-w))` and `v=1-r`, set

```
c3=m3; c4=m4-3; c5=m5-10*m3;
c6=m6-15*m4-10*m3^2+30.
8*r+c5/c3-c3^2/r^2=0;  t=c3/r^(3/2).
k4=c4-(t^2-2)*r^2;
k6=c6-(t^4-22*t^2+16)*r^3;
E[Z^4]=k4+3*v^2;
E[Z^6]=k6+15*k4*v+15*v^3.
```

For nonzero c3 the root equation is strictly increasing on r>0. The frozen
solver and admissible interval `0<r<=1` are unchanged. No alternative root
is chosen when a source row fails; such rows would prevent reporting a
valid-only covariance as if all 800 branches had been retained.

N900 measured `(m3,m4,m5,m6)` is
`(0.3385558273,1.7068063682,0.9630252133,3.8922101618)`.

| Required coordinate | Estimate | Full delete-one SE |
|---|---:|---:|
| Right-center weight w | 0.4045556 | 0.0259374 |
| Between-center variance r | 0.9116745 | 0.0208002 |
| Kernel variance v | 0.0883255 | 0.0208002 |
| Kernel fourth moment | 0.2667866 | 0.0443816 |
| Kernel sixth moment | **-2.1428116** | 0.4818003 |
| Moment determinant | **-0.2604400** | 0.0591488 |

The 800-LOO branch proportions are: admissible unique 100%, no admissible
`r<=1` branch 0%, degenerate zero-skew identification 0%, other failure 0%.
LOO sixth moments range from `-2.195096` to `-2.088424`; determinants range
from `-0.267131` to `-0.254826`. These are highly overlapping leave-one
estimates, not 800 independent tests or confidence limits.

## Mechanism consequence and boundary

The shape obstruction now survives the N900 extension, allowing the common
kernel's shape and variance to vary freely with size. Tracking two centers,
their weights and an arbitrary common symmetric blur is insufficient for
the empirical moment vector. Unequal kernels, asymmetric lobes, additional
components and signed geometric cancellation remain distinct possibilities;
this readout does not choose among them.

The earlier `3bacf19a` rank-three gap concerns the narrower
Gaussian-plus-two-center boundary. It explicitly supplied a common-symmetric
two-lobe counterexample with a nonzero Gaussian rank-three gap. The present
calculation applies the separate broader-class condition directly to N900;
it does not infer the broader result from that gap or from the 3.2% early
weight. These are dependent readouts of the same source.

All inference here is auxiliary and exploratory on the signed rank-profile
moments. Positivity is a property demanded of the candidate representation,
not assumed for the observed profile; source bins were not clipped. The two
negative margins are about -4.45 and -4.40 ordinary propagated SEs, respectively.
They are not two independent discoveries, a boundary-calibrated likelihood
test, or an exact production-confidence certificate. No conclusion about a
continuum field count, an asymptotic exponent or the primary width prediction
is made.

## Scientific card and reproducibility

- Changes: the N900 extension preserves the common-positive-symmetric
  two-translation moment obstruction; numerical branch ambiguity does not
  explain the result.
- Observer/sector: standardized moments of the signed odd topology
  rank-profile `D_A=Y_(4i)-Y_(2i)`; not a local insertion observable.
- Source: production `5f30397c`, `results/etop-n900-rank-width`, 32M per
  shape pair and 800 common batches. Shape archive `54430ea7`,
  `results/p267-max-gaussian-three-center-n900/score.json`.
- Dependencies: precisely the same N900 common-batch group as the width,
  three-center and boundary-gap outputs. No new source or combined-scale
  evidence score.
- Next useful discriminator: a physically named unequal-kernel or asymmetric
  readout that explains the odd/even remainder; the current analysis neither
  starts a family scan nor requests new samples.
- Lifecycle: fixed `ddf7d564` gate -> saved N900 point/800 LOO -> full aligned
  covariance -> auxiliary mechanism obstruction; zero new MC, zero source
  constructions and zero N100/N400 reruns.

The point moments are read directly from the frozen shape archive. Its LOO
records contain already-matched three-center coordinates, so m3..m6 are
recovered by forward algebraic moment evaluation, without reconstructing or
refitting that source. The source reports maximum matched-moment error
`8.22e-15`. This compressed representation is used only to recover its input
moments, not to assume Gaussian lobes in the model being tested.

`results/etop-n900-symmetric-two-lobe-moments/score.json` preserves the source
hashes, all 800 recovered moment vectors, branch statuses, derived coordinates
and their complete covariance. Reproduce from this repository's fixed objects:

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/etop_n900_symmetric_two_lobe_moments.py
```
