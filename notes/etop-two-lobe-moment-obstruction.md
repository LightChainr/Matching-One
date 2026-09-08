# The two thermal lobes are not translated copies of one symmetric kernel

The N100/N400 odd rank-clock difference D_A has a visibly two-lobed
canonical readout and strongly evolving standardized rank moments. A natural
mechanism is simply two peaks changing their weights and separation. The
new calculation goes beyond fitting that picture: **even allowing any common
symmetric positive lobe shape, its required sixth moment is negative.**

This is an exploratory inference from the existing source moments, not an
independent new stream or a claim that the signed D_A profile is a probability
density. A positive two-lobe mixture is the candidate being tested. Its
necessary moment inequalities can be applied to signed input without assuming
the input already has that representation.

## A low-cost geometric explanation that does not close

First take two Gaussian lobes with the same within-lobe variance, but arbitrary
relative weight, separation, common center and total width. Standardize the
whole profile to mean zero and variance one. Let w be the weight of the right
lobe, r the fraction of variance in the lobe centers, and

`t=(1-2w)/sqrt(w(1-w))`.

For the standardized central moments m_j, write the cumulant combinations

```
c3=m3; c4=m4-3; c5=m5-10m3;
c6=m6-15m4-10m3^2+30.
```

The two-point centers have cumulants

```
c3=t r^(3/2); c4=(t^2-2)r^2;
c5=t(t^2-8)r^(5/2); c6=(t^4-22t^2+16)r^3.
```

Gaussian noise contributes only its variance 1-r. Thus m3/m4 uniquely
identify the positive r root of `2r^3+c4 r-c3^2=0`, then w. The fifth and
sixth moments are overidentifying predictions, not additional fit parameters.

| Quantity | N100 | N400 |
|---|---:|---:|
| Implied right-lobe weight | .3807896 | .4359509 |
| Between-center variance fraction | .8834633 | .8584904 |
| Fifth-moment observed minus prediction | -.1534149 | -.1703871 |
| Sixth-moment observed minus prediction | -.1987226 | -.0100165 |
| Joint fifth/sixth residual, nominal chi2 / 2 | 848.387 | 177.594 |

Changing lobe weights can summarize part of the decreasing skewness, but it
does not explain the higher moments. The nearly unchanged fifth-moment defect
is descriptive, not a newly established universal invariant.

## Removing the Gaussian assumption leaves a sharper obstruction

Now let `X=B+Z`, where B has two possible centers and Z is **any independent,
mean-zero symmetric positive random variable** common to both lobes. Z can be
discrete or continuous, and its shape and scale may vary freely with N.
Only finite sixth moment is needed. This includes all mixtures of two
translated copies of a common symmetric kernel, not merely Gaussians.

Symmetry makes Z's third and fifth cumulants zero. The measured c3/c5
therefore identify the centers without using any assumption about the
kernel's fourth or sixth moments:

\[
8r+c_5/c_3-c_3^2/r^2=0,\qquad
t=c_3/r^{3/2},\qquad w=(1-t/\sqrt{t^2+4})/2.
\]

For c3 nonzero, the left side is strictly increasing for r>0, so this has
at most one admissible root in (0,1]. Both current mean moment vectors and
all their saved delete-one vectors admit such a root. Nonzero skewness is
needed for this particular identification; the exactly symmetric-total case
would need a different argument.

With v=1-r, the remaining kernel moments are forced:

```
k4 = c4 - (t^2-2)r^2
k6 = c6 - (t^4-22t^2+16)r^3
E[Z^4] = k4 + 3v^2
E[Z^6] = k6 + 15k4 v + 15v^3.
```

Every positive kernel must obey `E[Z^6]>=0` and the Cauchy--Schwarz
moment inequality `v E[Z^6] - E[Z^4]^2 >= 0`. The reconstruction instead gives:

| Necessary kernel quantity | N100 | SE | N400 | SE |
|---|---:|---:|---:|---:|
| v | .07199694 | .00425629 | .03930024 | .0272568 |
| E[Z^4] | .1859464 | .0066983 | .3817423 | .0924043 |
| E[Z^6] | **-1.757800** | .0741887 | **-3.849346** | 1.20234 |
| v E[Z^6] - E[Z^4]^2 | **-.1611323** | .0081322 | **-.2970074** | .0507798 |

The model fails a positivity requirement even before asking whether it
reproduces a pointwise curve. This is broader than the Gaussian residual.
N400's kernel-variance estimate is near the physical boundary; the table
reports full-source delete-one propagated errors, not a specially calibrated
boundary-model likelihood test or an exact confidence certificate.

## Mechanism consequence

Weights, separation and a common symmetric width/shape are not enough in
these source profiles. The remaining possibilities include unequal lobe
kernels, intrinsic lobe asymmetry, more than two components, or signed
geometric cancellation without a positive-mixture interpretation. The
calculation does not select one of these, count continuum fields, or determine
an asymptotic exponent. It does tell the next physical analysis which
information cannot be replaced by merely tracking two peak locations.

This also separates two questions: a quarter-power effective width can still
be a useful finite-size fingerprint even when this full-shape model fails.
The running N900 width experiment keeps its original primary prediction.
These shape formulas were developed from already revealed N100/N400 data;
they are not retroactively called preregistered source tests.

## Scientific card and source

- Changes: the positive common-symmetric-kernel two-lobe class is inconsistent
  with the current standardized moment estimates, even with free parameters
  at each N. It is a mechanism-space reduction, not another harmonic vote.
- Observer/sector: rank-step D_A=Y_(4i)-Y_(2i), normalized odd topology
  response, moments over the full p interval; not local birth insertion.
- Source: `6d8a3ed9d961c66889c3c1e4575485443fdd1c39`,
  `results/p267-standardized-rank-shape/score.json`.
- Dependency groups: same N100 2M-counter and N400 8M-counter blocks as the
  width and peak analyses. The two sizes are independent; the several scores
  within a size are not independent evidence.
- Estimation: recompute center-mixture reconstruction inside every original
  aligned delete-one vector; retain joint parameter/residual covariance.
- Next discriminant: whether the same moment-cone obstruction survives at
  N900; if so, distinguish asymmetric from unequal lobe structure using the
  already collected full profile, without launching a Gaussian-family scan.
- Output: `results/etop-two-lobe-moment-closure/score.json` and `REPORT.md`.
  One deterministic readout of saved moments, no new MC or reliability solves.

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/etop_two_lobe_moment_closure.py
```
