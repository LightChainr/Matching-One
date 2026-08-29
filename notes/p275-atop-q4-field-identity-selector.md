# From harmonic identity to field identity in the A_top channel

## The change of question

The quotient prism at `fc14817` ended the useful H4/H8/H12 competition for
the global matching channel.  H4 fits (`chi2=1.088/2`) while the two fixed
aliases fail.  The double projector at `ddf41aa` then places thermal
`Q4 epsilon` first in the listed singlet/Alexander-odd staircase, and `078bd61`
constructs an exact nonzero finite coupling between `A_top` and a local
rank-birth H4 source.

The missing step is no longer another angular vote.  It is the identity of the
H4 field that carries the lattice-to-continuum coupling.

## Why the committed archives stop one statistic short

The existing objects were inspected by observable semantics, not by filename.

- `fc14817` contains full threshold-rank histories for `A_top=P2-P0` and
  selects H4 on the square modulus.  It does not retain the birth line `ell`,
  a local H4 source, or the same-sample product `A_top*J_D4`.
- Gaussian norm-2/norm-4/norm-5 Hecke children add scale and spin-character
  information but remain covers of the square modulus.  Their marginal rank
  histograms cannot recover a connected marked response.
- `cc1d43c` has a genuine complex character and a complete 4D covariance, but
  its observable is the charged local response `O_chibar S_chi`, not the
  unlabelled global `A_top`.  Relabelling it would undo the double-projector
  semantics.
- `078bd61` proves the exact response identity and tiny nonzero controls, but
  explicitly records that the production archives lack `qJ` cross moments.

Therefore the field-identity score cannot be reconstructed.  More analysis of
the same marginals would only manufacture a normalization.

## The one observable that connects all three results

For every Bernoulli configuration use

```text
q = rank image[H1(black)->H1(T2)] - 1 = A_top,
J_D4 = sum_v chi4(ell_v) (I12-I01)_v,
B = E[sum_v(I01+I12)_v].
```

The exact exponential-source response is

```text
Cov(A_top,J_D4).
```

Both factors are complement odd, so their connected covariance is allowed.
Divide by the unmarked birth mass,

```text
gamma_D4 = Cov(A_top,J_D4)/B.
```

This removes the total thermal/rank-pivotal density while preserving the
complex H4 phase.  It is the lattice version of the missing matrix element
`g[A_top,Q4 epsilon]`, with the important boundary that a non-Q4 H4 field can
also contribute.

For a period frame rotated by the Gaussian representative `z=a+ib`, transport
back to the canonical modulus frame exactly:

```text
Gamma = conjugate[(z/abs z)^4] gamma_lab.
```

All phase coordinates in the machine artifact are rational complex numbers;
no fitted angle is introduced.

## A three-modulus Hecke/eta triangle

Use the physical moduli

```text
tau0 = i,
tau1 = 2i,
tau2 = 5i/2.
```

The first edge retains the exact degree-2 Hecke result

```text
E4hat(2i)/E4hat(i)=11/4,
E4hat(tau)=Im(tau)^2 E4(tau).
```

The third point is not modularly equivalent to either endpoint.  Direct
90-decimal q-products give

```text
E4hat(5i/2)/E4hat(i) = 4.2934368543749230211...,
2 log|eta(2i)/eta(i)| = -0.5198603854199589821...,
2 log|eta(5i/2)/eta(i)| = -0.7816530999014344284....
```

Thus the gauge-free energy-block cocycle is the genuine three-modulus relation

```text
[R(2i)-R(i)]/[R(5i/2)-R(i)]
  = 0.6650781343866131773...,
```

not a basis-transformed copy of the square point.

## Small cyclic geometries, not a quotient contest

The three sizes are `N=50,130,170`.  At each size the square, aspect-2 and
aspect-5/2 period matrices all have Smith invariants `(1,N)`.  This removes the
cyclic/noncyclic quotient change from the new selector.

| N | tau=i | tau=2i | tau=5i/2 |
|---:|:---|:---|:---|
| 50 | `[[7,-1],[1,7]]` | `[[4,-6],[3,8]]` | `[[4,-5],[2,10]]` |
| 130 | `[[11,-3],[3,11]]` | `[[8,-2],[1,16]]` | `[[6,-10],[4,15]]` |
| 170 | `[[13,-1],[1,13]]` | `[[9,-4],[2,18]]` | `[[8,-5],[2,20]]` |

All three geometries at a fixed N use one priority field, random-root counter,
replica interval and batch partition.  Different N use independent seeds.

## The field-identity model spaces

Set

```text
Y(N,tau)=N^(13/8) Gamma(N,tau),
F(tau)=E4hat(tau)/E4hat(i),
e(tau)=2 log|eta(tau)/eta(i)|.
```

The frozen subspaces are:

```text
ordinary thermal Q4:
  Y = c F(tau)                                  one complex c

thermal Q4 energy-Jordan completion:
  Y = F(tau) [c0+c1 log N+c2 e(tau)]            three complex c

other pure H4 completion:
  Y = c_tau                                     three complex c

other affine-log H4 completion:
  Y = a_tau+b_tau log N                         six complex c.
```

The Jordan model is not “ordinary plus a free log”.  Its log-slope modulus
vector is forced to be proportional to `E4hat`, and its normalized intercept
differences are forced onto the eta vector.  An arbitrary allowed H4 field is
given enough freedom to fit its own modulus intercepts/slopes; if that model
wins, the outcome is a different H4 completion, not a reason to reopen H8.

The observation vector contains nine complex means, hence 18 real coordinates.
The Q4 ordinary/Jordan/generic-pure/generic-affine models have respectively
2/6/6/12 real nuisance parameters.  They are nested mechanism subspaces, not
18 marginal votes.

## Complete sufficient statistics

One uniform random root per configuration, multiplied by N, is an unbiased
site-sum estimator.  Retain in every aligned batch and geometry:

```text
sum q, sum q^2,
sum I01, I12, I02,
sum Re/Im J_S4 and J_D4,
sum q Re/Im J_S4 and q Re/Im J_D4,
unmarked N*(I01+I12) birth mass,
priority-field digest and random-root digest.
```

The S source is the complement-parity control.  A landing-H4 version may be
stored in the same rows as a correlated diagnostic, but it is not a second
primary score.

Every connected covariance, division by B, exact phase transport, N scaling,
model coefficient and residual must be recomputed inside synchronized
delete-one batches.  Archive the full 18x18 real covariance in frozen order.
Marginal histograms, unmarked birth mass alone, or marked means without `qJ`
are not sufficient.

## Decision map

- Ordinary Q4 fits and Jordan gives no material log direction: identify the
  field as the ordinary thermal `Q4 epsilon` completion.
- Ordinary fails, Jordan fits, and the log direction is nonzero: identify the
  inherited thermal Q4 Jordan completion.
- Both Q4 subspaces fail but a generic H4 subspace fits: retain H4 and replace
  its proposed field identity.
- Every nonzero subspace fails: the flow is multi-field or `J_D4` does not
  isolate the global H4 coupling.

This branch freezes the selector and a 20M-per-geometry first acquisition but
does not implement the runner, launch production, open a PR, or add any new
harmonic hypothesis.
