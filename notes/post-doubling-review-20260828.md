# Post-doubling review: from finite-size fit to Gaussian-semigroup tests

Status: research decision update after server checkpoint
`6d2d68a62b433e337b97eadaf1870cb58e2f7666`.

This note records the strongest current evidence, the claims that remain open,
and the next parameter-free tests. It supersedes the execution order in the
older post-confirmation note but does not alter any raw result.

## 1. Evidence state

### Same-N signal and held-out radial model

The independent 100-million-replica confirmation resolves the prescribed
orientation sign at every declared size `N=65,85,130,145,170`. Pooling the
available seeds gives

\[
A_4 =
\frac{N^{13/8}\Delta M_N}{\Delta\cos(4\theta)}
=0.7885\pm0.0352,
\]

with a small five-size internal chi-square. The frozen Stage-1 model trained on
`N=65,85,130` predicts held-out `N=145,170`; zero effect, the declared power and
log corrections, a free radial exponent, and `cos4+cos8` do not improve the
held-out score.

This confirms a finite-size same-area/same-modulus orientation effect and gives
out-of-sample support to

\[
\Delta M_N\propto\Delta\cos(4\theta)N^{-13/8}
\]

over the present size range. It does not by itself prove that `13/8` is the
unique asymptotic exponent.

### Prospective Gaussian doubling

Multiplication of a Gaussian orientation by `1+i` maps

\[
N\mapsto2N,\qquad\theta\mapsto\theta+\frac{\pi}{4}.
\]

For a pure reflection-even spin-4 sector with site-count exponent `13/8`, the
parameter-free prediction is

\[
\frac{\Delta M_{2N}}{\Delta M_N}
=-2^{-13/8}
=-0.32420988866275241648\ldots.
\]

A fresh post-protocol 100-million-replica run gives

```text
65 -> 130: -0.31382 +/- 0.0908
85 -> 170: -0.34095 +/- 0.1118
```

with covariance-aware joint `chi2=0.03445/2`. Both exact `pi/4` sign flips and
both radial magnitudes pass without fitting an amplitude or exponent.

This is stronger than a conventional two-point exponent estimate because the
same algebraic operation fixes the sign and scale simultaneously.

### Root conversion

Threshold-rank curves give

\[
C_N=-\frac{\Delta p_N^*\,\overline{M'_N}}{\Delta M_N}
\]

between about `0.99984` and `1.00031` at all five pilot sizes. Direct and
linearized root gaps agree. The local conversion from matching residual to
root displacement is therefore validated at current precision.

The earlier test of bare `N^2 DeltaRoot` was not the correct angular invariant:
the Gaussian pairs have different `DeltaCos4`. The relevant amplitude is

\[
A_p(N)=-\frac{N^2\Delta p_N^*}{\Delta\cos4\theta}.
\]

The frozen high-stat target is `A_p=A_M/B`, where
`B=N^{-3/8}\overline{M'}`.

## 2. What the doubling sign does and does not identify

A rotation by `pi/4` sends

\[
\cos(4m\theta)\mapsto(-1)^m\cos(4m\theta).
\]

The observed sign reversal selects the odd-`m` harmonic class:

\[
\cos4\theta,\quad\cos12\theta,\quad\cos20\theta,\ldots
\]

and excludes dominance by the even-`m` class such as `cos8`. It does not alone
distinguish spin 4 from spin 12. Wider angular tests must therefore include an
explicit `H12` alternative rather than treating the doubling result as a
complete harmonic identification.

The existing multi-size design and the prospective `N=185,265` tests remain
the economical route. The four-orientation `N=1105` system becomes useful only
after its power and covariance gates are met.

## 3. Gaussian-integer semigroup spectroscopy

Write an orientation as

\[
g=a+ib,\qquad N=|g|^2,\qquad\theta=\arg g.
\]

For a Gaussian multiplier `h` with `q=|h|^2` and angle `phi`, a spin-`s`
correction with site-count exponent `alpha` transforms as

\[
\delta X_{hg}
=q^{-\alpha}
\operatorname{Re}\!\left(
e^{is\phi} C_s e^{is\theta}
\right).
\]

Thus Gaussian multiplication acts as an exactly known finite-size/angular
operator. Ordinary scaling fields appear as semigroup eigenvectors. A
logarithmic partner appears as a Jordan block.

For the candidate matching-odd field,

\[
s=4,\qquad\alpha=\frac{13}{8},\qquad x=2+2\alpha=\frac{21}{4}.
\]

The `1+i` experiment is its first measured semigroup eigenvalue.

### Exact angular arithmetic

Design files must use exact rational columns:

\[
\cos4\theta
=\frac{a^4-6a^2b^2+b^4}{(a^2+b^2)^2},
\]

\[
\sin4\theta
=\frac{4ab(a^2-b^2)}{(a^2+b^2)^2}.
\]

The new generator `scripts/gaussian_semigroup_design.py` emits these values,
the raw Gaussian product, the canonical D4 representative, lineage order, and
the no-fit scale prediction. It never reads percolation results.

## 4. Full-curve doubling triptych

The next high-value test uses one set of threshold-rank batches to score three
linked ratios:

\[
\frac{\Delta M_{2N}}{\Delta M_N}=-2^{-13/8},
\]

\[
\frac{\overline{M'_{2N}}}{\overline{M'_N}}=2^{3/8},
\]

\[
\frac{\Delta p_{2N}^*}{\Delta p_N^*}=-\frac14.
\]

The third equation follows from the first two. It is the parameter-free
root-level form of the proposed `L^-4` mechanism.

Run all three for

```text
65 -> 130
85 -> 170
145 -> 290
```

with lineage order preserved and full batch covariance. If the raw `-1/4`
ratio fails, also score the finite-slope target

\[
-2^{-13/8}
\frac{\overline{M'_N}}{\overline{M'_{2N}}}
\]

before changing any exponent.

## 5. Direct diagnostic for a logarithmic partner

If

\[
\Delta M_N
=\Delta c_4\,N^{-13/8}(A+B\log N),
\]

then

\[
R_N
=\Delta M_{2N}+2^{-13/8}\Delta M_N
=-2^{-13/8}B\log2\,
\Delta c_4\,N^{-13/8}.
\]

The doubling residual directly estimates the nilpotent/Jordan contribution.
The pure-power null remains primary. Only if several independent lineages show
a coherent residual should one shared `B/A` be fitted.

This diagnostic should be repeated on the thermal-even, matching-odd full-curve
projector, not only at one fixed probability.

## 6. Theory status

At `c=0,h=5/8`, quotienting the level-2 null submodule leaves a non-null
level-4 quasiprimary direction represented by

\[
Q_4=
\left(
40L_{-2}^2-60L_{-3}L_{-1}-9L_{-4}
\right)|h\rangle.
\]

Its bulk chiral combinations have `x=21/4` and spin `+/-4`. The successful
doubling eigenvalue strongly supports this mechanism class.

Still open:

- whether the square-site matching transformation couples to this direction
  with odd parity;
- whether it belongs to a higher-rank logarithmic multiplet;
- whether a lower or degenerate competing odd-harmonic sector is allowed;
- whether the same assignment predicts the derivative parity spectrum.

The self-dual square-bond and C4 self-matching site controls remain necessary.

## 7. Covariance and provenance are mandatory

The existing threshold-rank production reused aligned counter streams across
sizes without mixing `N`. This is legitimate only if treated as deliberate
cross-size coupling and analyzed with its measured covariance.

PR #46 audits the existing batches. Before its scores are authoritative it
must:

- assert equal batch weights or implement a weighted estimator;
- report covariance eigenvalues, condition number, and effective rank;
- use stable factorization/pseudoinverse for near-singular matrices;
- calibrate plug-in chi-square statistics with batch bootstrap or an explicit
  finite-batch approximation;
- retain both full-covariance and diagonal-only results.

Future production must use a clean committed source and executable hash and
must freeze either size-domain-separated RNG or deliberate coupled streams.

## 8. Variance-reduction status

Euler/local-motif controls now reduce individual-geometry matching variance by
about `2.34x` at `N=65` and `2.16x` at `N=85` on fresh evaluation data. That
is a real gain.

The first orientation-difference OLS basis overfit its pilot and remains
excluded. The next control design uses exact zero-mean differences of
equal-multiplicity motifs between the two same-N orientations. Production
adoption is gated on fresh-sample gain for the actual `Delta M` or root target.

## 9. Hardware decision

The fixed-probability confirmation and doubling tests remain CPU work. The
current 16-core implementation already processes hundred-million-scale runs in
minutes.

A GPU is justified only for the bidirectional threshold-rank workload after:

1. clean CPU provenance;
2. exact CPU/GPU equality;
3. stable integer histogram output;
4. measured end-to-end information-per-wall-time gain;
5. a powered sample matrix that is impractical on CPU.

The first GPU target is reusable full-curve sufficient statistics, not large
Pell scans or a speculative `N=1105` run.

## 10. Decision

The project has moved beyond a generic orientation correlation. It now has a
successful prospective, parameter-free semigroup-eigenvalue test:

\[
\boxed{
\Delta M_{2N}=-2^{-13/8}\Delta M_N
}
\]

on two independent Gaussian lineages.

The next objective is not another fitted exponent. It is to close the linked
full-curve triptych, test a third lineage, bound any logarithmic Jordan term,
and then score an overdetermined Gaussian-prime multiplication diagram.
