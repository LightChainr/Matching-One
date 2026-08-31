# P334: first-step direct absorption and the Doob clock innovation

Let a fixed rank-one prefix have d remaining labels and h current direct
absorbing sites. The next label is uniform. Let D mean that it is one of
those h labels, and p=h/d. The full suffix stopping time T is at least one.

For any terminal clock readout X=g(k0+T), write m=E[X|prefix] and
a=g(k0+1). All direct choices give X=a. If 0<p<1, conditioning only on the
binary direct event yields the exact identity

\[
B=\operatorname{Var}(E[X|D,\mathrm{prefix}]\mid\mathrm{prefix})
=\frac{p}{1-p}(a-m)^2.
\]

The complete next-label sigma-field is finer, so its Doob innovation
variance is at least B. Moreover

\[
\boxed{\operatorname{Var}(X|\mathrm{prefix})
=B+(1-p)\operatorname{Var}(X|D=0,\mathrm{prefix}).}
\]

Thus B is not just an arbitrary lower bound: it is the exact variance
explained by the already named direct-versus-safe split. The remainder
requires finer safe-label information and/or subsequent suffix information.
It is **not** a lower bound on variance remaining after the whole next label
is revealed, because the safe label's identity may explain some or all of it.

## Integrated clock and boundary cases

Put mu=E[T|prefix]. For the integrated matching-clock readout
X_int=1-(K1+k0+T)/(N+1),

\[
m_{int}=1-(K_1+k_0+\mu)/(N+1),\quad
B_{int}=\frac{h}{d-h}\frac{(\mu-1)^2}{(N+1)^2}.
\]

After revealing the next label, let mu_child be the expected *additional*
wait, with value zero on direct absorption. The tower property gives
E[mu_child|prefix]=mu-1. The next-label innovation is therefore
Var(mu_child|prefix)/(N+1)^2, yielding the same lower bound. The archived
integrated F2 readout differs from X_int only by the known prefix constant
K1/(N+1), so its conditional variance and innovation are identical.

- h=0: D is always safe and B=0. This does not force the full next-label
  innovation or suffix variance to be zero.
- h=d: every label absorbs, T=1 and all conditional variances vanish. Set
  B=0 directly, never evaluate the apparent 0/0, and do not form a prefix
  fraction. This also implies mu=1.
- If a nontrivial terminal readout is constant on its supported stopping
  times, its variance is zero; its prefix ratio is likewise undefined.

## Archive-only measurement definition

Use all forty original batches from `0d1e586d`, at N325/k0=193 and
N425/k0=252. Only already stored complete safe coefficient vectors are
consumed. For S(j)=f_j/binomial(d,j),

\[
E T=\sum_{j=0}^{d-1}S(j),\qquad
E T^2=\sum_{j=0}^{d-1}(2j+1)S(j).
\]

The canonical secondary readout is exactly the archive's F2 kernel,
g(k)=P(Bin(N,p_ref)>=k), p_ref=.59274605079. Its conditional moments are
evaluated from S(j-1)-S(j), without reconstructing a single child or network.

Per size and orientation, report sum B/sum V. This is a variance-weighted
fraction, not the average of per-prefix fractions. Completed univariate
clocks retained inside a whole-pair fallback may be used for this marginal
diagnostic; their coverage is explicit and no paired hybrid substitution
is changed. Unsolved rank-one prefixes are missing, not zero variance.
Delete-one-original-batch covariance retains both orientations and both
readouts within each size.

Two orientation bounds cannot simply be added into an H4 contrast bound.
That requires their common next-label event intersection and conditional
covariance. The present result makes no such paired claim.

## Completed forty-batch readout

The direct-versus-safe event explains a modest but definite part of the
clock's conditional suffix noise. Numbers are ratios of summed bound and
summed conditional variance; +/- denotes the original-batch standard error
in percentage points, not uncertainty in the exact clock coefficients.

| N / orientation | Stored clocks / rank-one rows | Integrated fraction | Canonical fraction |
|---|---:|---:|---:|
| 325 first | 8973 / 8997 | 8.2614% +/- .0272 | 22.3434% +/- .0351 |
| 325 second | 8932 / 8966 | 8.3588% +/- .0362 | 22.4428% +/- .0403 |
| 425 first | 8830 / 8910 | 6.9320% +/- .0276 | 19.3629% +/- .0357 |
| 425 second | 8959 / 9081 | 6.8622% +/- .0337 | 19.3264% +/- .0456 |

Mean direct probabilities are .062683/.063698 at N325 and
.052900/.052455 at N425. Mean waits are 11.7054/11.6167 and
13.8645/13.9620. The canonical readout gives more weight to the earliest
clock outcomes than the integrated affine readout; the same first-step
event therefore explains a substantially larger fraction of its noise.
This is a readout-dependent result, not a fitted scaling law.

What remains after observing this binary event is 91.64%-93.14% of the
integrated conditional noise and 77.56%-80.67% of the canonical noise. It
must be resolved inside the safe event, by the identity of the safe next
label and/or by later insertions. The stronger full-next-label innovation
cannot be recovered from these marginal clock coefficients alone.

For transparency, prefixes with h=0 number 665/619 at N325 and 649/708 at
N425. Their floor is zero, but they carry 7.98%-9.23% of the integrated
suffix variance and 3.18%-3.59% of the canonical suffix variance in the
four solved pools. No stored solved prefix has h=d.

The readout includes 12 and 38 completed first-orientation clocks retained
inside whole-pair fallbacks at N325 and N425. They are valid for this
single-orientation conditional calculation; no paired substitution rule
was changed. Count coverage ranges from 98.6565% to 99.7332%. This does
not bound the missing *variance mass*: unresolved prefixes are not imputed,
so the displayed fractions describe the solved archive pool, not a
guaranteed estimate for every rank-one prefix.

`results/p334-first-step-doob-clock-innovation/batches.csv` preserves all
40 source batch identities, counter intervals and orientation sums (80
rows). `score.json` preserves the full joint 12-statistic delete-one-batch
covariance per size, including both kernels and both orientations, and
SHA256 hashes of every unchanged source archive. No cross-size pairing
or H4 variance lower bound is manufactured.
