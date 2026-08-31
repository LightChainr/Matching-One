# Continuous center–lifetime moments from the common-label fork archive

## New product

The same forty original prefix batches now carry the unperturbed **continuous**
center/lifetime raw moments and both Euler-invisible common-label derivatives.
This exposes the missing joint birth-clock moment without new tails, sampling,
DP, a thermal-curve rerun, or an independently counted evidence block.

Data lock: `f4682eb379b5709a2840faf92beef44ff27f6f23`.
Reader: `82b715dc`, `scripts/p334_continuous_center_lifetime_moments.py`.
Input forks: `e32a85939279b8574278024d647b56d2d1485247`; contact marks:
`959a7fa26677c416b874d272f1ba66523fb38f73`. Every compressed input is read once;
its SHA256 is recorded. The policy is exactly `4db356e1`: a shared next label,
joint rank safety, equal joint contact degrees, and the two half-sum/difference
marks built only from the rank-zero loop count.

## The continuous law, including its sampling variance

Write `a=K1`, `b=K2`, `D2=(N+1)(N+2)`, and `D3=D2(N+3)`.
The common uniform order statistics obey

\[
E\tau_1={a\over N+1},\quad E\tau_2={b\over N+1},\quad
E\tau_1^2={a(a+1)\over D2},\quad
E\tau_2^2={b(b+1)\over D2},\quad
E\tau_1\tau_2={a(b+1)\over D2}.
\]

Consequently, for `C=(tau1+tau2)/2`, `W=tau2-tau1`, the saved five moments are

\[
EC={a+b\over2(N+1)},\quad EW={b-a\over N+1},\quad
ECW={b(b+1)-a(a+1)\over2D2},
\]
\[
EC^2={a(a+1)+2a(b+1)+b(b+1)\over4D2},\qquad
EW^2={(b-a)(b-a+1)\over D2}.
\]

For example, the conditional sampling contribution to lifetime variance is

\[
\operatorname{Var}(W\mid a,b)
={ (b-a)(N+1-b+a)\over (N+1)^2(N+2)}.
\]

It would be lost by squaring `(b-a)/(N+1)`. The two order statistics are not
independent Beta draws. Direct double births `a=b` retain `W=0` exactly.
The same pass also saves `E tau1^3=a(a+1)(a+2)/D3`, its second-clock analogue,
and `R1_plateau_M2=E(tau2^3-tau1^3)/3`.

## What is genuinely joint

`ECW=(E tau2^2-E tau1^2)/2` is reconstructible from endpoint marginals.
The additional joint coordinate is

\[
E(\tau_1\tau_2)=EC^2-\tfrac14EW^2,
\]

and, after population centering,
`Cov(tau1,tau2)=Var(C)-Var(W)/4`. This is why both `C2` and `W2` are retained
rather than replacing them by moments of the two separate marginal clocks.
For connected tangents, within each orientation the coordinator must compute

\[
H_{\mathrm{Cov}(C,W)}=H_{CW}-\mu_C H_W-\mu_W H_C,
\quad H_{\mathrm{Var}C}=H_{C^2}-2\mu_C H_C,
\quad H_{\mathrm{Var}W}=H_{W^2}-2\mu_W H_W.
\]

Pooled population means, and their delete-one-batch versions, belong in these
products. A mean of batchwise products is a different estimator. The p158
coordinator handles that shared covariance; this output does not add a second
joint significance analysis. In particular a visible raw `H_C2` cannot alone be
called a change in dispersion, because it includes center translation.

## Direct baseline readout

Each entry below is the original fullpopulation mean of its thirty-two saved
conditional tails per prefix; the uncertainty in parentheses is one batch SE.

| N | orientation | mean C | mean W | E W² |
|---|---|---:|---:|---:|
|325|first|0.59194303 (0.00022423)|0.04771577 (0.00025158)|0.00394726 (0.00004304)|
|325|second|0.59198266 (0.00024857)|0.04776412 (0.00017449)|0.00393758 (0.00002839)|
|425|first|0.59241390 (0.00024195)|0.04324783 (0.00016205)|0.00322576 (0.00002289)|
|425|second|0.59182383 (0.00027277)|0.04389325 (0.00024332)|0.00331287 (0.00003593)|

These are linear descriptive readings, not a new cross-size or orientation
test. The common geometry and suffix dependence remain encoded by the same
twenty original batches for each size.

## Platform and thermal crosswalk

For the rank-one plateau `[tau1,tau2)`, its integrated mass is `EW`, its first
thermal moment is `ECW`, and its second is `R1_plateau_M2`. Thus its
lifetime-weighted centroid and squared width use `ECW/EW` and
`R1_plateau_M2/EW-(ECW/EW)^2`, respectively. They are not the unweighted `EC`
and `Var(C)`. Root handles those ratio readings with the same batch source.
Equivalently, with `A=F1+F2`, `E=F2-F1`, the derivative identities are
`H_C=-integral(H_A)/2`, `H_W=-integral(H_E)`,
`H_CW=-integral(p H_E)`, and `H_R1_plateau_M2=-integral(p^2 H_E)`.

## Machine interface and normalization

`results/p334-continuous-center-lifetime-moments/batch_moments.json` contains
`sizes[N].labels` and `joint_20_batch_means` (20 by 288), plus integer numerator
arrays and the exact moment denominators. Example labels:

- `all.baseline.first.C`, `all.baseline.second.CW`;
- `all.H.plus.first.C2`, `all.H.minus.second.W2`;
- `all.baseline.first.R1_plateau_M2`, `all.H.plus.first.R1_plateau_M2`.

Replace `all` by `00,01,02,10,20` for zero-padded cell contributions. Every cell
uses the full 1000-prefix denominator per batch. The five active cells exhaust
the tangent; **they do not exhaust the unperturbed baseline**, whose `all`
includes all nine rank cells. Baseline integer sums divide by `32000` times
the moment denominator. Tangent sums divide by the original `64000` times the
moment denominator; its mark half-factor must not be applied again.

No raw orientation moments are formed by multiplying paired S/D averages.
S/D projections, with the stored signed `delta_cos4`, follow only after the
orientation-specific population/connected operations.

Reproduction (one input pass; no simulation):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
/Users/lc/python-envs/research-py311/bin/python \
scripts/p334_continuous_center_lifetime_moments.py --output /tmp/p334-continuous-moments-new
```
