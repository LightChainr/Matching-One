# P154 lag-one: an exact birth-posterior clock and a conditional N scale

**New result:** the lag-one source response divided by a birth slope is
exactly `1/N` times a positive-posterior average of the birth-conditioned
source mark. This converts the already measurable rank-population response
into a dimensionless *relative birth-mark* coordinate. It requires no rigid
translation ansatz and does not determine the angular transmission to U.

This note uses only the frozen kernel `4daae57e`, the prior algebra
`6909d9c2`, and the root decision card `d0a9daf1`. No paths, histograms,
new experiment outputs, or old covariance were reprocessed. The fresh
N85/N340 production owned by another team was not inspected or changed.

## 1. The discrete birth/Russo identity

Fix geometry g. Let K_i be birth i=1,2, `b_i(k)=Pr(K_i=k)`, and
`epsilon_j=s_j−E[s_j | j,rank_j,g]`, in the original bulk-count units.
Write `w_{n,k}(p)=Bin(n,p)[k]` and

\[
a_i(k)=\mathbb E[\epsilon_{k-1}\mid K_i=k],\qquad
F_i(p)=\mathbb E\Pr\{\operatorname{Bin}(N,p)\ge K_i\}.
\]

Rank conditioning gives
`E[epsilon_(k−1) 1{K_i<=k−1}]=0`. Consequently the actual lag-one
response, not an auxiliary equilibrium source, is

\[
J_i(p)=\sum_{k=1}^N w_{N,k}(p)b_i(k)a_i(k),\qquad
F_i'(p)=N\sum_{k=1}^N w_{N-1,k-1}(p)b_i(k).
\]

Direct 0→2 events appear in **both** birth laws. For `0<p<1` and
`F_i'>0`, define the probability measure

\[
\nu_{i,p}(k)=\frac{Nw_{N-1,k-1}(p)b_i(k)}{F_i'(p)}.
\]

The elementary identity `w_{N,k}=(Np/k)w_{N−1,k−1}` proves

\[
\boxed{\quad \Lambda_i(p):=N\frac{J_i(p)}{F_i'(p)}
 =\mathbb E_{\nu_{i,p}}\!\left[\frac{Np}{K_i}a_i(K_i)\right].\quad}
\tag{1}
\]

For a pooled geometry mean, replace ν by the probability measure on
`(g,k)` with its original geometry weight. Equation (1) is unchanged,
with J_i and F_i' pooled **before** taking their ratio.

There is also a local pivotal-site form. If `G_i(A_j)` counts vacant
sites that cause birth i at the next update, then, for k=j+1,

\[
b_i(k)=\frac{\mathbb E G_i}{N-j},\quad
b_i(k)a_i(k)=\frac{\mathbb E[\epsilon_jG_i]}{N-j},\quad
a_i(k)=\frac{\sum_r\Pr(r_j=r)\operatorname{Cov}(s_j,G_i\mid j,r,g)}
 {\mathbb E G_i}.
\tag{2}
\]

Thus a_i measures source enrichment among birth-causing configurations,
not the unconditional magnitude or variance of the bulk source.

## 2. What the factor 1/N proves—and what it does not

If `|a_i(k)|<=M` on the relevant posterior support, `K_i>=cN`, and
`E_nu|K_i−Np|<=C N^alpha`, with alpha<1, then exactly

\[
|\Lambda_i-\mathbb E_\nu a_i|
 \le\mathbb E_\nu\!\left[|a_i|\frac{|Np-K_i|}{K_i}\right]
 \le (MC/c)N^{\alpha-1}.
\tag{3}
\]

Posterior tails require the corresponding weighted-integrability bound;
they cannot simply be discarded. Under these explicit assumptions,
`J_i/F_i'=O(1/N)`. The bulk score itself can grow with N, so its conditional
birth enrichment being O(1) is a **mechanism assumption**, not a consequence
of rank centering. If that enrichment grows as N^beta, the clock ratio may
instead scale as N^(beta−1).

A nonzero a_i cannot be literally constant on the complete finite birth
support: at the latest possible birth, every still-unborn configuration
must give birth at the next update, so the rank-centered mark averages to
zero. A nonzero constant-mark model is a critical-window approximation
or a formal auxiliary kernel, not an exact global implementation of this
source. Root's separate constant-mark resolvent retains this distinction.

Even the pointwise clock has a generally nonzero derivative. With the
source mark and conditional centering held fixed as in the contract,

\[
\Lambda_i'=\Lambda_i/p+
 \frac{N}{1-p}\operatorname{Cov}_{\nu_{i,p}}(a_i(K_i)/K_i,K_i),\qquad
J_i'=\frac{\Lambda_iF_i''+\Lambda_i'F_i'}N.
\tag{4}
\]

This identifies precisely the term lost by treating an observed clock
ratio as a constant shift. Neither (1) nor its 1/N condition fixes the
geometry dependence or p-derivative needed for U; no U power law follows
without those extra assumptions.

## 3. A single next-experiment discriminator

At each N's prescribed pooled root let
`D=Fbar1'+Fbar2'`, `c=(Fbar2'−Fbar1')/D`, and let R be the existing
root-comoving rank-1 response. Define just one dimensionless coefficient

\[
\boxed{\Xi_N=\Lambda_1-\Lambda_2
 =N\left(\frac{\overline J_1}{\overline F_1'}-
          \frac{\overline J_2}{\overline F_2'}\right)
 =\frac{2NR}{D(1-c^2)}.}
\tag{5}
\]

To prove the last equality, use
`R=−(Jbar2−Jbar1)+(Jbar1+Jbar2)(Fbar2'−Fbar1')/D`.
This is an identity for arbitrary lag-one kernels. The earlier note used
the same algebraic coefficient as a relative **rigid** shift and required
that hypothesis to predict U. Equation (5) now identifies what the
coefficient actually measures without that hypothesis: the difference
between the two slope-weighted birth-source enrichments.

Also exactly, `−N*rootdot` is the slope-weighted mean of Λ1 and Λ2.
A common nonzero enrichment can therefore move the root while a relative
enrichment is visible in R. Neither automatically gives angular U.

For the forthcoming same-lineage N85/N340 experiment the single proposed
comparison is **Xi340−Xi85**, with the source's bulk normalization
unchanged and original paired-batch omissions propagated through (5).
Zero is the conditional prediction of a scale-stable birth-enrichment
difference; it is not a universal identity. Under the additional project
critical-slope law `D~N^(3/8)` and stable c, it predicts
`R340/R85~4^(−5/8)=0.420448`, not `1/4`. The exact comparison uses the
measured D and c through (5), without fitting a new exponent.

For context only, multiplying the **already published** lag-one
relative-clock planning values by N gives Xi260=−.131115 and
Xi340=−.130116. These are descriptive centers, with no new uncertainty
calculation or parameter fitting, and do not freeze a numerical Xi85
prediction. This auxiliary discriminator needs no new observer or
simulation: J1/J2, D, c and R are existing scorer outputs. It does not
change another team's production contract or promote a weak U response
to a resolved transmission mechanism.

**Scientific boundary:** the new content is the positive birth-posterior
representation, its explicit conditional scale bound, and the exact
non-rigid meaning of Xi. Rigid common-clock cancellation and the
constant-mark resolvent/U correction belong to the separate prior/root
notes. No new samples, tests, fitting, or external comments were used here.
