# Intrinsic birth covariance response: prefix, label, and suffix mechanisms

A nonzero total intrinsic rank-covariance response can come entirely from
different original prefixes moving by different constant clock shifts.
The **within-prefix covariance derivative** removes that explanation.
It can then be separated into selection of child means and selection of
child suffix covariances, using the existing two independent suffixes per
next label. All products needed for the separation have unbiased finite-M
estimators with the eight independent quartets already collected.

## Fixed population, one physical orientation, exact census score

Set `X=K1/(N+1)` and `Y=K2/(N+1)`. Their covariance is exactly the
`Cov(K1,K2)/(N+1)^2` intrinsic coordinate in e2ef9983; no uniform
order-statistic covariance is included here. Work in one orientation before
forming its S/D response with the other orientation.

Z is the original **complete ordered prefix**, including both geometries
and their ranks. Its population law is unchanged by the policy. The
common next label U is uniform at t=0. Its score for either frozen source
mark is supplied by the complete label census:

\[
 s(Z,u)=\pi_a\{g(u)-\bar g_a\},\quad u\in A_a,
 \qquad s=0\text{ otherwise},\qquad E_U[s\mid Z]=0.
\]

Here A_a is the jointly rank-preserving joint-degree class of the common
policy, and pi_a and bar g_a are exact full-label values. After U the
uniform suffix law is unchanged. Thus, for any output f,

\[
 h_f(Z):=H E[f\mid Z]=E_U[s(Z,U)m_f(Z,U)\mid Z],
 \qquad m_f(Z,u)=E[f\mid Z,u].
\]

This round uses the **exact-score difference estimator**, with no
same-class pair mask. For independent uniform U,V, allowing U=V, and
two conditionally independent suffixes 0,1 per label, put

```
f_U = (f_U0+f_U1)/2,       f_V = (f_V0+f_V1)/2,
b_q,f = (f_U+f_V)/2,
h_q,f = (s_U-s_V)(f_U-f_V)/2.
```

Then `E[b_q,f|Z]=mu_f(Z)` and `E[h_q,f|Z]=h_f(Z)`. Indeed, the difference
expectation is `Cov_U(s,m_f)=E_U(s m_f)` because the score is centered.
Labels from different classes, or one active and one inactive label,
must not be removed from this exact-score estimator. This is the frozen
estimator in the current census-based production, not the earlier
same-class mark-difference estimator.

## First split: what changes within a fixed original prefix?

Write `mu_X(Z)=E[X|Z]`, `mu_Y(Z)=E[Y|Z]`, and
`C_Z=Cov(X,Y|Z)`. Since the Z law is fixed,

\[
 \boxed{H\operatorname{Cov}(X,Y)
 =E_Z[H C_Z]
  +\operatorname{Cov}_Z(h_X,\mu_Y)
  +\operatorname{Cov}_Z(\mu_X,h_Y).}
\]

The within-prefix derivative is

\[
 H C_Z=h_{XY}-h_X\mu_Y-\mu_Xh_Y.
\]

The other two terms describe how conditional mean responses align with
the original conditional means across prefixes. There is no derivative
of the prefix prevalence in this decomposition.

### Eight-quartet products, with no self-pair bias

Let M=8 independent quartets conditional on the same Z. For any two
quartet vectors a_q,b_q, define the ordered distinct-quartet product

\[
 \mathcal P_M(a,b)
 ={(\sum_q a_q)(\sum_q b_q)-\sum_q a_qb_q\over M(M-1)}.
\]

It is unbiased for `E[a_q|Z] E[b_q|Z]`, even when a_q and b_q within the
same quartet are dependent. Equivalently it is
`bar a bar b - sampleCov_q(a,b)/M`, where sampleCov uses denominator M-1.
Accordingly the finite-M unbiased within estimator is

\[
 \widehat{H C_Z}=\bar h_{XY}
       -\mathcal P_M(b_X,h_Y)-\mathcal P_M(b_Y,h_X).
\]

Multiplying the two quartet averages without the diagonal subtraction
retains shared-tail/label noise. Multiplying b and h within one quartet
does not estimate a product of conditional expectations.

For L original independent prefix records, use the analogous distinct-
prefix product `P_L` on their quartet means. Put

```
W_hat = average_i [hbar_i,XY - P_M(b_i,X,h_i,Y) - P_M(b_i,Y,h_i,X)],
B_XY_hat = average_i P_M(b_i,X,h_i,Y) - P_L(bbar_X,hbar_Y),
B_YX_hat = average_i P_M(b_i,Y,h_i,X) - P_L(bbar_Y,hbar_X),
T_hat = average_i hbar_i,XY
        - P_L(bbar_X,hbar_Y) - P_L(bbar_Y,hbar_X).
```

These estimate the source-population within, two between, and total
responses, respectively, and have the **exact samplewise closure**
`T_hat=W_hat+B_XY_hat+B_YX_hat`. The within self-products are corrected at
the quartet level; the global products are corrected at the prefix level.
This is the source-population U-product convention in f34bcd6f. A different
target conditioned on the exact empirical prefix collection would instead
retain its diagonal population terms with their quartet-noise correction;
the two finite-L targets should not be silently mixed.

Keep the original full20k denominator and20 original batches per size.
Compute P_M separately within each prefix before batch pooling. In each
delete-one-batch replicate retain those prefix quantities, re-pool the
remaining prefixes and re-form P_L with the new L. Quartets, labels and
suffixes do not become independent population replicates.

## Which translation mechanisms does within response rule out?

Consider the conditional moment-level translation null
`X_t=X+t a(Z), Y_t=Y+t b(Z)`. Each birth can have its own arbitrary
prefix-dependent constant shift. Every conditional covariance and variance
is unchanged, so **H C_Z=0 for every Z**. A resolved nonzero E_Z[H C_Z]
therefore rules out this entire class, not merely one global translation.
A zero average does not prove the null; different conditional responses
could cancel.

For the special center shift a=b=v(Z), h_X=h_Y=v(Z), and the surviving
total covariance response is precisely

\[
 H\operatorname{Cov}(X,Y)
 =2\operatorname{Cov}_Z\!\left(v(Z),{\mu_X(Z)+\mu_Y(Z)\over2}\right).
\]

By contrast, a label- or suffix-dependent common shift
`X_t=X+t v(path), Y_t=Y+t v(path)` can change within covariance while
leaving every path's lifetime Y-X unchanged:

\[
 H C_Z=\operatorname{Cov}(v,X\mid Z)+\operatorname{Cov}(v,Y\mid Z)
      =2\operatorname{Cov}\!\left(v,{X+Y\over2}\mid Z\right).
\]

Consequently within covariance alone does not exclude all heterogeneous
center transport. A useful companion already available from XX/XY/YY is

\[
 H\operatorname{Var}(Y-X\mid Z)
 =H\operatorname{Var}(X\mid Z)+H\operatorname{Var}(Y\mid Z)-2H C_Z.
\]

It vanishes for *every* pathwise common center shift, even a suffix-dependent
one; a resolved nonzero value would rule out that broader lifetime-fixed
class. The same quartet-product formulas give the two diagonal variance
derivatives by using `bar h_XX-2 P_M(b_X,h_X)` and its Y analogue. These
are necessary moment predictions, not assertions about the current data.

## Second split: label means versus suffix covariance selection

At a fixed prefix write
`m_X(u)=E[X|Z,u]`, `m_Y(u)=E[Y|Z,u]`, and
`v_XY(u)=Cov_suffix(X,Y|Z,u)`. Conditional total covariance gives

\[
 C_Z=E_U[v_{XY}(U)]+\operatorname{Cov}_U(m_X(U),m_Y(U)).
\]

The suffix law at a specified Z,u is not perturbed, so differentiation yields

\[
 H C_Z=\underbrace{E_U[s\,v_{XY}]}_{H_{\rm suffix\ selection}}
  +\underbrace{E_U[s\,m_Xm_Y]-h_X\mu_Y-\mu_Xh_Y}
                _{H_{\rm label\ means}}.
\]

In particular, the first term selects labels with different *existing*
suffix covariances. It is not a change of the suffix dynamics conditional
on the exact child. Neither term is generically positive: both are signed
responses of covariance quantities.

The two saved independent suffixes provide the required products at each
label, without estimating a noisy child mean and squaring it:

```
raw_u,XY   = (x_u0*y_u0+x_u1*y_u1)/2,
cross_u,XY = (x_u0*y_u1+x_u1*y_u0)/2,
suffix_u,XY = raw_u,XY-cross_u,XY
            = (x_u0-x_u1)(y_u0-y_u1)/2.
```

Their conditional expectations are `E[XY|Z,u]`, `m_X(u)m_Y(u)` and
`v_XY(u)`. Apply the same b/h exact-score formula to each child quantity.
There is a **samplewise** identity `h_raw,XY=h_cross,XY+h_suffix,XY`.
Thus the new within estimators are

```
suffix_hat(Z) = hbar_raw,XY-hbar_cross,XY,
label_hat(Z)  = hbar_cross,XY
                - P_M(b_X,h_Y)-P_M(b_Y,h_X),
within_hat(Z) = suffix_hat(Z)+label_hat(Z).
```

Together with the two between-prefix pieces they close exactly to T_hat,
before and after the original-batch aggregation. Crucially, the mixed
suffix product uses two replicas of the **same next label**. Cross-quartet
products estimate mu_X mu_Y over the whole prefix and cannot substitute
for the child product m_X(u)m_Y(u).

## Minimal once-extracted coordinates and handoff

The current per-prefix/per-quartet baseline
`b=(X,Y,XX,XY,YY)` and two-mark h vectors suffice for total/within/between.
For the label/suffix split, retain additionally b/h for

```
cross_XX=x0*x1,
cross_XY=(x0*y1+x1*y0)/2,
cross_YY=y0*y1,
```

averaged or score-differenced across U,V as defined above. Only h_cross_XY
is strictly needed for the covariance derivative; the requested full
three-entry cross block also supplies baseline and diagonal mechanisms.
The suffix block is raw minus cross and need not be separately stored.
No higher suffix replication or new label draw is necessary.

Everything is evaluated within a physical orientation first; only then
form S=(first+second)/2 or D=(first-second)/delta_cos4 of the **complete
covariance derivatives**. Covariance of paired X/Y differences is a different
quantity with cross-geometry terms. All source marks and outputs retain
the same20-batch covariance, not independent-evidence labels.

Scientific card: the new hierarchy distinguishes prefix mean organization,
next-label conditional-mean organization and selection of suffix covariance,
and specifies which lifetime-preserving transport classes each coordinate
can exclude. It refines the intrinsic-rank result e2ef9983, after the
continuous-timing separation in06abeeae. Exact census scores and the
current f34bcd6f source-population U-product convention are used throughout.
This note supplies formulas and minimal coordinates only: no replay, new
MC/DP, empirical result, retrospective validation campaign, PR or comment.
