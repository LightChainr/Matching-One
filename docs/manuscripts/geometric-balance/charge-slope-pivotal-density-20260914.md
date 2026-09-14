# Charge-free-energy slope as a safe-sector pivotal density

2026-09-14.  Exact finite-width identity plus a critical-scaling interpretation.

This note connects the fixed-width charge transfer to the repository's existing pivotal/arm-analysis machinery.  The algebraic identity is exact.  The final `w^-1/4` arm interpretation uses critical scaling and is not promoted as a rigorous square-site exponent theorem.

## 1. Finite strip event and Russo derivative

For graph `G=G4` or `G8`, circumference `w`, and open strip height `m`, let

\[
A_{w,m}^G=\{\text{no horizontally essential occupied component}\},
\qquad
Q_{w,m}^G(p)=P_p(A_{w,m}^G).                                 \tag{1.1}
\]

This is a decreasing event of the `wm` independent site variables.  Russo's formula therefore gives

\[
\boxed{
\frac{d}{dp}Q_{w,m}^G(p)
=-\sum_{v}P_p(v\text{ is pivotal for }A_{w,m}^G).}             \tag{1.2}
\]

Consequently

\[
-\frac1m\frac{d}{dp}\log Q_{w,m}^G(p)
=rac1m\sum_v
\frac{P_p(v\text{ pivotal})}{Q_{w,m}^G(p)}.                  \tag{1.3}
\]

The finite-state Perron representation makes the `m->infinity` limit analytic, so

\[
\boxed{
(I^0_{G,w})'(p)
=\lim_{m\to\infty}
\frac1m\sum_v
\frac{P_p(v\text{ pivotal})}{Q_{w,m}^G(p)}.}                 \tag{1.4}
\]

Thus the derivative of the topological void free energy is literally a pivotal intensity per transfer row, normalized by survival in the safe sector.

## 2. Safe-conditioned form

The pivotal event for site `v` depends only on the other site variables.  Since `A` is decreasing, on a pivotal outside configuration the event `A` holds exactly when `v` is closed.  Therefore

\[
P_p(A\cap\{v\text{ pivotal}\})
=(1-p)P_p(v\text{ pivotal}).                                 \tag{2.1}
\]

Writing `R^G_{w,m}` for the number of pivotal sites in the strip,

\[
\boxed{
-\partial_p\log Q_{w,m}^G(p)
=\frac1{1-p}E_p[R^G_{w,m}\mid A_{w,m}^G].}                   \tag{2.2}
\]

Hence the infinite-strip safe/Q-process has a well-defined pivotal density per row

\[
r^0_{G,w}(p)
:=\lim_{m\to\infty}\frac1mE[R^G_{w,m}\mid A_{w,m}^G],        \tag{2.3}
\]

and

\[
\boxed{(I^0_{G,w})'(p)=\frac{r^0_{G,w}(p)}{1-p}.}             \tag{2.4}
\]

This is the pivotal version of the Perron-eigenvector derivative.

## 3. Equivalence with the safe-row occupation deficit

Let `bar K^0_{G,w}(p)` be the mean number of occupied sites in one added row under the Perron/Doob safe phase.  Differentiating Bernoulli row weights gives independently

\[
(I^0_{G,w})'(p)
=\frac{wp-\bar K^0_{G,w}(p)}{p(1-p)}.                         \tag{3.1}
\]

Equating (2.4) and (3.1) yields the exact identity

\[
\boxed{
wp-\bar K^0_{G,w}(p)=p\,r^0_{G,w}(p).}                       \tag{3.2}
\]

So the conditioned safe phase is depleted relative to an ordinary Bernoulli row by exactly `p` times its topological pivotal density.

This gives a useful semantic check on the left/right Perron vectors: their occupation bias is not an arbitrary spectral statistic; it measures how many sites per row could destroy safe homology if opened.

## 4. Charge slope at coexistence

At black density `p` and complementary matching density `q=1-p`,

\[
\Theta_w(p)=I^0_{4,w}(p)-I^0_{8,w}(q).                       \tag{4.1}
\]

Differentiating in black `p`,

\[
\Theta'_w(p)
=(I^0_{4,w})'(p)+(I^0_{8,w})'(q).                             \tag{4.2}
\]

Using (2.4) separately in the two safe phases,

\[
\boxed{
\Theta'_w(p)
=\frac{r^0_{4,w}(p)}{q}
+\frac{r^0_{8,w}(q)}{p}.}                                    \tag{4.3}
\]

Equivalently, from (3.2),

\[
\boxed{
pq\,\Theta'_w(p)
=p\,r^0_{4,w}(p)+q\,r^0_{8,w}(q)
=w-\bar K^0_{4,w}(p)-\bar K^0_{8,w}(q).}                     \tag{4.4}
\]

At the charge root this is the exact thermal response that divides the critical sector mismatch to produce the pseudo-critical shift.

## 5. Critical arm interpretation

A site that is pivotal for creation of horizontal homology must connect macroscopically distinct occupied/vacant topological channels.  In an ordinary critical bulk picture this is a four-arm-type local event.  The standard percolation four-arm exponent is

\[
x_4=5/4.                                                       \tag{5.1}
\]

A row contains `w` possible pivotal locations, suggesting

\[
r^0_{G,w}(p_c)\asymp w\,\pi_4(w)
\asymp w^{1-5/4}=w^{-1/4}.                                   \tag{5.2}
\]

This is exactly the power observed in the safe-transfer control,

\[
\Theta'_w(p_w^{ch})w^{1/4}
=3.4851,3.4479,3.4262,3.4124,3.4031                         \tag{5.3}
\]

for `w=4,...,8`.

Equation (5.2) should be read as a scaling mechanism, not as a completed square-site proof.  The conditioning on survival in the magnetic/topological safe sector must be shown not to change the bulk thermal exponent; in CFT language this is precisely the statement that differentiating the magnetic-sector energy couples to the ordinary thermal field.

## 6. Practical analysis route

This identity suggests a much cheaper way to interrogate the charge slope than numerical finite differences of eigenvalues.

1. Work in the Perron/Doob safe chain at the charge root.
2. Mark whether each site of the new row is topologically pivotal for safe survival.
3. Measure the two conditional pivotal densities `r^0_4,r^0_8`.
4. Verify (4.3) against the eigenvalue derivative.
5. Study the local multi-arm geometry of those pivotal samples.

The repository already contains several C4/pivotal analysis tools; the missing step is to retarget them to the **safe-sector Q-process** rather than an unrelated unconditioned local source.

## 7. Claim boundary

Equations (1.2)--(4.4) are exact consequences of Russo's formula, Bernoulli score differentiation, and the fixed-width Perron limit.  The identification of the asymptotic power with a four-arm exponent is a critical-scaling hypothesis for square-site percolation, albeit one that matches both the transfer data and the standard thermal exponent `y_t=3/4`.
