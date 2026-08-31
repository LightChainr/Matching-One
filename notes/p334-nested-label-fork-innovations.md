# A nested shared-label fork separates first-site information from later clock noise

The next-label experiment has a direct mechanism coordinate:

\[
 \boxed{\Gamma_{\rm next}
 =\tfrac14[B_{AA}-B_{EE}]
 =\operatorname{Cov}_U\{E[F_1|Z,U],E[F_2|Z,U]\mid Z\}}.
\]

It measures whether the next label moves the conditional first/second-birth
responses together or trades one against the other. It is exactly zero when
the first birth is already known, and can be activated by R0 prefixes. This
is a signed clock-coupling readout, not another positive variance component
or a field identification.

## Sampling semantics and the three matrix estimands

Z is the common **ordered** prefix, with both embedded geometries and their
checkpoint ranks. Draw U,V independently and uniformly from its remaining
labels, **allowing U=V**. Under each one-label child draw two fresh uniform
suffix permutations a,b. The four tails are conditionally independent given
Z,U,V. Within each tail, the two orientations use the same label sequence.
Thus a tail returns one complete paired vector X, including both F1/F2 and
their fixed canonical/integrated linear readouts. Cross-branch common random
numbers beyond Z,U,V would change the following estimands.

Set `a=X_Ua-X_Va`, `b=X_Ub-X_Vb`, and record

\[
 \widehat V=(aa^T+bb^T)/4,\quad
 \widehat B=(ab^T+ba^T)/4,\quad
 \widehat W=(a-b)(a-b)^T/4.
\]

There is a **samplewise** identity `Vhat=Bhat+What`. With
`m_U=E[X|Z,U]`, `Q_U=Cov(X|Z,U)`, put

\[
 B(Z)=\operatorname{Cov}_U(m_U|Z),\quad
 W(Z)=E_U(Q_U|Z),\quad V(Z)=\operatorname{Cov}(X|Z)=B(Z)+W(Z).
\]

Conditioning on U,V, the difference means are `delta=m_U-m_V` and their
noise covariance is `Q_U+Q_V`; a,b are conditionally independent. Hence
`E[Bhat|Z,U,V]=delta*delta^T/2`,
`E[What|Z,U,V]=(Q_U+Q_V)/2`, proving respectively
`E[Bhat|Z]=B(Z)`, `E[What|Z]=W(Z)` and `E[Vhat|Z]=V(Z)`.

Vhat and What are positive semidefinite per quartet. Bhat need not be: its
negative sample eigenvalues or scalar projections must not be clipped.
Clipping would bias next-label information and break samplewise closure.
Their expectations, including the full-population average E[B(Z)], are PSD.
Forcing U and V to be different would instead multiply the B term by
`d/(d-1)` and also change Vhat's target; independent draws avoid that change.

## Which physical direction receives the new information?

Use the even topology convention explicitly:
`A=P2-P0=F1+F2-1`, **`E=P2+P0=1-F1+F2`**. The plateau is `1-E`, not E.
The two birth directions are

\[
 (A-E)/2=F_1-1,\qquad (A+E)/2=F_2.
\]

For an R1 or R2 orientation, K1 is prefix-measurable. Its F1 differences
vanish in every tail, so the entire quartet's A/E differences lie on `(1,1)`;
for R2 both births are known and all differences vanish. If both orientations
have rank at least one, this remains exact for their paired difference.
An R0 orientation can supply the first-birth `(1,-1)` direction.

In either an individual orientation or the fixed H4-normalized paired
contrast, the corresponding next-label information coordinates are

\[
 B_{11}=\tfrac14(B_{AA}+B_{EE}-2B_{AE}),\quad
 B_{22}=\tfrac14(B_{AA}+B_{EE}+2B_{AE}),\quad
 B_{12}=\tfrac14(B_{AA}-B_{EE})=\Gamma_{\rm next}.
\]

The last coordinate has an immediate quartet estimator
`Gammahat=(a_A*b_A-a_E*b_E)/8`. It is exactly zero on both-ranks-at-least-one
prefixes, without averaging away Monte Carlo noise. On R0-containing
prefixes its sign tests aligned versus opposing first/completion innovation.
At the expectation level `Gamma^2<=B11*B22`; a single indefinite Bhat is not
subject to that PSD constraint. The same signed projection of What measures
the later-suffix cross-birth covariance, and the two sum to the total one.

For the integrated clock coordinates `C=(K1+K2)/2`, `Wclock=K2-K1`,

\[
 A_{\rm int}=1-2C/(N+1),\qquad E_{\rm int}=1-W_{\rm clock}/(N+1),
\]
\[
 \Gamma_{\rm next}^{\rm int}
 =\frac{4\operatorname{Var}_U E[C|Z,U]
       -\operatorname{Var}_U E[W_{\rm clock}|Z,U]}{4(N+1)^2}
 =\frac{\operatorname{Cov}_U(E[K_1|Z,U],E[K_2|Z,U])}{(N+1)^2}.
\]

For a paired contrast, replace each clock by the same paired normalized
contrast. When K1 is known, `delta Wclock=2 delta C`; the signed coordinate
vanishes. R0 makes their difference an observable physical mechanism rather
than a duplicated completion statistic. W(Z), the matrix of later noise,
is distinct from the scalar lifetime Wclock throughout.

## A finite fork average is not an exact conditional mean

Suppose L independent next-label groups each have t independent tails.
Their mean Xbar has

\[
 E[\bar X|Z]=m(Z),\qquad
 \operatorname{Cov}(\bar X|Z)=B(Z)/L+W(Z)/(Lt).
\]

Thus **32 tails in eight independent quartets** mean L=16,t=2 and give
`B/16+W/32`, not `(B+W)/32`. If only two next labels were reused for16 tails
each, the result would be `B/2+W/32`. The number of independent next-label
draws matters, not the number of distinct labels observed or tails alone.

Let X0 be the original baseline suffix observation and assume the fresh
fork RNG is independent of that suffix conditional on Z. Then

\[
 \operatorname{Cov}(X_0,\bar X)=\operatorname{Cov}(m(Z)).
\]

Both readouts share the original prefix variability and are not independent
replicates. For the eight-quartet design the actual population covariance
reduction is

\[
 \operatorname{Cov}(X_0)-\operatorname{Cov}(\bar X)
 =\tfrac{15}{16}E B+\tfrac{31}{32}E W.
\]

By contrast,

\[
 E[(X_0-\bar X)(X_0-\bar X)^T]
 =\tfrac{17}{16}E B+\tfrac{33}{32}E W.
\]

The baseline-minus-fork residual measures **added** conditional noise, not
removed noise. The old exact-DP residual identity must therefore not be
reused for finite fork averages. The already frozen Bhat/What estimates give
the appropriate gain directly; no regression or new kernel is required.

## Population and scientific handoff

All grouping is measurable from Z, particularly the nine `(Rf,Rs)` cells.
Keep all20000 original prefixes and20 original1000-prefix batches per size.
Average quartets within a prefix before forming those batch means;32 tails
do not become32 new prefix replicates. Known both-R2 prefixes contribute
their exact full response and zero fork-noise matrices, not missing rows.
Any prefix-measurable restriction must retain its full-population weight.

The most direct new output is the R0-supported first/completion innovation
plane `(B11,B22,Gamma)`, alongside later-suffix W and their total closure,
for both canonical and integrated paired responses. Report the same-batch
covariance of these named projections and the mean readouts; no PSD clipping
or high-dimensional inverse is part of the estimand.

Scientific card: this defines a genuine shared-label intervention separating
next-site information from suffix noise and yields a new signed cross-birth
coordinate with an exact R1/R2 control. It also gives the correct variance
gain of finite nested averaging. These are finite-permutation identities,
not numerical findings, independent data, field identification or new DP.
Production and all joint covariance remain with the root and its single
coordinator. The source population and full A/E semantics stay unchanged.
