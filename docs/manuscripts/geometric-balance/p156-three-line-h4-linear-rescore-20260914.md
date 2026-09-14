# The archived P156 physical three-line H4 is a fixed linear transform of C,Q,S

2026-09-14.  Zero-new-sampling algebra for the frozen #156 pilot.  This note does not create an additional evidence block: it is a deterministic linear view of the same three sector probabilities and the same 100 batch partitions.

## 1. Frozen C3 coordinates

For the three declared primitive-line probabilities `(P_0,P_1,P_2)`, #156 uses

\[
C=P_0-\frac{P_1+P_2}{2},                                      \tag{1.1}
\]

\[
Q=\frac{\sqrt3}{2}(P_2-P_1),                                 \tag{1.2}
\]

\[
S=P_0+P_1+P_2.                                                 \tag{1.3}
\]

These are an invertible real coordinate system:

\[
P_0=\frac{2C+S}{3},                                           \tag{1.4}
\]

\[
P_1+P_2=\frac{2(S-C)}{3},                                     \tag{1.5}
\]

\[
P_2-P_1=\frac{2Q}{\sqrt3}.                                   \tag{1.6}
\]

## 2. Physical embedded spin-four on the same three lines

For `tau=1/2+i y`, let

\[
Z_4(\ell_1)
=\left(\frac{\tau}{|\tau|}\right)^4
=c+i s,                                                        \tag{2.1}
\]

while reflection gives

\[
Z_4(\ell_2)=c-i s,\qquad Z_4(\ell_0)=1.                       \tag{2.2}
\]

The physically embedded three-line harmonic is

\[
A_4^{(3)}=P_0+(c+i s)P_1+(c-i s)P_2.                          \tag{2.3}
\]

Substitute (1.4)--(1.6):

\[
\boxed{
\Re A_4^{(3)}
=\frac{2(1-c)}{3}C
+\frac{1+2c}{3}S,}                                            \tag{2.4}
\]

\[
\boxed{
\Im A_4^{(3)}
=-\frac{2s}{\sqrt3}Q.}                                       \tag{2.5}
\]

Therefore physical three-line H4 is a **fixed declared linear transform of the existing frozen coordinates**.  No reclassification, extra sector lookup or new sample is needed.

The same formulas apply to continuum-subtracted residuals because the transform is linear:

\[
\Delta A_4^{(3)}
=L_\tau(\Delta C,\Delta Q,\Delta S).                          \tag{2.6}
\]

The covariance transforms by the same fixed matrix.

## 3. Hexagonal fixed point recovers C exactly

At the exact Eisenstein modulus,

\[
c=-1/2,
\qquad |s|=\sqrt3/2.                                          \tag{3.1}
\]

Equation (2.4) becomes

\[
\Re A_4^{(3)}=C,                                              \tag{3.2}
\]

while (2.5) is, up to the orientation sign convention,

\[
\Im A_4^{(3)}=\pm Q.                                          \tag{3.3}
\]

Thus the complex C3 character `C+iQ` is literally the embedded physical spin-four three-line readout at the hexagonal point.  Away from that point the linear map changes with the physical modulus.

## 4. N30 and N56 transform coefficients

### N30: `tau=1/2+5i/6`

\[
c=-0.5570934256055362\ldots,
\qquad
s=-0.8304497873350802\ldots.                                  \tag{4.1}
\]

Hence

\[
\Re A_4^{(3)}
=1.038062283737024\ldots\,C
-0.038062283737024\ldots\,S,                                  \tag{4.2}
\]

\[
\Im A_4^{(3)}
=0.958929\ldots\,Q.                                           \tag{4.3}
\]

Using only the rounded residual point estimates printed in the original pilot note (`Delta C=0.00754883`, `Delta Q=-0.00132089`, `Delta S about 0.013`) gives the orientation only:

\[
\Delta A_4^{(3)}\approx 0.00734-0.00127i,                     \tag{4.4}
\]

with the explicit warning that the `S` number in that prose was rounded.  The committed batch CSV should be used for the authoritative value and covariance.

### N56: `tau=1/2+7i/8`

\[
c=-0.4844970414201184\ldots,
\qquad
s=-0.8747927906331950\ldots.                                  \tag{4.5}
\]

Thus

\[
\Re A_4^{(3)}
=0.989664694280079\ldots\,C
+0.010335305719921\ldots\,S,                                  \tag{4.6}
\]

\[
\Im A_4^{(3)}
=1.01013\ldots\,Q.                                            \tag{4.7}
\]

The rounded original residuals (`Delta C=0.00175134`, `Delta Q=0.00112146`, `Delta S about 0.0008`) give roughly

\[
\Delta A_4^{(3)}\approx0.00174+0.00113i.                      \tag{4.8}
\]

Again the batch-rescore script, not these rounded prose inputs, is the authoritative route.

## 5. Covariance and evidence accounting

The archived file

`results/local-20260829/P156-square-bond-primitive-pilot/result.batches.csv`

contains `l0,l1,l2,rank1_other` for each original batch.  Therefore the script

`scripts/rescore_p156_projective_h4.py`

can compute the full `2 x 2` covariance of `(Re Delta A4_three, Im Delta A4_three)` under the **same** batch partition.

This transformed statistic is not independent evidence from `(C,Q,S)`.  It is a derived view of exactly the same three counts and must remain in the same dependency block under the repository governance rules.

## 6. What cannot be recovered

The archive stores only the total `rank1_other` count, not the primitive line inside that category.  Therefore

\[
\boxed{A_4^{full}\text{ cannot be reconstructed from the old pilot}.} \tag{6.1}
\]

No average phase may be imputed to `rank1_other` after reveal.  A future full projective harmonic measurement would have to retain the primitive `(u,v)/+/-` label (or directly accumulate the declared complex harmonic) prospectively.

This information-loss statement is as important as the recoverable three-line transform: it prevents a convenient but invalid post-hoc full-H4 reconstruction.
