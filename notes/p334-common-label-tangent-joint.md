# One intervention, two mark directions and two response directions

This readout appends the final common-next-label tangent archive to the
`56b38332` shared twenty-batch covariance. It does not create new prefix or
suffix samples. The source is the original `e32a8593` fork block.

## Fixed response chart

Both orientations must preserve checkpoint rank. Within each joint contact
degree class `(e_first,e_second)`, the next-label intervention preserves class
mass. Define R0-only loop marks `Lf,Ls` and

```
g_plus=(Lf+Ls)/2, g_minus=(Lf-Ls)/2,
S(Y)=(Y_first+Y_second)/2,
D(Y)=(Y_first-Y_second)/delta_cos4.
```

For any chosen birth observer Y, the four entries are the linear responses
`plus->S`, `plus->D`, `minus->S`, `minus->D`. They retain both signs and their
common covariance. The named birth observers are F1/F2 at the original p_ref
and their p-integrals; complete A and E follow by F1+F2 and F2-F1. K1/K2/C/W
retain the same scalar-intervention semantics and the existing linear map.

Let H_ab denote the response of Y_b to mark L_a under this **same joint policy**.
The chart is an invertible linear re-expression:

```
plus ->S = (Hff+Hfs+Hsf+Hss)/4
plus ->D = (Hff-Hfs+Hsf-Hss)/(2 delta)
minus->S = (Hff+Hfs-Hsf-Hss)/4
minus->D = (Hff-Hfs-Hsf+Hss)/(2 delta).
```

A nonzero plus->S is a common clock response. A nonzero plus->D would show a
geometry-differential response to one common intervention; its uncertainty is
measured directly, not inferred from the positivity of individual responses.
The minus source is a separately named antisymmetric intervention, so its D
response does not establish that the plus source distinguishes geometry.

## Relation to the previous own-policy readout

The old contact result fixes only the responding orientation's own R0-safe
domain and its own degree e. Its partner can change rank, and the partner's
degree is not held fixed. The new policy restricts both ranks and preserves
the joint degree-class distribution. Therefore neither the two old positive
clock slopes nor their pooled average determines any entry of the new matrix.
Their shared covariance can be preserved without identifying their estimands.
No subtraction between these different interventions is called an unexplained
physical residual.

## Covariance contract

Consume only a final committed source, with original batch IDs 0..19 and
denominator 1000 prefixes per N/batch. Append its centered common-batch
factor to the existing factor. Linear coordinates and uncertainty are retained
without a high-dimensional inverse, omnibus fit, variance clipping or scanning
for the strongest contrast. If plus->D is weak, report its estimate and SE
directly. The 2x2 chart is a response description, not a field count.

Status: interface and definitions prepared before reading the new final score.
