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

The chart above was committed as `03d4f58c` before the final score was read.

## Actual shared-policy response

Final source `4db356e1b026853468f94d59d938895a2367ceb7` has now been joined.
The near-reference response matrix is source-selective: plus primarily changes
the common A/E response, while minus gives a clear A difference. The plus->D
entries remain unresolved. Table entries are in units of 1e-5, with original
twenty-batch SE, and share their full covariance.

| N / observer | plus -> S | plus -> D | minus -> S | minus -> D |
|---|---:|---:|---:|---:|
| 325 / A | −5.416 ± .772 | 3.676 ± 2.698 | −.776 ± 1.110 | 15.366 ± 3.467 |
| 425 / A | −5.873 ± 1.266 | −.511 ± 1.901 | −.315 ± 1.343 | 14.839 ± 2.421 |
| 325 / E | 1.884 ± .582 | −.064 ± 1.357 | .580 ± .601 | −1.937 ± 1.845 |
| 425 / E | 1.933 ± .592 | −1.261 ± 1.781 | 1.720 ± .694 | −5.229 ± 1.568 |

The two common clock responses are observable under one intervention. There is
no resolved source-even geometry difference at p_ref in this block. A strong
source-odd A difference is compatible with two primarily local geometry
responses; it does not identify an unperturbed H4 field.

## Direct cross-orientation readout

The original-scale cross responses can be recovered without new statistics:

```
first_loop -> second = R_S+ + R_S- - delta*(R_D+ + R_D-)/2
second_loop -> first = R_S+ - R_S- + delta*(R_D+ - R_D-)/2.
```

F1 cross responses have support only in checkpoint cell 00. In 01+10 one
geometry's first birth is already fixed, while that geometry's R0-only input
mark is zero. Thus the F1 cross response in this group is an exact zero,
not a negative statistical result. F2 can carry a cross response: a safe R0
mark could affect the partner's future R1 completion.

The measured 01+10 canonical F2 cross responses are

| Direction | N325 | N425 |
|---|---:|---:|
| first loop -> second | 1.56587e-5 ± 1.13867e-5 | 1.28997e-5 ± 1.13644e-5 |
| second loop -> first | −9.05532e-7 ± 1.05733e-5 | −3.13374e-6 ± 1.37113e-5 |

Neither is resolved. The predefined 00 cross responses are likewise weak;
all their canonical and integrated values are retained in the score. Absence
of a resolved cross effect is not an exact decoupling theorem.

If these cross responses were zero, even with unequal local gains, the chart
would obey `R_D-=2 R_S+/delta` and `R_D+=2 R_S-/delta`. This makes explicit why
large minus->D can arise without cross-orientation transmission. The new common
policy permits this local-response interpretation directly; the older two own
policies alone could not establish it.

## Two ensemble mean-response directions, one covariance

The parent's already-computed rank readout `73608ba9d3eef34c6980cb5a049f726cfebdd72d`
is appended as supplied LOO columns. Its all-population A Jacobian determinants
are `1.22719e-8 ± 4.07190e-9` / `1.55912e-8 ± 4.69216e-9` at p_ref and
`4.20205e-10 ± 1.62825e-10` / `5.14308e-10 ± 1.28723e-10` after integration.
They are the nonlinear restatement of this same measured response matrix,
with no new data or independent confirmation. In physical orientation
coordinates the matrix is near diagonal, consistent with the weak cross terms.

This gives evidence against one fixed ensemble mean-response direction for
both controls. It is not a per-prefix or per-stratum rank theorem: mixing
different single-direction prefix responses can itself produce ensemble rank
two, as can the 01+10 mixture of two unidirectional controls. General
prefix-dependent scalar states and continuum field counts remain outside this
readout. The rank calculation is post-reveal, as documented by its producer.

Reproduce the thin final join with
`scripts/p334_common_label_tangent_joint.py` in the managed research Python
environment. Outputs are under `results/p334-common-label-tangent-joint/`.
They include all original common-label batch columns, the sixteen named cross
responses, the supplied rank LOO columns and the prior trigger/contact
covariance factor. No raw paths were replayed and no tests, fits or sampling
were added.

Scientific card: instantaneous rank/Euler-preserving controls have a measured
source-even common response and source-odd geometry-difference response.
The ensemble 2x2 response is consistent with primarily local gains; genuine
cross-orientation transmission is not yet resolved. Source/geometry/dependency
remain the same N325/N425 `e32a8593` block, with one common-batch coordinator.
