# Full-root one-defect gate: a single physical flip moves rank and source

## Decision

The frozen existing-fibre gate returns

```text
TWO_INDEPENDENT_DEFECT_GAIN_REJECTED
```

No new population was generated.  The calculation consumes the complete N25
axis/tilted radius-one collar fibres, imports the pooled matching root and all
componentwise Schur coefficients from the full population, and never
re-estimates a counterterm inside an edge class.

The broad Bell/rank gate contains 135,253 exact raw row classes and 8,250,462
physical pair fibres in which the same alternating site flip changes both
digital rank and canonical source Bell state.  Its first lexicographic Bell
transition has exact nonzero Schur weight `+9.6599270326e-10`.  Both sparse
kernel values happen to be zero, so this first result proves labelled Bell
motion but is not the strongest numerical-source witness.

## Stronger kernel-changing certificate

The extension was frozen before its output and additionally requires
`g16(bell0) != g16(bell1)`.  It retains 6,846 raw classes and 740,950 physical
pair fibres.  The first lexicographic nonzero transition is

```text
geometry          axis
rank              0 -> 1
source orbit       axial2
contact mask       3  (both local black arms contacted)
corner word        7
Bell               9021064 -> 2430024
g16                8 -> 0
P0                 3.2977760475e-5
P1                 4.7982147558e-5
S0                -1.5642115001e-7
S1                +1.9419147797e-9
S0+S1             -1.5447923523e-7
```

The outward rational interval for `S0+S1` is strictly negative.  Thus one
physical thermal defect simultaneously moves the topological slow variable
and the numerical canonical source.  The proposed decomposition into two
independent defects, and the automatic six-arm gain inferred from it, is
false on the exact finite state graph.

## The surviving signed functional is sharply localized

Summing every kernel-changing diagonal edge gives total signed mass

```text
-4.9488399165e-6.
```

Its source-orbit decomposition is

| source orbit | signed mass | share of total signed mass |
|---|---:|---:|
| NN other | `-3.7052834188e-6` | `74.87%` |
| diagonal | `-9.9747361552e-7` | `20.16%` |
| far | `-1.5234702796e-7` | `3.08%` |
| knight B | `-1.2215786505e-7` | `2.47%` |
| knight A | `-6.4285415345e-8` | `1.30%` |
| axial2 | `+9.2707426182e-8` | `-1.87%` |

The rank-stage decomposition is even simpler:

```text
first birth  0->1   -5.8210906659e-6   117.63% of final total
second birth 1->2   +8.7225074943e-7   -17.63% cancellation
```

Every one of the 12 rank-transition by source-orbit cells excludes zero, as
do both `S*1` row sums and all six `1^T*S` column sums.  The surviving leading
functional is therefore not a sparse accidental edge.  It is a contact-local,
first-birth-dominated four-arm channel with a smaller completion cancellation.
This is the concrete mechanism that replaces the retired independent-defect
picture.

## Exact contact support rule

The pre-frozen four-mask decomposition gives an even sharper structural
result.  Among all 6,846 kernel-changing diagonal row classes, the no-contact
mask is not merely cancelled: it is absent.

| local contact mask | exact row classes | pair fibres | signed mass |
|---:|---:|---:|---:|
| 0, no arm | 0 | 0 | exactly `0` |
| 1, first arm | 940 | 117,870 | `-4.0937024059e-6` |
| 2, second arm | 952 | 119,663 | `-4.1126658895e-6` |
| 3, both arms | 4,954 | 503,417 | `+3.2575283789e-6` |

Thus, on both complete N25 geometries under the radius-one collar semantics,

```text
kernel-changing diagonal edge  =>  local source-to-thermal-arm contact.
```

The two one-arm channels agree to about `0.46%`; together they produce
`-8.20637e-6`, while the both-arm sector cancels about `39.7%` of that
magnitude.  This is a local collision/OPE selection rule, not evidence for a
separated two-defect object.  The frozen decision string for an empty mask-0
residual is `CONTACT_ZERO_RESIDUAL_UNRESOLVED`; the post-output scientific
content is stronger and simpler: there is no radius-one contact-free edge to
bound.

## Contact fusion and birth stage do not factorize

The exact follow-up tensor keeps contact mask, rank transition and source orbit
jointly, rather than comparing their marginal sums.  It reveals a sharper
support rule:

```text
one-arm contact (masks 1 or 2)  =>  NN source orbit only;
every non-NN source orbit        =>  both-arm contact (mask 3).
```

This is literal support, not cancellation.  The two one-arm masks contain
1,892 raw classes and 237,533 physical fibres, all in `nn_other`.  Mask 3
contains all six source orbits.

Collapse masks 1 and 2 to the exchange-even single-contact column and retain
mask 3 as the fused double-contact column.  The signed full-root Schur tensor
is

\[
\begin{array}{c|rr}
 & \text{single contact} & \text{double contact}\\ \hline
0\mathbin{\to}1 & -2.8838028012142906\,10^{-6}
                 & -2.9372878646696404\,10^{-6}\\
1\mathbin{\to}2 & -5.3225654941777728\,10^{-6}
                 & +6.1948162436108913\,10^{-6}
\end{array}
\]

Its determinant has the strict outward enclosure

```text
-3.3498535471290615e-11 < det < -3.3498535471290614e-11.
```

The normalized sign determinant is exactly `-1` at displayed precision:
the two determinant products have opposite signs.  Thus contact fusion and
birth stage are maximally non-factorizing at N25.  At first birth the fused
and single-contact masses agree within `1.85%` and have the same negative
sign.  At second birth the fused sector reverses sign and is `116.39%` of the
single-contact magnitude, turning the completion sum positive.

The sign reversal is not explained by exposure frequency.  After division by
the positive state mass, the single-contact conditional Schur densities are
approximately `-2.00e-4` and `-8.13e-4` at first and second birth, whereas the
double-contact densities are `-8.98e-5` and `+2.147e-4`.  The conditional
completion amplitude itself changes sign.

Within mask 3, the NN orbit is positive at both stages.  The combined non-NN
orbits rotate from `-3.8159420314e-6` at first birth to
`+2.5723855338e-6` at completion.  The surviving microscopic transmission map
therefore has two coupled operations:

```text
contact fusion  x  topological completion,
```

not one scalar contact counterterm multiplied by one scalar birth amplitude.
The full tensor, positive masses and exact support counts are retained in
[`contact-stage-tensor.json`](contact-stage-tensor.json).

## Consequence

The next theoretical object is no longer a generic six-arm upper bound.  It is
the signed OPE/contact counterterm carried by NN and diagonal source orbits and
its transmission into full original-`U`.  Only after enlarging the collar can
one ask whether a genuinely separated residual is created; at radius one that
residual is exactly empty.

## Boundary

The old raw fibres retain rank, radius-one `C/B/W`, contact word and Bell
before/after, but not a canonical joint `x+y+z` component map or a single
background mask.  The nonzero transition-class sum rigorously proves that at
least one physical edge exists, but does not display one configuration-level
witness or its annular placement.  A future formal joint-map certificate can
make the witness constructive; it is not needed to preserve the mechanism
elimination above.
