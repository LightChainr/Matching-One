# P250 two-sector annihilator bridge freeze

## Object

The model-free Hankel score retained rank five separately in the plus and minus
two-charge blocks.  For each hand, stack the charge-1 and charge-2 degree-two
Hankel matrices into a `12 x 6` matrix and extract its minimum right singular
covector.  This is the truncated quadratic relation

```text
q00 + q10 Tx + q01 Ty + q20 Tx^2 + q11 Tx Ty + q02 Ty^2 = 0.
```

Every delete-one batch re-extracts the line.  Comparisons are projective: a
shared resolved coefficient is normalized to one and the other five complex
coefficients are scored with their complete paired covariance.

## Parameter-free map family

The raw section has already been moved to the exact C4-covariant gauge.  Map
the plus staircase to the minus staircase by either an orientation-preserving
quarter-turn or the Alexander/Gaussian reflection `(a,b)->(a,-b)`, compose
with every `R^k`, and test both linear and coefficient-conjugating transport.

The primary family is reflection plus complex conjugation for `k=0,1,2,3`.
All four are reported; none is selected by its observed p-value.  The family
is rejected at `alpha=0.01` only if every member is rejected.

The deck-generator choices introduce no additional fit.  Relabeling charges
or multiplying charged rows by nonzero fifth-root phases acts invertibly on
the stacked Hankel rows and therefore leaves its common right null line
unchanged.  The exact deck result is invariance, not four duplicate votes.

## Exact limitation

The C4 fiber multipliers `3` and `2` are inverse modulo five, as orientation
reversal requires.  But the same-parent N505 children `(10+i)(2+i)` and
`(10+i)(2-i)` are not D4-equivalent Gaussian quotients: conjugating the first
also conjugates the parent.  A surviving bridge would therefore identify a
truncated Hecke-sector morphism, not a microscopic graph isomorphism.

If every parameter-free bridge fails, the minimal new row is the degree-five
boundary `(5,0)..(0,5)`.  It alone extends every degree-three shift of the
quadratic relation; C4 closure requires the 20-point Manhattan radius-five
shell.  A full next flat-extension matrix additionally needs degree six.
