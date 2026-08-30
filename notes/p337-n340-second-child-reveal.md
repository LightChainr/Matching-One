# P337 N340 same-lineage second-child reveal

Preregistration `5369c21` froze three H4-amplitude targets before any N340 data:
nominal area decay, the already observed N85-to-N170 effective transfer, and
scale-neutral.  It also froze the projective scalar at zero.

Huawei `DevEnvC_HZsCM6` generated 12M samples per shape in 80 aligned batches
with a new counter block.  Runtime was 51.50 seconds; stderr was empty.  Every
metadata freeze gate and the exact Smith `(2,170)` contract passed.

The heldout amplitude is `-0.00485726 +/- 0.00124889`.  The closest frozen target
is nominal H4 (`-1.005` measurement SE); scale-neutral is excluded by `5.008`
measurement SE (`3.946` predictive SE after N170 source uncertainty).  The
secondary effective-transfer target is `2.270` measurement SE away but only
`0.942` predictive SE after retaining its large N85 source uncertainty.

The pair direction flips back positive as exact H4 requires, at `3.889 sigma`
rather than the projected 5 sigma.  The projective scalar remains null at
`0.506 sigma`.  This is evidence for scale curvature returning toward nominal
H4 decay, not evidence for a persistent scale-neutral charged state.  It does
not select a unique correction law or fit a new exponent.
