# Round-2 work record: attempted closure of r_linear(w=8) and what it revealed

Small honest record.  Attempting to *close* the open numerical item
"exact `r_linear` at w = 8" (#593 gives only `>= 150`) led to a definitional
finding that is worth writing down before anyone quotes that table.

## The attempt

The natural continuous-time object behind "linear realization dimension" of a
task (sources `mu_a`, readouts `f_j`) is the rank of the block-Hankel matrix
built from the Markov parameters `M_k = U^T G^k V` (`U` columns = sources,
`V` columns = readouts):

```
r_Hankel = rank [ M_{i+j} ]_{i,j >= 0}  =  McMillan degree of U^T(sI-G)^{-1}V.
```

I computed this exactly at widths 4–8 (float SVD with threshold; n up to 1430,
blocks 32-wide) — this is not the object #593 counted:

| w | block-Hankel rank (raw / centred) | reported r_linear (#593) |
|---|-----------------------------------|--------------------------|
| 4 | 9 / 8 | 10 |
| 5 | 4 / 5 | 26 |
| 6 | 2 / 2 | 72 |
| 7 | (unstable power iteration) | 218 |

The reported sequence 10, 26, 72, 218 is also *not* reproduced by the raw
(no centring) version, nor by any single-source/readout subset tried in the
quick scan.  Two conclusions:

1. `r_linear` in #593 is not the task block-Hankel/McMillan degree of the
   (4 sources × 8 readouts) impulse response — it is some other linear
   object.  Candidates that would need the #593 script to disambiguate:
   rank of the joint Krylov span of `{G^k f_j}` and `{G^k mu_a}` inside the
   function space; rank of an exact-lumping-type linear (not positive)
   refinement; or the dimension of the smallest *linear* lumping carrying the
   declared dictionary.  Note w4–7 satisfy `r_linear <= r_pos` with equality at
   w4, w5 and strict inequality at w6, w7 (72 < 76, 218 < 232) — consistent
   with *some* quotient/refinement object, not with a full Hankel rank.
2. The atlas item "w8 exact r_linear" therefore stays **open, blocked on a
   definition**, not on compute.  This is recorded rather than papered over:
   quoting `>= 150` or the ratio `125@w8` without the defining matrix would
   invite the same ambiguity.

## Where this leaves the round-2 claims

* The finite-group theorem (main note) does not depend on `r_linear` at all;
  it is a statement about response selection and isotypic factorization, both
  of which are verified exactly.
* Recommended cheap next action for the repository: publish the exact matrix
  whose rank #593 reports as `r_linear` (or its `mod p` elimination input);
  one page of code settles w8 and makes the w4–7 column reproducible.

Boundary: this note makes no claim about which notion is "right"; it only
insists that the name `r_linear` currently denotes a quantity whose definition
is not recoverable from the numbers alone.
