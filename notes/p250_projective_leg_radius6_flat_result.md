# P250 radius-six flat-extension result

## Decision

The five-state flat-extension hypothesis fails decisively in both hands.  The
first compatible rank in the frozen 5--9 ladder is rank eight for both plus and
R2-gauged minus.

| hand | rank null | chi-square / resolved modes | finite-batch p | decision |
|---|---|---:|---:|---|
| plus | rank <= 5 | `2753.86 / 150` | `1.88e-61` | reject |
| plus | rank <= 6 | `683.42 / 112` | `2.88e-24` | reject |
| plus | rank <= 7 | `395.95 / 78` | `2.58e-19` | reject |
| plus | rank <= 8 | `64.46 / 48` | `0.1978` | do not reject |
| plus | rank <= 9 | `13.18 / 22` | `0.9435` | do not reject |
| minus, R2 gauge | rank <= 5 | `2045.02 / 150` | `5.43e-49` | reject |
| minus, R2 gauge | rank <= 6 | `970.21 / 112` | `3.03e-36` | reject |
| minus, R2 gauge | rank <= 7 | `250.46 / 78` | `2.55e-9` | reject |
| minus, R2 gauge | rank <= 8 | `68.86 / 48` | `0.1214` | do not reject |
| minus, R2 gauge | rank <= 9 | `17.50 / 22` | `0.7825` | do not reject |

The frozen scorer therefore returns `rank5_flat_extension_rejected` and a
truncated rank lower bound of eight in each hand.

## Scientific change

Radius five selected `Alexander R2 + conjugation` as the only surviving map of
the degree-two annihilator line.  Radius six now shows why that result had to
remain a line statement: the line does not generate a five-state flat quotient
when the complete degree-three Hankel block is exposed.

This is not a weak boundary effect.  Ranks five, six, and seven are rejected in
both sectors by many orders of magnitude, while rank eight is comfortably
compatible in both.  Any path-independent realization of these observed
two-charge `H3` blocks consequently needs at least eight states under the
frozen alpha `0.01` rule.

The result does not say that the true dimension equals eight.  `rank(H3)<=8`
surviving is compatibility at this truncation; flatness of an eight-state
system would require stability at the next moment order.

## R2 bridge remains locked

The preregistration allowed the full R2-conjugate kernel-projector comparison
only after both rank-five flatness gates passed.  They did not.  The scorer
therefore emitted `LOCKED_RANK5_FLAT_EXTENSION_FAILED` and did not reinterpret
the two five-dimensional kernels.

A new, separately frozen existing-data analysis could compare the rank-eight
truncated relation spaces in the already selected R2 gauge.  That is a new
question; it must not be smuggled into this reveal as a replacement bridge.

## Execution

HZsCM6 produced 1,200,000 fresh replicas in 400 batches, seed
`25060610120261250`, counters `[0,1200000)`, in 788 seconds.  Exit was zero,
production stderr was empty, and remote/local hashes match.  Raw artifacts were
committed in `f047606` before the one completed frozen score.

An earlier zero-second invocation was rejected by the manifest because it used
the runner default sample count.  It generated no scientific samples and is
preserved only as protocol evidence.

## Scientific card

- Mechanism space changed: the five-state flat quotient is removed; both
  truncated sectors require at least eight states.
- Not proved: exact dimension eight, rank-eight flatness, the R2 full-kernel
  bridge, or ordered translation commutation.
- Observer/sector/source/geometry: fixed-p projective-leg Z5 two-charge endpoint
  moments; plus and R2-gauged minus; minimal degree-six blocks on norm 505.
- Dependency group: 80k degree-four, independent radius-five 1.2M, and fresh
  radius-six 1.2M streams.
- Next lift: freeze an existing-data rank-eight projector bridge if cross-hand
  identity is the target; acquire the next moment order only if exact
  rank-eight flatness is worth distinguishing.
