# Frozen stronger diagonal-edge extension

Status: `FROZEN_BEFORE_OUTPUT`

The primary Bell/rank gate found a nonzero first Bell transition whose sparse
kernel values were both zero.  This extension keeps every other definition,
root, Schur coefficient, ordering rule, and stop rule unchanged, but further
requires `g16(bell0) != g16(bell1)`.  It tests whether one physical site flip
simultaneously moves rank and the numerical canonical source, rather than only
its labelled Bell state.  It consumes the same existing exact fibres and does
not add a population or inspect a new descriptor.
