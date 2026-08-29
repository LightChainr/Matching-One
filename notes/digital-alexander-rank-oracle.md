# Digital Alexander rank identity: finite exact oracle

Status: Phase A finite theorem search for Issue 269.  The oracle is exact on
the declared quotients and deliberately stops short of a general Alexander-
duality proof.

## Algebraic reduction

Let `r_b,r_w` be the ambient torus-homology ranks of the black-primal and
white-matching complements.  Since each rank lies in `{0,1,2}`, define

```text
q_either = 1[r_b>0] - 1[r_w>0],
q_cross  = 1[r_b=2] - 1[r_w=2].
```

Exhausting the nine abstract rank pairs proves the elementary implication

```text
q_either = q_cross  =>  2 q_either = r_b-r_w.
```

The only rank pairs satisfying the premise are

```text
(0,0), (0,2), (1,1), (2,0), (2,2).
```

Thus the genuinely geometric step is the already observed equality of the
`either` and `cross` matching differences.  Once that premise holds, the
factor-two identity is an exhaustive rank lemma rather than a new conjectural
fit.

## Finite checks

The additive oracle reuses the canonical general-period homology engine and
does not depend on the open configuration-Betti PR.  It enumerates every mask
on axis L=2/3, Gaussian `(2,1)`, diamond L=2, and the self-matching C4 N=10
control.  For every mask it archives

```text
(r_black, r_white, q_either, q_cross),
q_either-q_cross,
2q_either-(r_black-r_white),
r_black+r_white-2.
```

The last residual tests the stronger rank-sum proposal separately; it is not
silently folded into the weak identity.

The exact joint tables contain only the three rank pairs

```text
(0,2), (1,1), (2,0).
```

Across the five geometries this gives zero failures for both
`2q=r_black-r_white` and the stronger `r_black+r_white=2`.  It also resolves
the older raw-rank-equality counts: `q=r_black-r_white` occurs exactly on the
`q=0`, `(1,1)` configurations, never for `q=+/-1`.

## Deterministic general-period search

A fixed-seed scan samples configurations on 160 nonsingular integer period
matrices with orders 5 through 31.  Empty and full configurations are always
included.  This is a reproducible counterexample search, not a confidence
statement: a null scan cannot replace a proof on arbitrary periodic digital
complexes.

The frozen run evaluated 85,152 configurations and again saw only `(0,2)`,
`(1,1)`, `(2,0)`, with zero weak-identity, common-channel, or strong-rank-sum
counterexamples.

## Boundary and next proof obligation

If the finite oracle remains null, Phase B must show why the 4-connected
foreground and 8-connected matching complement have equal rank-one indicators
under periodic digital Alexander/relative homology.  The proof must handle
nonseparating cycles and distinguish the weak identity from the stronger
rank-sum statement.

No continuum field, amplitude, or value of the square-site threshold follows
from this finite oracle alone.
