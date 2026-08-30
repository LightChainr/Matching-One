# P333/P321/P370 one-mark endpoint-jet gate

The terminal mark is not a disguised scalar fugacity: its detach Q velocity has rank-one mark-coordinate action on the extended radical at every site. Nevertheless one mark is insufficient. At both widths three and four the marked affine jet has five moduli and endpoint/radical normalization leaves three, but the Gram restriction is inconsistent with coefficient rank zero and augmented rank one: none of those three moduli changes the offending skew form. The smallest counterexample is width three, where dim(W)=6.

| width | dim V | dim W | rank G0 | dim radical | Q mark gate | affine jet | + endpoint/radical | + Gram | + source | decision |
|---:|---:|---:|---:|---:|:---:|---:|---:|---:|---:|---|
| 3 | 5 | 6 | 3 | 3 | pass | 5 | 3 | empty | empty | `one_mark_insufficient` |
| 4 | 14 | 15 | 3 | 12 | pass | 5 | 3 | empty | empty | `one_mark_insufficient` |

## Exact interpretation

- Width 3: `one_mark_insufficient`. Canonical restricted Gram-skew rank 2; within the scalar endpoint-mark family this requires at least 1 further independent mark covectors (at least 2 total marks) if the canonical line is retained.
  First empty restriction: `endpoint_radical_normalized` -> `gram_self_adjoint` with exact witness `y^T C=0 but y^T b=1`.
- Width 4: `one_mark_insufficient`. Canonical restricted Gram-skew rank 4; within the scalar endpoint-mark family this requires at least 2 further independent mark covectors (at least 3 total marks) if the canonical line is retained.
  First empty restriction: `endpoint_radical_normalized` -> `gram_self_adjoint` with exact witness `y^T C=0 but y^T b=1`.

## Boundary

- Exact rational first-jet algebra only at widths three and four.
- The one-mark state is a terminal response quotient; it does not retain marked-cluster geometry after emission.
- No continuum LCFT, physical transfer matrix, or formal-semigroup K identification is made.
- Failure would motivate the smallest higher marked quotient allowed by the exact skew-rank lower bound, not an unrestricted matrix extension.
