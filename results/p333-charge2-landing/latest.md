# P333/P321/P370 width-four charge-two landing gate

The frozen C4 charge-two landing character fails exactly at the Gram gate: the unique endpoint/radical-normalized affine direction has no Gram-self-adjoint extension. Together with the prior scalar and C4 charge-one results, this exhausts all individual terminal landing irreps at width four.

| dim V | dim mark | dim W | rank G0 | dim radical | affine jet | + endpoint/radical | + Gram | + source/landing | decision |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 14 | 1 | 15 | 3 | 12 | 3 | 1 | empty | empty | `charge2_fails` |

## Exact gates

- New non-scalar rank beyond all scalar columns: 1.
- Translation residual ranks: G0=0, G1=0.
- Canonical restricted Gram-skew rank: 4.
- First empty restriction `endpoint_radical_normalized -> gram_self_adjoint` with `y^T C=0 but y^T b=1`.
- Required semantic change: Retain rooted/landing connectivity after emission, or couple multiple landing irreps in one registry, instead of collapsing each irrep to an independent terminal endpoint mark.

## Boundary

- Exact rational Q=1 first-jet algebra at width four only.
- The single mark is the fixed C4 charge-two landing character, not a fitted endpoint column.
- Success would identify a minimal closed finite-width module only; it would not identify an LCFT field or physical transfer matrix.
- Failure exhausts individual terminal landing irreps, not arbitrary rooted-connectivity or coupled-irrep modules.
