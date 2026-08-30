# P398 width-four rooted-connectivity grade-one completion

The AP/OP/DP rooted connectivity registry gives a minimal exact completion of the width-four Q-adic grade one. It adds precisely 2 trivial, 2 charge-one and 3 charge-two directions to the previously tested rank-six responses; the resulting 13-by-13 B-coordinate matrix has full rank and determinant 3072. Thus the seven-dimensional gap is resolved by multiplicity-bearing source-to-landing connectivity inside the existing C4 sectors, without another terminal character.

| sector | old rank | rooted raw rank | incremental rank | target | complete |
|---|---:|---:|---:|---:|---|
| trivial | 3 | 2 | 2 | 5 | `true` |
| charge1_rational | 2 | 2 | 2 | 4 | `true` |
| charge2 | 1 | 3 | 3 | 4 | `true` |

## Exact completion gate

- New rooted coordinates: 7; dimension lower bound: 7.
- Old tested rank: 6; rooted raw rank: 7; combined rank: 13/13.
- Exact determinant of the combined B-coordinate matrix: 3072.
- All raw responses and H-dual grade-one vectors have zero exact C4 translation residual.

## Rooted registry

- `AP`: source-to-adjacent landing pair, with the remaining sites singleton; supplies trivial, charge-one and charge-two projections.
- `OP`: source-to-opposite landing pair; supplies an additional trivial and charge-two projection.
- `DP`: source-adjacent pair plus complementary landing pair; supplies the third missing charge-two copy.

## Boundary and next gate

- Promote this seven-coordinate rooted registry to a declared extended module and test the affine/endpoint/radical/Gram/source intersection. The present result proves only that it is the smallest response space with enough associated-grade information.
- Width four and the already-certified grade-one layer only.
- This is an exact span/completion certificate, not yet an affine endpoint/Gram/source closure test for a new extended module.
- No terminal-character family, random production, continuum field or physical Jordan interpretation is added.
