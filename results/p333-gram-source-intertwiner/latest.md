# P333/P370 Gram-source affine intertwiner intersection

## Exact dimensions

| width | dim V | affine Hom | endpoint-normalized | + Gram | + source | X=T in final | decision |
|---:|---:|---:|---:|---:|---:|:---:|---|
| 2 | 2 | 2 | 1 | 1 | 0 | yes | `canonical_unique` |
| 3 | 5 | 3 | 2 | 1 | — | no | `empty_intersection` |
| 4 | 14 | 9 | 8 | 5 | — | no | `empty_intersection` |

The reported normalized dimensions are affine tangent dimensions. Exact rational particular solutions and primitive-integer tangent bases are stored in `latest.json`.

## Interpretation

Width two selects X=T only because translation is identity in the two-state quotient. At the first nondegenerate widths three and four, T is not first-jet-Gram self-adjoint, and the entire Gram-compatible affine sigma-Hom slice is incompatible with fixing the all-singleton source. The frozen physical intersection is empty there.

- Width 2: canonical `T` has Gram residual rank 0; final decision `canonical_unique`.
- Width 3: canonical `T` has Gram residual rank 2; final decision `empty_intersection`.
  Exact inconsistency witness: `y^T C=0 but y^T b=1` with nonzero state coefficients `[{'state': [0, 0, 0], 'coefficient': -1}, {'state': [0, 0, 1], 'coefficient': -2}]`.
- Width 4: canonical `T` has Gram residual rank 4; final decision `empty_intersection`.
  Exact inconsistency witness: `y^T C=0 but y^T b=1` with nonzero state coefficients `[{'state': [0, 0, 0, 0], 'coefficient': -1}, {'state': [0, 0, 0, 1], 'coefficient': -2}, {'state': [0, 0, 1, 1], 'coefficient': -2}, {'state': [0, 0, 1, 2], 'coefficient': -1}]`.

This intersection is distinct from the signed-history nilpotent `K=(D-J)^2`: K lives in the three-mark endpoint radical, whereas this calculation selects off-diagonal maps between crossed/trivial affine closure modules.
It is also distinct from PR #385: an exact selected line does not by itself make Jordan behavior statistically separated from near-colliding ordinary models.

## Boundary

- This is exact rational algebra only in the width-2,3,4 noncrossing join representation.
- The first-jet Gram is restricted from the full connectivity pairing; no generic-Q detach or loop-weight generator is added.
- Source normalization fixes the declared all-singleton source and must not be reinterpreted as a measured thermal insertion.
- Uniqueness here would not overcome PR385 finite-noise nonseparation without an independently justified physical spectral or symmetry condition.
