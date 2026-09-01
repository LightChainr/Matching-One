# Stabilizer-free generic-modulus homology-vector preflight

This is a deterministic theory/design calculation. It generated zero random samples.

## Selected design

- Base period matrix: `[[10, 3], [0, 10]]`, tau=`3/10+i*1`, area 100.
- The base and all four signed-stencil cells have exact orientation-preserving stabilizer `{+I,-I}`.
- Directions: an exact shear pair and an independent oblique-aspect pair, all at equal area/cost.
- Full output schema has 19 theory channels here and requires sparse retention of every finite-sample primitive winding line.

## Frozen vectors and information score

`mu_KdV`, `mu_Q4_Jordan`, and continuum-subtracted zero `mu_embedding` are stored in `latest.json`.
After profiling the two continuum-tangent leakage columns, KdV versus Q4/Jordan has unit-Fisher shape D^2 = `0.66366006` (absolute cosine `0.66816997`).
The timing-based four-cell lower-bound cost is `0.0112135142` CPU-hours per one million samples/cell-equivalent, giving a conditional shape-value proxy `59.1839497` per CPU-hour.
The strict three-model maximin V remains `0`: the embedding prediction is zero after continuum subtraction and no nonzero operator-amplitude floor is frozen. This does not erase the nondegenerate KdV/Q4 shape result.

## Frozen finite-pipeline contract

Save P0, the sparse rank-1 primitive-winding vector, P2 and all thermal derivatives. Form q=P2-P0, E=P0+P2, refind the pooled q=0 root, use D=d_p q as the common physical normalizer, continuum-subtract each directional homology vector, and form the channelwise U vector at that root.

## Decision

The exact geometry and theory templates are ready for one bounded N100 full-vector covariance pilot. Large production is not started by this result.

Execution provenance: local fallback after TV2N0X rejected the existing SSH key; no key reset was attempted.
