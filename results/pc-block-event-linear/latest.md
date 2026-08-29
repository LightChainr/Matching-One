# Linear open-boundary block-event evaluator

The deterministic evaluator uses two disjoint-set forests: one for connectivity in the
full rectangle and one with cross-half edges suppressed for largest-cluster selection.
It has `O(s^2)` time and memory complexity.

## Differential validation

- configuration/graph pairs compared: `520`;
- exact-oracle mismatches: `0`;
- covered sizes: `s=1,2`; graphs: square NN and matching NN+diagonals.

## Deterministic large controls

| graph | s | selected sites per half | full edge checks | half edge checks | result |
|---|---:|---:|---:|---:|---|
| square | 64 | 4096 | 16192 | 16128 | success |
| matching | 64 | 4096 | 32194 | 32004 | success |

The frozen square-graph nonmonotonicity witness is also preserved: mask `6` succeeds,
while opening one additional site to form mask `22` fails with `left_tie`.

## Boundary

No stochastic sampler is included. This result does not establish independent trials,
an event probability near `p_c`, a runtime guarantee, or a new certified threshold bound.
