# Tiny exact open-boundary block-event oracle

The rectangle consists of adjacent `s x s` cells with open boundaries. A success requires
one unique largest nonempty open cluster in each half and connection of those selected
clusters in the full `2s x s` rectangle.

## Exact enumeration

| graph | s | sites | edges | successes/configurations | P(E) at p=1/2 | coefficients c[k] |
|---|---:|---:|---:|---:|---:|---|
| square | 1 | 2 | 1 | 1/4 | `1/4` | `[0, 0, 1]` |
| matching | 1 | 2 | 1 | 1/4 | `1/4` | `[0, 0, 1]` |
| square | 2 | 8 | 10 | 82/256 | `41/128` | `[0, 0, 2, 8, 19, 24, 20, 8, 1]` |
| matching | 2 | 8 | 16 | 144/256 | `9/16` | `[0, 0, 4, 20, 41, 44, 26, 8, 1]` |

The coefficient vector defines the exact reliability polynomial
`sum_k c[k] p^k (1-p)^(2s^2-k)`.

## Frozen semantics

- an empty half fails because it has no largest open cluster;
- a tie for largest cluster fails;
- largest clusters are chosen inside each half, then connectivity is tested in the union;
- the matching graph adds both diagonals of each unit square;
- no torus wrap or boundary identification is used.

## Boundary

This oracle certifies event semantics only at `s=1,2`. It does not estimate a block-event
probability near `pc`, validate a random sampler, or produce a critical-probability bound.

## Reproduction

```bash
python scripts/pc_block_event_exact.py --format json
python scripts/pc_block_event_exact.py --format markdown
python -m unittest tests.test_pc_block_event_exact
```
