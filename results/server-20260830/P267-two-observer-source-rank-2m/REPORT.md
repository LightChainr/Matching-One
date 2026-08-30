# P267 Target 1: two-observer x two-source rank production

## Integrity

- Frozen source and scorer commit: `9f2033a2d205a01d3554bc2e9c93bed5760f5b29`.
- N325: 2,000,000 replicas, 100 batches, seed `202608303252`, counter
  interval `[19000000000,19002000000)` on Zy.
- N425: 2,000,000 replicas, 100 batches, seed `202608304252`, counter
  interval `[21000000000,21002000000)` on XP.
- The ARM64 binary hash is identical on both hosts.  Each host's generated
  `SHA256SUMS.txt` passes in place; all locally retained files and decoded
  gzip streams reproduce those remote hashes.  The roughly 0.5 GB per-size
  marked-birth streams remain on the execution hosts as frozen by the
  preregistration, with their hashes retained in `SHA256SUMS.txt`.
- Both 200-row complement audits have zero endpoint, site, line, local-mark,
  index and separated-mark failures.  The regenerated exact mapping gate has
  the preregistered SHA256 `5afeef4026...` and rank-two response determinant
  `-2`.

## Frozen reveal

The primary statistic is the complex determinant of the P4-projected
two-observer x two-source coupling matrix, with delete-one recomputation of
the intrinsic root and JD/JS Gram orthogonalization.

| size | determinant | chi-square/df | p | normalized wedge |
|---:|---:|---:|---:|---:|
| 325 | `0.08083 + 0.16172 i` | 1.263/2 | 0.532 | 0.139 |
| 425 | `0.11570 - 0.33244 i` | 3.016/2 | 0.221 | 0.493 |

Because the two sizes use disjoint counter domains, the frozen joint score is
their block-diagonal sum: chi-square `4.279/4`, p `0.370`.  This is well below
the preregistered rejection threshold `18.4668`; the rank-one/common-projective
lane is not rejected.  Neither determinant is individually resolved, so the
conditional N425/N325 transfer was not scored.

This is an informative null for the mechanism selector.  Adding a second,
exactly non-aliased spatial-H4 observer did not expose a two-dimensional
JD-perp/JS source plane at these sizes and precision.  It does not prove that
only one source exists, identify either column with a continuum field, or
exclude a second direction that is weak or projectively aligned in this
observer basis.

## Five-line scientific card

1. **Question:** do O_far and the non-aliased O_sep4 observer see a genuinely
   two-dimensional JD-perp/JS source plane?
2. **Exact gate:** both direction orbits are disjoint, source-separated and
   rank two before data; all production audits pass.
3. **Result:** joint determinant chi-square `4.279/4`, p `0.370`, versus the
   frozen rejection threshold `18.4668`.
4. **Interpretation:** the simplest common projective lane survives; no
   individual matrix entry is promoted to a field identity.
5. **Next selector:** a future experiment must rotate the source or observer
   basis by an independently motivated mechanism, rather than repeat these
   two rows or fit a post-hoc mixture.

