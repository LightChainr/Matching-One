# P334 N9 mark-fiber theorem

## Exact result

All six frozen matching/layer-4/Y=0 rows have the same normalized target-fiber signature.
The old one-carrier/one-mark image, grouped by the normalized first ordered M face, is

| partners in anchor fiber | number of anchors |
|---:|---:|
| 0 | 12 |
| 67 | 8 |
| 70 | 4 |
| 164 | 24 |

Therefore `12*0 + 8*67 + 4*70 + 24*164 = 4752 = 11M`. The coarse source demand is `6912=16M`, so the exact image deficit is `2160=5M`, or `5/16` of demand.

This is a target-image counting identity; it does not invoke maximum flow.

## Two-output-mark closure

After keeping both bases fixed and releasing exactly two output marks, every one of the 48 normalized first-M anchors has all 432 partners. Thus

`M^2/N = 432^2/9 = 20736 = 48M`.

The new image is `15984=37M` beyond the old image.

## Structural correction

The old 11M image is not a union of eleven complete M-sized first-anchor fibers. The 11/16 fraction is an exact global sum over nonuniform fibers, not a literal eleven-of-sixteen local-slot rule.
The finer nonuniform histogram is the theorem: 12 anchors are absent, 8 have 67 partners, 4 have 70, and 24 have 164.

## Crosswalk to the birth-age result

The birth-age review concerns path-time memory of (tau1,ell,tau2) and the direct 0-to-2 collision mass D_N. This certificate concerns the image capacity of a bounded N9 D x F switching reservoir. Neither identity implies the other, and the exact N10 1/57 witness is not rerun here.
The next production statistic should reuse the existing same-batch `(tau1,ell,tau2,DIRECT_RANK2)` archive: test equality of exit hazards across `tau1` within each `(k,ell)` risk set, estimate `D_N=mean(DIRECT_RANK2)` in the same batches, and retain their full covariance. Direct births are a separate no-plateau channel, not missing risk-set rows.

## Boundary

- Exact for the six frozen N9 matching/layer4/Y=0 rows.
- The target-image deficit follows without max flow; source injection still uses the separately certified Hall result.
- No arbitrary-HNF formula or saturation theorem is claimed.
- No slot, source, phase or provenance label contributes target capacity.
