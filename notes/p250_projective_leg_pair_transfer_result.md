# P250 fresh projective-leg pair-transfer result

## Result

The pair-only turn succeeded.  A fresh 40k counter stream on the XP Huawei
node resolves every hand-charge real transfer row through `d=3`; the weakest
one is 7.995 standard errors.  This establishes the projective-leg insertion
as a practical mesoscopic charged propagator even though its prior cubic row
was not a useful OPE observable.

The frozen amplitude-free scores over `d=1,2,3` are:

| shape | joint score | reading |
|---|---:|---|
| one complex transfer eigenvalue | 19.124 / 8, p=0.0142 | appreciable tension |
| within-channel power shape | 9.297 / 8, p=0.318 | compatible |
| nonzero deck phase | 5.846 / 8, p=0.664 | not resolved |
| constant deck-phase step | 3.690 / 4, p=0.450 | compatible with zero/constant step |

The clean scientific statement is therefore not an OPE claim: **the charged
pair row propagates, its first three separations look more power-like than a
single exponential transfer eigenvalue, and the complex data do not resolve a
deck-character phase.**  The power score is channelwise; it does not yet claim
one universal exponent shared by both hands and charges.

The individual effective powers are also useful coordinates for the next
design.  From `d=1->2` and `d=2->3`, respectively, they are `1.211/0.956`
(plus r1), `1.259/1.132` (plus r2), `1.163/1.547` (minus r1), and
`1.123/1.361` (minus r2).  A larger parent, rather than more replicas at N325,
is the natural next discriminator between a common scaling dimension and a
short-distance state mixture.

## Periodic-distance boundary

The parent geometry is the norm-65 Gaussian torus with periods
`((8,-1),(1,8))`.  Thus `d=4,5,6` already probe the antipodal/periodic-image
regime; the rise in real-row resolution at `d=5,6` is not a long-distance
nonmonotonic tail.  Those rows and the full 96-coordinate covariance are kept
as finite-torus transfer data, while the frozen shape decision intentionally
uses only `d=1,2,3`.

## Reproduction capsule

- freeze commit: `f35a0c0cd8de1656efa064b03a4aeb21d6dae41e`
- runner/scorer commit: `47d086f35f71ec70e66d4bce0610d72de0b447bf`
- seed/counters: `25025033720260931`, `[0,40000)`
- host/path: `Huawei-CodeBuddy-XPk2PZ`,
  `/workspace/p250-projective-leg-pair-transfer-40k`
- raw batch SHA-256:
  `6a5ab14f7a26accdd1329d3ee7704e84400685ab6634738bc9065b2487e1da4b`
- response SHA-256:
  `4af91032bdd4c489a29511cd4305fb9c9793519ecc7aa5f970fb30dfd2d66c6b`

The raw response contains complex `T` and axis-difference `A` rows for
`d=1..6`; the CSV batches supply their complete joint covariance.  Cubic
fields are absent by construction.
