# Flat-twist projective representation

All exact permutation, projector and irreducibility gates pass.

## Exact decompositions

- F2 under D4 image: `2 A1 + B1`; under S,T: `1 + 2D standard irreducible of PSL(2,F2)=S3`.
- F3 under D4 image: `2 A1 + B1 + B2`; under S,T: `1 + 3D standard irreducible of PSL(2,F3)=A4`.

The balanced axes-minus-diagonals H4 alias is a second A1 copy under D4, not a symmetry-distinct scalar irrep. At F3 the axis-odd and diagonal-odd lines are B1 and B2. The full modular image mixes all three non-scalar D4 coordinates into the irreducible 3D standard representation of A4.

Restricting that triplet to the order-3 shear gives one real neutral line and one real 2D charged rotation block. The diagonal-odd coordinate is one axis of this charged doublet, not a standalone modular field.

## Placement of the N65 contrast

D4 basis `[H4,axis-odd,diagonal-odd]`: `[0.0011029842067728016, 0.005129985113886593, 0.0020775612521830144]`.

T-cyclic basis `[neutral,charged-v,charged-w]`: `[0.004825423534390718, -0.002061215452793925, 0.0020775612521830144]`.

The covariance-aware charged-doublet diagnostic is `5.905028 / 2 df`; the neutral coordinate is `z=1.602`.

The approximately two-sigma diagonal-odd marginal is therefore inseparable from its charged partner under a realizable modular shear. It is not promoted to a discovery.

## Minimal no-fit prediction

For a covariantly transported T-shear source, the exact mixing is

```text
H' = (axis_odd - diagonal_odd)/sqrt(2)
A' = H/sqrt(2) + (axis_odd + diagonal_odd)/2
D' = H/sqrt(2) - (axis_odd + diagonal_odd)/2
```

so the observed N65 vector predicts `[0.0021583896116662287, 0.004383700795185517, -0.0028238455708840906]` with no fitted amplitude. A future charged/projective-source run should freeze this whole vector and its transport, not score diagonal-odd alone.
