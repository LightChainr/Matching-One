# P333/P398: second jet closes, but the projected gate is now automatic

The nine-mark success is stable through second order, but repeating the same projected gate cannot reveal new structure: charge-one removal has turned it into an involution identity, and J2 is already zero.

## The second-order family is defined

The predecessor protocols define the retained family, not just numerical jets: `Dbar_i(Q)=[[D_i(1)+(Q-1)P_i,0],[(Q-1)E_i,I]]` and `Gbar(Q)=[[G(Q),d_Q G(Q) S_ref],[transpose,0]]`. Commit 5389200 fixes the nine reference columns and the constant row `E_i=C0^T P_i`. Keeping that family makes `Dbar_i2=0` and fixes every Gram coefficient. We do not replace the emission by a Q-dependent row or introduce a new second-order observable.

Taylor coefficients use `t=Q-1`: `G_k=binomial(b,k)` and the cross block is `(k+1)G_(k+1)S_ref`. In particular the second coefficient is half the second derivative.

## Unique second jet

With `X0=Tbar` and `X1=0`, the second affine equation is homogeneous: `X2 A0=B0 X2`. Filtration and fixed mark transport set 207 of 529 entries. The remaining 322 entries have an exact sufficient homogeneous subsystem of rank 312; the same source normalization supplies 10 independent equations, giving rank 322 and uniquely `X2=0`. The projected second Gram residual has rank zero.

| Taylor degree | rank of fixed-radical Gram coefficient | projected Gram-skew rank |
|---:|---:|---:|
| 0 | 0 | 0 |
| 1 | 4 | 0 |
| 2 | 4 | 0 |
| 3 | 4 | 0 |
| 4 | 1 | 0 |

## No next Jantzen layer

The nine-mark Gram has valuations `0^19, 1^4`: `dim J1=4`, `dim J2=0`. The exact leading radical determinant is `14`. This follows from the nondegenerate 4x4 leading form; no second-order coefficient can create a hidden higher layer at Q=1.

## Why more of the same gate cannot discriminate

The induced translation on the surviving radical has three trivial directions and one sign direction, hence `R_B^2=I`. C4 covariance gives `R_B^T H(Q)R_B=H(Q)`; multiplying by the involution immediately gives `H(Q)R_B=R_B^T H(Q)` for the entire polynomial. All five coefficient checks through degree four confirm the identity exactly. Affine covariance likewise keeps `X(Q)=Tbar` constant in the exact linear emission family.

This is the new stop/redirect result: second-jet survival is real, but it adds no independent selection beyond eliminating the radical charge-one block. The next identifying datum must couple to the quotient/physical emission or otherwise impose a genuinely different constraint; repeating higher projected Gram jets cannot do it.

## Boundary

No known unprojected Gram failure was rerun or repaired. The all-order statement is only the same fixed-Q=1-radical projection, not an unprojected Gram module, a moving-radical prescription, or a physical transfer/Jordan realization. If one discards the inherited exact `(Q-1)E_i` family and keeps only its first jet, arbitrary `t^2` emission deformations would be additional unrecorded data; they are not silently admitted here.

```bash
python3 scripts/p398_rooted_second_jet.py
python3 -m unittest discover -s tests -p test_p398_rooted_second_jet.py
```
