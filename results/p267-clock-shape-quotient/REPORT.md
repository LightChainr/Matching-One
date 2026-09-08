# N100 clock-anchored shear deformation

Zero new samples; post-reveal reuse of the existing 200 paired batches.

Clock shape r_C = -0.277981748 +/- 0.019377.

| residual | estimate | SE | marginal z |
|---|---:|---:|---:|
| A_top | -0.00089960894 | 0.00035015 | -2.5692 |
| E_top | 9.52915132e-05 | 0.0002448 | 0.38926 |
| W | 0.000157451869 | 6.5876e-05 | 2.3901 |

Joint clock-closure diagnostic: chi2=17.8234/3, nominal p=0.000478321.
This localizes the coordinator's common-secant failure; it is not an independent test.

The exact integral identities give integral R_A(p) dp=0 and integral R_E(p) dp=-R_W.
The odd zero is imposed by the clock gauge. The even area estimate is nonzero in sign but only 2.39 SE; it does not independently reject zero at alpha .01.

Under fixed E4, C contributes 0.218432 chi-square and conditional A/E/W contribute 58.0503/3. Conditioning on a true E4 clock is stronger than using an empirical clock.

## Which deformation direction?

Marginal A is negative, W positive and E unresolved. The complete covariance matters: E/W residual correlation is about -.763.
The descriptive single-axis nuisance fits are:

- A_only: chi2=17.4445/2, nominal p=0.000162924.
- E_only: chi2=9.48118/2, nominal p=0.00873351.
- W_only: chi2=6.71686/2, nominal p=0.0347899.

A W-only deformation remains compatible at .01; a unique A/E/W loading is not identified. These source-selected diagnostics are not corrected independent model elections.
The source-selected covariance-matched readout is Psi=-0.022459784 R_A + 0.22250609 R_E + R_W. Freeze it for future data; its in-sample SNR is not a new significance result.

## Minimal representation and next falsifiers

After subtracting Y2, use clock profile [0,1,r_C] and shear profile [0,0,1]. Their loadings are D and R. Three shapes make this saturated: it is a transparent coordinate choice, not two-field evidence.
A next independent same-semantics block can test the source R direction via two cross-products, retaining source uncertainty. A stronger optional law keeps R/D_C fixed; its three source values and full covariance are in score.json. No area exponent, same-lineage homogeneity or new production is assumed.

## Scientific card

- Changed space: a readable C clock does not close the fixed-p A/E/W morphology; the remainder has an exact zero-area odd interpretation and an even lifetime-area coordinate.
- Not proved: independent fields, physical Jordan identity, universal two-profile closure, or a separately significant nonzero even integral.
- Observer/source/geometry: signed normalized A/E/C/W rank/clock contrasts at p=.59274605079, N100, three fixed moduli and the same cyclic/noncyclic O map.
- Dependency: one common-random 2M/shape, 200-batch block; all new readouts are correlated post-reveal reuse.
- Next upgrade: an independent geometry/scale block that preserves the frozen deformation direction or rejects it; whole-p curves separately determine how the zero-area redistribution is carried.
