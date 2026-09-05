# New N100 production: shape splitting and a hidden thermal redistribution

**The three-modulus experiment has produced a new response, not just another
model ranking.** At this finite area the four observables cannot all be an
offset plus an amplitude times one shared scalar shape function. Yet their
integrated clock C follows the E4 shape closely. Resolving the existing stream
along p shows how these statements coexist: the clock-calibrated odd response
has exactly zero area but a strongly resolved thermal dipole and higher moments.

This combines a pre-acquisition fixed-p test with an explicitly post-reveal
functional analysis. They share data and are not independent discoveries.

## 1. The actual experiment

The three shapes are `tau=2i,4i,1/2+i`, each at N100, with the same rational
rotation O=(1/5)[[4,-3],[3,4]] and Smith pair (1,100) -> (5,20). Period matrices,
signed orientation normalization and counter allocation were frozen in
`experiments/etop_n100_three_modulus_20260831.json` at **4c1ec50**. Runner
**295ac79** precedes all new target generation.

Each shape pair used 2,000,000 new permutations and 200 aligned batches. The
same seed and counter block align the permutations across all three shapes;
these are six million pair evaluations, **not six million independent common
counter draws**. Every shape pair compares two orientations. The entire
experiment is one correlated dependency group, with its 12x12 covariance saved.
The three local single-thread processes completed in approximately 20--29 s.
No Huawei connection, paid GPU, new engine or repeated old test campaign was
required. Raw threshold histograms, rank moments and receipts are included.

Write Y=(A_top,E_top,C,W), D=Y(4i)-Y(2i), U=Y(1/2+i)-Y(2i).
Any common affine shape profile implies U=rD. The E4 candidate fixes
r=-0.287083852577789; height-only E4 gives -0.212358001359808 and y^2 gives -1/4.
These exact three-point relations eliminate a separate offset and amplitude
for every coordinate. They do not require old N50 calibration or an area law.

| N100 hypothesis | chi-square / df | Gaussian reference p |
|---|---:|---:|
| affine actual-modulus E4, frozen primary | 58.26869 / 4 | 6.70e-12 |
| affine height-only E4 | 24.42330 / 4 | 6.57e-5 |
| affine y^2 | 23.89369 / 4 | 8.39e-5 |
| **any free common scalar secant**, retrospective relaxation | **19.21808 / 3** | **2.46e-4** |

The useful conclusion is not that y^2 wins a ranking. **Even the flexible
single-profile class misses the joint response.** A second finite shape
direction is useful; three shapes then permit a saturated representation, so
this is not a count of continuum fields. A nonlinear vector curve is also
outside the rejected affine class and remains possible.

## 2. C is a useful shape coordinate despite that joint failure

For the fixed E4 residual U-r_E4 D, the separate marginal standardized
residuals are A=-2.468, E=0.381, **C=-0.467**, W=2.627. By comparison, the
height-only E4 residual for C is 3.515 standard errors. These marginal views
are correlated and do not sum to the joint statistic.

Instead of assigning another field immediately, define the empirical clock
coordinate r_C=U_C/D_C. The new stream gives

```text
r_C = -0.277981748 +/- 0.01937.
R_j = U_j-r_C D_j.
```

This quotient retains a fixed physical interpretation: how much shear response
remains after its integrated-clock change is explained by the 2i -> 4i change.
It is estimated from the same stream. Ratio uncertainty and all cross-shape
covariance are retained; it is not a new frozen independent test.

## 3. The exact integral explains the apparent tension

Let K1,K2 be the first/second topological threshold ranks and
F_k(p)=Pr[Binomial(N,p)>=k]. Then

```text
integral_0^1 F_k(p) dp = 1-k/(N+1),
A_top(p) = F_K1(p)+F_K2(p)-1,
E_top(p) = 1+F_K2(p)-F_K1(p),
C = (E[K1]+E[K2])/(2(N+1)),
W = (E[K2]-E[K1])/(N+1).
```

After taking the normalized orientation contrast P4, constants cancel:

```text
integral P4[A_top(p)] dp = -2C,
integral P4[E_top(p)] dp = -W.
```

Consequently the clock-quotient functions satisfy

```text
integral R_A(p) dp = 0 exactly,
integral R_E(p) dp = -R_W.
```

**A matching clock therefore does not imply a matching thermal profile.**
Positive and negative pieces of R_A may cancel in the integral while remaining
large throughout p. In birth coordinates R_F1=(R_A-R_E)/2 and
R_F2=(R_A+R_E)/2: their signed areas are +R_W/2 and -R_W/2. Clock cancellation
can coexist with opposite first-/second-birth area shifts.

## 4. The new stream resolves the redistribution

The full threshold histograms permit analytic thermal moments without new
sampling or quadrature. With z=N^(3/8)(p-p_ref) as a convenient finite-N
coordinate (not a fitted scaling exponent):

| clock-quotient moment | estimate | SE | estimate/SE |
|---|---:|---:|---:|
| integral R_A | 0 by construction | -- | -- |
| integral z R_A | -0.0002936353 | 0.000046967 | -6.25 |
| integral z^2 R_A | 0.0008485571 | 0.00010592 | 8.01 |
| integral R_E | -0.0001574519 | 0.000065885 | -2.39 |
| integral z R_E | 0.0004207693 | 0.000071136 | 5.92 |

The empirical odd curve has positive and negative lobe areas approximately
+0.000347135 and -0.000347135. Core crossings are near p=0.42428 and 0.77193.
The largest mean lobe is near p=0.32125, not near p_ref=0.59274605079.
**This is a finite-geometry thermal redistribution, not yet a critical-window
operator identification.** Roots, peak locations and lobe areas are descriptive
post-reveal summaries; the moments and pointwise errors use the correlated
influence of the fitted clock ratio. No p-grid significance votes are pooled.

## 5. Scientific card and next discovery

- **Changed mechanism space:** a common scalar affine modulus response is too
  small for this N100 A/E/C/W vector. An integral-clock component plus an odd
  zero-area thermal redistribution and an even lifetime-area deformation is a
  more informative working description.
- **Not established:** exact E4 field identity, number of continuum operators,
  scaling limit, or critical-window origin of the largest residual.
- **Observer/sector/source/geometry:** normalized quotient-orientation A_top
  odd and E_top even; clock C and lifetime W; new threshold-rank histograms;
  the three frozen N100 period pairs.
- **Dependency:** all fixed-p, all-p, model and quotient readouts share seed
  20260831125401, offset 267100000000. Old N50 transfer is a separately modeled
  source; source and target uncertainty are both included there.
- **Next discriminant:** keep the three homothetic shapes and distinguish
  critical-window moments from the full-p dipole at a larger common scale.
  Does the zero-area odd redistribution move into the critical window, or does
  it remain a thin-geometry/off-critical effect? A blind sign/shape prediction
  from the present clock-quotient curve is more informative than another vote
  among scalar harmonics. N400 by uniform period doubling preserves all three
  Smith lineages; naive Gaussian norm-two scaling to N200 does not preserve a
  common Smith pair across these three shapes.

Artifacts: `results/etop-n100-three-modulus/` contains the new production and
joint scores. `results/etop-thermal-redistribution/` contains exact integrated
moments and full compact-grid covariance. Scripts are
`run_etop_n100_three_modulus.py` and `etop_thermal_redistribution.py`.
