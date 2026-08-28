# Q42 — Amplitude closure on C05 histograms

Issue 35. **Existing C05 data only.** Seed `0xC0100001` (`3222274049`),
**2e6 CRN samples** per same-N Gaussian pair, 40 batches, cross channel.
This is **not** a 1e8 independent confirmation. N=1105 was not started.

## Estimators (same samples, batch covariance retained)

- Reconstruct occupation-space `M(p)=P(K_+<=m)-P(K_->m)`, `m~Bin(N,p)`.
- `p_ref = 0.59274605079210` (Jacobsen 2015 coordinate, not a threshold claim).
- Direct roots: zeros of each orientation's `M_i`.
- Linearized roots: `p_ref - M_i(p_ref)/M_i'(p_ref)`.
- Quadratic roots: solve `M + M' x + M'' x^2 / 2 = 0` at `p_ref`.
- `A_M = N^{13/8} ΔM / Δcos4`, `B = N^{-3/8} mean(M')`, `A_p = -N^2 Δroot / Δcos4`.
- Closure `C_N = -Δroot * mean(M') / ΔM` with all three factors from the **same batch**.
- Train sizes 65, 85, 130; held-out 145, 170. N=145 is retained.

## Per-size closure (batch mean ± SE, B=40)

| N | split | ΔM | mean(M') | Δroot | C_N direct | C_N lin | A_M | B | A_p | A_M/B |
|---|---|---|---|---|---|---|---|---|---|---|
| 65 | train | 0.000218787±0.000403792 | 8.37812±0.00158319 | -2.61881e-05±4.81889e-05 | 0.9998±0.00166639 | 1.00007±0.00167656 | 0.141713±0.261545 | 1.75107±0.000330895 | 0.0811588±0.149341 | 0.0810182 |
| 85 | train | 0.00164377±0.000370388 | 9.25354±0.00166311 | -0.000177676±4.00479e-05 | 1.00559±0.00734902 | 1.00577±0.00733547 | 1.40777±0.317211 | 1.74894±0.000314333 | 0.805103±0.181469 | 0.805097 |
| 130 | train | 9.67594e-05±0.000460184 | 10.8434±0.00229461 | -8.83833e-06±4.24318e-05 | 1.00014±0.0025671 | 1.00024±0.00256839 | 0.193311±0.919377 | 1.74758±0.000369812 | 0.109562±0.525996 | 0.110129 |
| 145 | holdout | 0.000226332±0.00052912 | 11.2897±0.00257252 | -2.0096e-05±4.68477e-05 | 0.997565±0.00432825 | 0.99767±0.0043249 | 0.383871±0.897414 | 1.74651±0.000397965 | 0.220324±0.513618 | 0.219918 |
| 170 | holdout | 0.00091405±0.000456834 | 11.9799±0.0023461 | -7.63791e-05±3.8121e-05 | 0.997291±0.00380765 | 0.997391±0.00380586 | 2.41454±1.20677 | 1.74596±0.000341923 | 1.38439±0.690952 | 1.38265 |

Plug-in C_N from batch-mean (Δroot, M', ΔM) with delta-method SE:

- N=65: 1.00284 ± 0.00642286
- N=85: 1.00022 ± 0.000357264
- N=130: 0.990472 ± 0.0451413
- N=145: 1.00241 ± 0.00792517
- N=170: 1.00105 ± 0.00143296

## Direct vs linearized vs quadratic roots

| N | Δroot direct | Δroot lin | Δroot quad | direct−lin | (direct−lin)/SE_root |
|---|---|---|---|---|---|
| 65 | -2.61881e-05±4.81889e-05 | -2.6202e-05±4.8202e-05 | -2.61951e-05±4.82016e-05 | 1.38981e-08±1.34642e-08 | 0.00028841 |
| 85 | -0.000177676±4.00479e-05 | -0.000177709±4.00553e-05 | -0.00017771±4.00555e-05 | 3.32441e-08±7.87924e-09 | 0.000830109 |
| 130 | -8.83833e-06±4.24318e-05 | -8.8394e-06±4.24367e-05 | -8.83952e-06±4.24365e-05 | 1.06538e-09±6.06877e-09 | 2.5108e-05 |
| 145 | -2.0096e-05±4.68477e-05 | -2.00949e-05±4.68513e-05 | -2.00977e-05±4.68523e-05 | -1.14117e-09±4.58926e-09 | -2.4359e-05 |
| 170 | -7.63791e-05±3.8121e-05 | -7.63831e-05±3.81247e-05 | -7.63851e-05±3.81241e-05 | 4.07541e-09±4.33523e-09 | 0.000106907 |

## Orientation slope difference M'_1 − M'_2

| N | ΔM' | SE | ΔM'/mean(M') | N^{-3/8} ΔM' / Δcos4 |
|---|---|---|---|---|
| 65 | -0.0106403±0.00319799 | 0.00319799 | -0.00127002 | -0.00163124 |
| 85 | -0.00906324±0.00265319 | 0.00265319 | -0.000979435 | -0.00107433 |
| 130 | 0.000597988±0.00396685 | 0.00396685 | 5.51477e-05 | 7.06916e-05 |
| 145 | 0.00156458±0.00509673 | 0.00509673 | 0.000138585 | 0.000126212 |
| 170 | 0.000785207±0.00583195 | 0.00583195 | 6.55438e-05 | 7.17714e-05 |

## Power note on A_M vs C_N

At 2e6 samples, `ΔM` is below one batch SE at N=65, 130 and 145.
`C_N≈1` is therefore a **linearity / same-sample closure** test: batchwise
`Δroot` and `ΔM` share Monte Carlo noise, so their ratio is tight even when
the orientation amplitude itself is not resolved. The slope amplitude `B`
is the well-measured quantity (all five sizes agree at the 10^{-3} level).
`A_M` and `A_p` are **not** stable across sizes at this sample count;
that is a power limitation, not a closure failure.

## Acceptance (honest)

- Independent-seed reproducibility: **not tested here; this is the C05 seed 0xC0100001 at 2e6, not a 1e8 confirmation**
- C_N compatible with ~1 on held-out (within 2 SE): **True**
- Direct and linearized roots agree within 1 root-SE: **True**
- Batchwise covariance propagated: **True**
- No size dropped (N=145 retained): **True**
- Advance the root-shift *amplitude* (stable A_M, A_p): **fail at 2e6** (underpowered);
  advance the *linearization* (C_N, direct≈lin): **pass on holdout**.

Holdout C_N vs 1:

- N=145: C_N=0.997565±0.00432825, z=-0.562688, compatible=True
- N=170: C_N=0.997291±0.00380765, z=-0.711379, compatible=True

## Pell (7,5) secondary calibration

Pell C_N uses linearized roots from a three-point fixed-p scan, so C_N tests slope-asymmetry rather than direct-vs-linear root agreement. Do not mix Pell covariance with C05 batches.

- a7_d5.analysis.json seed=202608280705 n=5000000: C_N(lin)=1.00076, A_p(L^4)=0.566804, A_M=0.988355, A_M/B=0.566375
- a7_d5_rep2.analysis.json seed=202608290705 n=5000000: C_N(lin)=0.999186, A_p(L^4)=0.495132, A_M=0.871121, A_M/B=0.495536
- a7_d5_h0005.analysis.json seed=202608280705 n=5000000: C_N(lin)=1.00073, A_p(L^4)=0.567772, A_M=0.988355, A_M/B=0.567357

## Provenance

- C05 source commit: `f89191a5468ec8417bbeac373334cecb6b5833a7`
- analysis branch: `agent/q42-q43-closure-harmonics`
- wall_time_s: 1.2042782306671143
- machine: cursor Linux 8 cores / 16 GB, no GPU

