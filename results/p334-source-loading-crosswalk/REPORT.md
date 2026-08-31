# Shared-covariance microscopic sources of prevalence and clock loading

| N / readout / source | C | SE | L | SE | D=C+L | SE |
|---|---:|---:|---:|---:|---:|---:|
| 325 / canonical / original_H2_direct | -0.000295243141 | 0.00123567 | 0.00114380046 | 0.000509862 | 0.000848557321 | 0.00140994 |
| 325 / integrated / original_H2_direct | -0.000505491844 | 0.00211571 | 0.00152157974 | 0.000744441 | 0.0010160879 | 0.00241412 |
| 325 / canonical / collective | -6.30851549e-05 | 0.000263817 | -0.000203100942 | 0.000183026 | -0.000266186097 | 0.000254216 |
| 325 / integrated / collective | -0.00024723003 | 0.00103396 | -0.00131605199 | 0.000705854 | -0.00156328202 | 0.00106359 |
| 325 / canonical / unclassified_original_Y | -7.55057255e-07 | 3.14463e-06 | 2.83618947e-05 | 6.65662e-05 | 2.76068375e-05 | 6.65611e-05 |
| 325 / integrated / unclassified_original_Y | -2.78525858e-06 | 1.16815e-05 | -3.81972723e-05 | 0.000126575 | -4.09825309e-05 | 0.000128227 |
| 425 / canonical / original_H2_direct | 0.00137388466 | 0.000836097 | -0.000542004434 | 0.000550899 | 0.000831880227 | 0.00108073 |
| 425 / integrated / original_H2_direct | 0.00234579056 | 0.00142679 | -0.00066690238 | 0.000851322 | 0.00167888818 | 0.00184335 |
| 425 / canonical / collective | 0.000301559534 | 0.000183252 | -4.84464478e-05 | 0.000211249 | 0.000253113086 | 0.000270153 |
| 425 / integrated / collective | 0.00120461082 | 0.000732564 | 0.000516719711 | 0.000744439 | 0.00172133053 | 0.000920604 |
| 425 / canonical / unclassified_original_Y | 1.3990494e-05 | 8.56166e-06 | 3.68959321e-05 | 0.000157614 | 5.0886426e-05 | 0.000157779 |
| 425 / integrated / unclassified_original_Y | 4.6258431e-05 | 2.75872e-05 | 4.47043742e-05 | 0.000199933 | 9.09628052e-05 | 0.000197766 |

Same 20 original batches, not independent source tests. Full joint covariance is saved, but never inverted and no high-dimensional omnibus is computed. Direct means a final site in the original checkpoint H2 set; later-created triggers belong to collective. All whole-pair fallbacks remain unclassified. Point signs or cancellations are not asserted as statistically established. This is gated R1 only, not full A_top.
