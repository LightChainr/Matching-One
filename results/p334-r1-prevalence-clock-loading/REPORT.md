# R1 prevalence versus conditional-clock loading

| N / readout | C: prevalence | SE | L: conditional clock | SE | C+L: H4 difference | SE |
|---|---:|---:|---:|---:|---:|---:|
| 325 / p_ref | -0.000359083353 | 0.00150262 | 0.0009690614143 | 0.000396122 | 0.0006099780613 | 0.00156001 |
| 325 / p_integral | -0.0007555071328 | 0.00316134 | 0.0001673304792 | 8.07387e-05 | -0.0005881766536 | 0.00315435 |
| 425 / p_ref | 0.001689434689 | 0.0010278 | -0.00055355495 | 0.000431888 | 0.00113587974 | 0.00118417 |
| 425 / p_integral | 0.003596659808 | 0.00218668 | -0.0001054782945 | 7.51176e-05 | 0.003491181513 | 0.00220881 |

## Variance explained by the four R1 pair states

| N / readout | Between-state variance | Within-state variance | Between fraction | Batch SE |
|---|---:|---:|---:|---:|
| 325 / p_ref | 0.02651655489 | 0.004961207143 | 84.23901% | 0.2844 percentage points |
| 325 / p_integral | 0.1180752768 | 0.000178619844 | 99.84895% | 0.0037498 percentage points |
| 425 / p_ref | 0.01959475538 | 0.003992823162 | 83.07235% | 0.30556 percentage points |
| 425 / p_integral | 0.0880216622 | 0.0001200477664 | 99.8638% | 0.0037413 percentage points |

These variance terms describe individual hybrid contrasts, not their mean-estimator covariance. The 20-batch joint uncertainty for both decompositions is saved separately.

R1-layer weighted hybrid Y only: conditional suffix means on solved pairs and original suffix readouts on whole-pair fallbacks. Not full F2/A_top, not a causal decomposition, not an independent random block, and no attribution percentages when C and L cancel. Sizes are reported separately; no cross-size independence is assumed here.
