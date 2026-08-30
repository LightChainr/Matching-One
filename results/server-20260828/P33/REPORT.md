# P33 threshold-rank thermal tomography

The production C++17 kernel implements bidirectional Newman--Ziff threshold
ranks for the five frozen same-N orientation pairs (`N=65,85,130,145,170`).
It stores integer `K_minus` and `K_plus` marginal histograms for every
orientation and batch, plus the joint first and second rank moments.  The two
orientations share the same counter-keyed Fisher--Yates permutation.

Exact and reproducibility checks cover all 120 permutations of the primitive
Gaussian `N=5` torus, bin-by-bin equality with the Python reference on shared
counters, and byte-identical output with one and two worker threads.

## Local one-million pilot

An Apple ARM single-core pilot at `N=65` used 1,000,000 paired permutations,
100 batches, seed `2026093301`, and counters beginning at 2,000,000,000.  It
finished the connectivity calculation in about seven seconds.  At the frozen
central probability:

```text
p_ref                         0.592746050790
Delta M (first - second)      +1.9810e-3 +/- 5.5575e-4
Delta M' (first - second)     -8.4709e-3 +/- 3.6780e-3
root first                    0.5926749318571
root second                   0.5929113566735
root gap (first - second)     -2.36425e-4 +/- 6.63596e-5
```

At symmetric offsets `p_ref +/- 0.001`, the two orientation differences are
`1.98962e-3` and `1.97268e-3`.  Their central-reflection even component is
`1.98115e-3`; the odd component is `-8.47e-6`, consistent with the independently
reconstructed slope difference times the offset.  This demonstrates the
thermal-parity pipeline and leaves a central orientation amplitude after the
common thermal coordinate is removed.  The pilot is not the frozen production
confirmation: its Monte Carlo uncertainty is still about 28% of the signal.

Sparse joint histograms were deliberately omitted to keep output bounded; the
retained marginal histograms reconstruct `M`, `M'`, and roots, while joint
moments preserve threshold-rank covariance diagnostics.

## Huawei five-size 10-million pilot

The 16-core ARM server then ran all five frozen pairs with 10,000,000 paired
permutations per size and 100 batches.  The connectivity phase took 65.2
seconds total: 767,154 paired permutations/s, 182.6 million nominal
orientation-sites/s (counting two orientations), or 365.2 million nominal
bidirectional site-ranks/s.  The retained histograms, joint moments, and
metadata occupy 4.19 MB, about 41.9 kB per cross-size batch or 8.4 kB per
geometry/batch.  Relative to the Python L=8 reference throughput, the measured
aggregate paired-permutation rate is over 600x larger, despite including the
larger N=85--170 geometries.

| N | central Delta M | SE | root gap | SE | N^2 root gap |
|---:|---:|---:|---:|---:|---:|
| 65 | +7.16535e-4 | 2.03e-4 | -8.55350e-5 | 2.42e-5 | -0.361 |
| 85 | +7.30300e-4 | 2.02e-4 | -7.88704e-5 | 2.19e-5 | -0.570 |
| 130 | +4.49692e-4 | 2.36e-4 | -4.14783e-5 | 2.18e-5 | -0.701 |
| 145 | +6.40070e-4 | 2.09e-4 | -5.66827e-5 | 1.85e-5 | -1.192 |
| 170 | +5.82797e-4 | 2.02e-4 | -4.86397e-5 | 1.69e-5 | -1.406 |

All central amplitudes have the sign predicted by `Delta cos(4 theta)`, and
all root gaps have the corresponding negative sign.  Across every size the
central-reflection odd component from `p_ref +/- 0.001` agrees numerically with
`0.001 * Delta M'`; the much larger even component therefore survives removal
of the common thermal coordinate.

The frozen root-gap acceptance is not yet met.  A constant `N^2 Delta p`
trained on N=65,85,130 gives held-out standardized residuals -1.89 and -1.95
at N=145,170 (joint chi-square 7.40 for two held-out rows).  The ten-million
pilot resolves the signs and thermal parity, but shows substantial radial
drift; a production root-gap scaling claim needs more statistics and/or the
declared correction models.
