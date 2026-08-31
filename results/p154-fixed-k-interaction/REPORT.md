# Norm-4 topology under a fixed-occupation local interaction

## The new interaction response remains unresolved at 20k

Neither the full local interaction response nor its fixed-K part is resolved on these inherited blocks. The fixed-K z-scores are 0.152 and 0.891; this supplies no positive evidence for a second physical direction, and does not show that the response is explained by the thermal-count term.

This reanalysis recovers the old Phase-E configurations and adds the missing `E_top × edge interaction` moment. It does not rerun the former connectivity-B scorer or claim a new energy field.

For each of the two archived sizes, the reported direction difference is divided only by the exact cosine-four difference. Values below are the response per mean NN edge; uncertainty is aligned delete-one-batch standard error.

### N=65

- `J_Q` = 2.27488706e-05 ± 8.42e-05 (z=0.270).
- `J_R` = 1.22528773e-05 ± 8.08e-05 (z=0.152).
- `J_H` = 1.04959933e-05 ± 1.2e-05 (z=0.872).
- `J_R0` = 3.52848747e-05 ± 5.57e-05 (z=0.634).
- `J_R2` = -2.30319975e-05 ± 6.15e-05 (z=-0.375).

### N=130

- `J_Q` = 3.55115467e-05 ± 4.9e-05 (z=0.725).
- `J_R` = 4.27567326e-05 ± 4.8e-05 (z=0.891).
- `J_H` = -7.24518595e-06 ± 6.72e-06 (z=-1.079).
- `J_R0` = 4.06555753e-05 ± 4.23e-05 (z=0.960).
- `J_R2` = 2.1011573e-06 ± 4.34e-05 (z=0.048).

The two-size fixed-K mixed-response joint zero statistic is `0.8173/2`, nominal p=`0.664547`. This uses an estimated covariance and is not an exact finite-sample model certificate.

## What distinguishes the three sources

Write K for occupied sites, T for occupied NN edges and p=.59274605079. There are 2N simple NN edges in each archived torus.

```text
Q = (T − 4pK + 2Np²)/(2N)
R = (T − 2K(K−1)/(N−1))/(2N)
H = Q−R
J_X = Cov(E_top,X), E_top=I0+I2
J_Q = J_R + J_H
```

Q is a local pair interaction with positive finite-volume measure proportional to Bernoulli(p) × exp(λQ). Its λ derivative of E_top is J_Q. R removes the entire occupation-count conditional mean: E[R|K]=0. Conditional on K, Q and R generate exactly the same interaction. R includes a global K counterterm and is not itself a strictly local canonical field.

For any p-independent observable O, `Cov(O,H)=p²(1−p)²/[N(N−1)] × d²E_p[O]/dp²`. Thus H is an explicit second-thermal-score projection, whereas J_R measures how the sector distinguishes local edge arrangements at the same occupied count. This finite-product-measure separation does not prove vanishing RG thermal overlap.

J_R0 and J_R2 retain the rank-zero and rank-two contributions separately. Their sum is J_R. Neither a cancellation nor an unresolved sum means the underlying two mixed responses vanish.

## What this adds to the existing pilot

The old B asks whether nearby endpoints are connected through the complete torus; it is a legitimate global-connectivity readout, not a finite local occupation-cylinder function. Its unresolved mixed H4 response did not prove algebraic dependence. P40 stores the q × fixed-K motif second moment, where q=I2−I0 is configurationwise. E_top=q² requires a third moment absent from its Gram archive; it is not the square of the expected matching function. The present E_top × R measurement is therefore a different scientific question, not a renamed variance-reduction score.

Matching parity is stated for the fixed Euclidean stencil under the full graph-pair/complement transformation `(G_black,G_white,η,p)→(G_white,G_black,1−η,1−p)`. It does not identify the NN edge set with that of whichever graph becomes black.

## Sampling and reproducibility

The original 0578105 backend blob is unchanged. The original N65/N130 seeds and 20000 counter ranges are specified in the manifest. All 400 batch/orientation rows reproduce the original sample counts, K1/K2 sums and I0/I1/I2 sums exactly. Only local edge sufficient statistics were added. Parent and child results share their random streams; they are not independent evidence.

Pooled unbiased sample covariances are recomputed after removing each aligned batch from both orientations. The JSON retains all leave-one-out vectors and the full covariance, including the exact J_Q/J_R/J_H linear dependence. No inverse of the redundant full matrix is used. Sizes have different original seeds and remain separate covariance blocks.

## Interpretation and next physical question

Use the measured split to decide whether an observed local interaction response is explained by the occupation-count second derivative or requires fixed-K geometry. An unresolved result is a completed measurement at these inherited 20k blocks, not a prohibition on other singlets or scales. A resolved fixed-K row would still need transport across geometry/scale and overlap with the original norm-4 residual before a field could be named. This analysis does not fit an exponent, select a continuum operator or start additional production.

Reproduce by compiling `src/p154_fixed_k_interaction_replay.cpp`, replaying N65 and N130 to `raw/n65.csv` and `raw/n130.csv`, then running this script. The runner refuses to overwrite outputs. Input/output and code hashes are in `latest.json`.
