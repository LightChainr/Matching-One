# P250 N505 existing-data state increment

This analysis adds no samples. It cross-fits the archived 200k/400-batch
defined-event table, with whole batches (and therefore shared replicas and
both child views) kept together.

## Result

- M_spec held-out MSE: `0.70172768`.
- M_state held-out MSE: `0.48438168`.
- Row-weighted relative gain: `30.973%`.
- Equal-batch loss reduction (M_spec - M_state): `0.21736439 +/- 0.0126`; one-sided cluster t p=`1.192e-50`.
- Frozen decision: **pre-outcome antisymmetric state block adds held-out information**.

| hand | rows | M_spec MSE | M_state MSE | relative gain | batch p |
|---|---:|---:|---:|---:|---:|
| plus | 2653 | 0.7082069 | 0.47543481 | 32.868% | 6.629e-36 |
| minus | 2698 | 0.69535653 | 0.49317934 | 29.075% | 3.697e-33 |

## What the two models contain

M_spec is fit separately in each hand. In every outer fold it selects
two residue Fourier frequencies using training rows only, then fits those
four sine/cosine coordinates together with orientation and fibre controls.
M_state uses the identical selected frequencies and adds only Sminus:
the antisymmetric change of rank, canonical basis class, support site phase,
and component/site-change indicators between the two ordered intermediate
supports. No L or R coordinate is a predictor.

Selected training-only frequencies by fold:

- plus: k=1 (10/10 folds), k=9 (2/10 folds), k=10 (1/10 folds), k=15 (1/10 folds), k=36 (2/10 folds), k=43 (1/10 folds), k=46 (3/10 folds)
- minus: k=1 (10/10 folds), k=10 (10/10 folds)

## Boundary

The archive contains neither endpoint full fields nor endpoint periodograms. Delta Pminus=P_DJ-P_JD cannot be reconstructed. M_spec below is only a low-dimensional conditional-mean model for Rminus over sampled nonzero residues; it is not a terminal spectral-response test and cannot repair the unobserved residue-zero completion.

A held-out gain would establish incremental predictive information in the retained pre-outcome typed support state for this finite adaptive intervention. It would not identify a CFT field, count physical states, reconstruct Delta Pminus, prove spontaneous path memory, or exclude alternative nonlinear/gauge descriptions.
