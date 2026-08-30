# Issue #154 Phase-E local-singlet pilot

The frozen 20k same-stream pilot completed on N65 and N130. The new rows are local
connectivity observables not determined by the occupation count K; they are not named as
continuum energy fields.

| N | P4[J_black] | P4[J_white] | P4[J_even] | P4[J_odd] |
|---:|---:|---:|---:|---:|
| 65 | -0.0004090712 +/- 0.00021 (z=-1.97) | 2.087674e-05 +/- 0.00014 (z=0.14) | -0.0001940972 +/- 0.00013 (z=-1.46) | -0.000214974 +/- 0.00012 (z=-1.80) |
| 130 | 1.777344e-05 +/- 0.00019 (z=0.10) | -0.0001694119 +/- 0.0001 (z=-1.65) | -7.581923e-05 +/- 0.00011 (z=-0.69) | 9.359266e-05 +/- 0.0001 (z=0.91) |

## Frozen common-ray comparison

| candidate | scale N130/N65 | chi2 / 2 df | p |
|---|---:|---:|---:|
| A/E/C | -1.7458 | 1.415 | 0.493 |
| A/E/J_even | -20 | 2.29 | 0.318 |
| A/E/J_odd | -2.3734 | 2.695 | 0.26 |

## Frozen gate

Decision: **stop_at_20k**. Qualifying fixed rows: `[]`.
Neither fixed local row resolves at both sizes, and neither improves the A/E/C
common-ray score. This pilot therefore gives no production evidence that the
radius-1 connectivity row replaces C or pins the E_top plane; the frozen extension
is not run.
All three candidates are reported; none was selected by its p value. N65 and N130
are one parent-child lineage and constitute one scale comparison, not independent
geometric votes. Full within-size batch covariance is retained in `score.json`.
