# Same-prefix exact completion-pair refinement

All 140000 rows aligned across 14 original batches; every first-five-field tuple matched.

Early minus late; fixed a,b and original overlap weights. No new random prefixes.

| Metric | Estimate | Batch jackknife SE |
|---|---:|---:|
| RB_deltaY | -0.00344938566765506 | 0.00077225276269743 |
| original_probe_deltaY | 0.00668182660775089 | 0.0120870551771217 |
| paired_RB_minus_probe | -0.0101312122754059 | 0.0123956612572906 |
| RB_deltaE | -183.963984858095 | 41.2052909135037 |
| RB_deltaZ2 | 3.22818100446661e-08 | 7.23066188814621e-09 |
| RB_earlyY | 0.0955745614429761 | 0.000450216499320379 |
| RB_lateY | 0.0990239471106312 | 0.000525721971141274 |
| secondary_deltaNu | -0.0398651161631196 | 0.569877945474596 |

Descriptive SE_RB/SE_probe: 0.0638908941326869. Variance ratio: 0.004082046353074206.

| Support / cohort | Rank-one | Eligible | Supported |
|---|---:|---:|---:|
| primary / early | 24201 | 24201 | 24109 |
| primary / late | 30827 | 30827 | 30354 |
| primary / combined | 55028 | 55028 | 54463 |
| secondary / early | 24201 | 24201 | 24201 |
| secondary / late | 30827 | 30827 | 30821 |
| secondary / combined | 55028 | 55028 | 55022 |

Original committed primary and secondary point/SE reproduced (absolute tolerances 1e-14 / 1e-12).
This is a source/estimand identity check, not scientific replication.

JSON retains exact integer cell sums, mean_e, meanY, all 14 aligned deletion vectors/supports, full 8x8 covariance, and input/source hashes and commits. First seven metrics share exact (D,c) cells; the secondary uses all rank-one prefixes within D, including c=m.

Y=2e/(m-c); original prefix probe mean=(Y0+Y1+Y2+Y3)/4. RB_deltaE averages cell mean-edge differences; RB_deltaZ2=-2*RB_deltaE/[m(m-1)]. c=m has no conditional safe mean. Singular covariance is expected; no inverse is used.

- Estimator refinement after the original readout; same 140k prefixes and batch blocks, zero new samples.
- No independent confirmation: reproducing old primary/secondary is source/estimand identity, not scientific replication.
- Exact means conditional mean 2e/(m-c); finite-prefix uncertainty remains. Four actual old probes are averaged per prefix.
- SE_RB/SE_probe is descriptive; variance reduction does not guarantee empirical delete-batch jackknife SE shrinkage.
- No automatic Markov claim: one weighted moment can cancel across cells; finite L512 count clock only, no continuum claim.
