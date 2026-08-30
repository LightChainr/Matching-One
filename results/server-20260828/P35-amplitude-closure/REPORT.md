# P35 same-batch amplitude closure

The closure statistics use the Huawei ten-million-per-size threshold-rank
histograms from P33.  `M`, both analytic slopes, both direct roots and every
delete-one-batch replicate are reconstructed from the same integer sufficient
statistics at `p=0.592746050790`.

| N | C = -Delta p* mean(M') / Delta M | jackknife SE | direct gap | linearized gap |
|---:|---:|---:|---:|---:|
| 65 | 1.0003073 | 1.41e-4 | -8.55350e-5 | -8.55087e-5 |
| 85 | 0.9998355 | 1.06e-4 | -7.88704e-5 | -7.88834e-5 |
| 130 | 1.0000003 | 1.96e-5 | -4.14783e-5 | -4.14783e-5 |
| 145 | 1.0000077 | 7.90e-5 | -5.66827e-5 | -5.66822e-5 |
| 170 | 0.9999611 | 5.61e-5 | -4.86397e-5 | -4.86416e-5 |

The linearized relation closes at the `3e-4` level or better for every size;
the root lies well inside the local linear regime.  The slope amplitude
`B=N^-3/8 mean(M')` is also stable, decreasing only from 1.7514 to 1.7462.
Thus the mechanism linking a central matching residual to the root shift is
confirmed for these finite curves.

This does not by itself establish the radial exponents: the ten-million P33
pilot has sizeable Monte Carlo drift in `A_M` and `A_p`, and its held-out
constant-`N^2 Delta p` check fails the frozen threshold.  The result isolates
that remaining failure to radial amplitude/correction modeling rather than to
nonlinearity of the root conversion.  Full precision values and all jackknife
inputs remain embedded in `P33/p33_all_sizes_arm_10m.analysis.json`.
