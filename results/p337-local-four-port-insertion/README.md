# Fixed-origin four-port insertion: raw cross-moments

The two complete N25 outputs were produced once each after explicit GO
under contract `d7f15e68d593719e274fe9f0d5c7841a553beade`, using unchanged
producer `97f3efee905ed2d634d1be5245efe2a9cfe070bd`. Both accepted theory
gates, full hashes and build/run commands are recorded in [run.json](run.json).

Each file has K=0,...,25 and columns
`k,count,sum_q,sum_e,sum_s2,sum_qs2,sum_es2`.
The final three columns use the integer source `s2=-2*t_origin`; divide
them by2 exactly once. The full2^25 population is retained, including
origin-occupied configurations whose source is zero. There is no added
binomial multiplicity and no factor25 in these site-average source units.

Origin vacancy, occupied ports, exact two-component pairings and the
C4-averaged weights follow the frozen contract. The original black NN,
white matching/Alexander rank and q/E definitions are unchanged. The
rotation average is not asserted to be one fixed-cut idempotent projector.

Axis runtime:1.10562s; tilted runtime:1.03976s, as concurrent local
single-thread enumerations. The root coordinator owns the single full
covariance/moving-root response score. No score, root search, old-source
readout, extra statistic, MC, cloud task or test suite was run here.
