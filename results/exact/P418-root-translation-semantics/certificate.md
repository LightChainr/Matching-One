# P418 root-translation semantic certificate

No new Monte Carlo was run. Four deterministic nonzero period-loop witnesses per N505 child were translated through all 101 parent positions.

| gate | failures |
|---|---:|
| root component signature | 0 |
| real root scalar | 0 |
| BFS/CRT section and stored gauge | 0 |
| full translation-orbit factorization | 0 |
| historical archive section/gauge provenance | 0 |

There are `16` exact counterexamples to pulling the mask through a *single fixed configuration*. The first is `plus/black_period_column_0/r1/d=27` with canonical Phi5 gap `[492, 72, -372, -292]`.

The radius4/radius5/radius6 archives were produced at different runner commits, but all six files that define the section, gauge, cover, root observable and DFT have identical SHA-256 hashes at those commits; the radius4 stored gauge array also equals the reconstruction exactly. This localizes the Issue 418 failure away from root selection, away from the CRT section/gauge, and away from a historical old4 gauge mismatch. The mask factorization is restored exactly on the full 101 configuration translations times 101 anchors. The remaining target is the archive one-anchor/covariance/scorer assembly.

Boundary: the fixed-configuration counterexample is expected and is not evidence against stationarity; this gate does not decide sampling variance versus covariance/scorer assembly.
