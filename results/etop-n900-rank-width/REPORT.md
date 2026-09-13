# N900: third-scale intrinsic rank-clock width

32M new shared counters, 800 aligned batches, two modulus pairs. This block is independent of N100/N400; predictions retain their common N400 anchor covariance.

Measured rank-step centered z variance: 2.339461729 +/- 0.120385.

| conditional prediction | expected Vz | observed minus expected | total SE | z | nominal p |
|---|---:|---:|---:|---:|---:|
| quarter_power_width | 2.56553539 | -0.22607366 | 0.15124 | -1.4948 | 0.134967 |
| fixed_critical_width_profile | 2.09475087 | 0.244710856 | 0.141704 | 1.72692 | 0.0841819 |

The two prediction comparisons share both the target and N400 anchor; they are not independent tests. The complete comparison covariance is saved, and no forced model ranking is reported.

Finite N400-to-N900 effective width: {'estimate': 0.30687685023405237, 'se': 0.038610153900379376, 'scope': 'Finite N400-to-N900 signed-profile effective width; not an asymptotic exponent.'}.

Both candidates are conditional finite-regime shape predictions selected from N100/N400, not established critical exponents. Include target sampling uncertainty and the perfectly correlated shared N400 anchor. Either, both, or neither may be incompatible; do not force a winner or call an unresolved comparison model recovery.
