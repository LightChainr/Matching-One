### Science card — Phase-E local-singlet mixed plane

Branch `analysis/p154-phase-e-mixed-plane-pilot-20260830`; freeze `0578105`,
result commit pending. The same-stream pilot records the matching-even local
connectivity row `B` together with `B^2`, `I0*B`, and `I2*B`, making
`J_top=Cov(I2-I0,B)` and `J_bulk=Cov(I2+I0,B)` directly scoreable rather than
inferring energy from the topological `E_top` coordinate.

At N65/N130, `P4[J_bulk]` is respectively `7.82e-5 +/- 1.07e-4` and
`-2.64e-5 +/- 6.47e-5`; it is unresolved and changes sign. Replacing C by
J_bulk lowers the two-size common-plane chi-square from `2.2740/2` to
`.29264/2`, but the improvement `1.9814` misses the frozen threshold 4 and
J_bulk is not resolved at both sizes. `J_top` is negative at both sizes but
only `2.08` and `1.39` standard errors from zero.

Mechanism consequence: this explicit radius-one matching-even mixed local
singlet does not pin or replace C at 20k and stops without extension. This is
not a no-go for every local energy definition and does not identify a continuum
field. A future Phase-E restart needs a genuinely different microscopic
singlet, not more samples of this row.

Huawei provenance: `DevEnvC_ZyTrST` N65 and `DevEnvC_XPk2PZ` N130, each 20k / 100
batches; remote/local hashes agree.
