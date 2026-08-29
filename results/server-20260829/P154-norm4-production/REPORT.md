# P154 norm-4 production reveal

Issue #154's frozen dyadic norm-4 experiment is complete. Each source size
combines the pre-existing aligned 100M block with a new aligned 1.9B block;
N260 and N340 each contribute an independent 1B target block. All new blocks
use generation commit `bfab0330f5f56ca4d746b45d737f1607e3d229a0` and the
counter intervals frozen in `experiments/norm4_variance_pilot_20260829.yaml`.

## Frozen result

| score | chi-square / df | p-value | decision |
|---|---:|---:|---|
| analytic q=2 scalar | 20.897 / 2 | 2.899e-5 | rejected |
| rank-2 Jordan scalar | 5.397 / 2 | 0.06731 | survives at 5%, with tension |
| q=2 common thermal-jet generator | 65.040 / 10 | 3.982e-10 | rejected |
| Jordan common thermal-jet generator | 18.043 / 10 | 0.05425 | survives at 5%, with tension |
| central-D 4-to-1 character | 4.117 / 2 | 0.1277 | not rejected; correlated secondary |
| root +1/16 character | 4.130 / 2 | 0.1268 | not rejected; correlated secondary |

The scalar q=2 residuals are `0.27250, 0.69306`; the corresponding Jordan
residuals are `0.01079, 0.21007`. The second Jordan lineage is the visible
scalar tension (`2.32` marginal standard errors). The target values are
`U260=1.66079` and `U340=1.98451`.

The scalar and thermal-jet rows reuse the same curves and are not additive
evidence. Their agreement is nevertheless structurally informative: the
fixed q=2 multiplier fails both the scalar projection and the single common
ranks-2-through-6 jet, while the fixed Jordan multiplier lands just outside
the conventional 5% rejection boundary in both views.

## Mechanism inference

The ordinary analytic `A+C/N` closure is no longer a live leading model for
these two noncyclic dyadic lineages. A rank-2 Jordan direction captures the
dominant radial transfer far better, and its survival across the full frozen
jet makes the scalar result unlikely to be a one-coordinate accident.

The remaining lineage-dependent tension says that a *scalar* Jordan law is
not the whole transfer operator. The most economical next model is a
two-dimensional conjugation-even transfer: one Jordan generalized eigenvector
plus one subleading even mode, fitted on revealed source generations and then
held out on a new Gaussian multiplier. This is a mechanism inference, not a
claim that the present data identify a unique 2x2 matrix.

## Exploratory next prediction

Do not rescue q=2 with a free scalar exponent. Freeze a source-only 2x2
even-generator transfer, then use norm 10 as a quadrature target. The pair
`2+i -> 3+i` has the already-recorded invertible spin-4 phase matrix, so it can
separate the conjugation-even correction exposed here from the odd generator
that vanishes at the norm-4 phase node. A successful held-out norm-10 matrix
composition would upgrade the current "Jordan plus one missing mode"
interpretation into a genuine Gaussian-semigroup mechanism; failure would
falsify that economical completion.

## Execution and provenance

N260/N340 and N65/N130 were generated on Huawei DevEnv
`f415a4bcbd9a438b85f5f29e4a507ea4`. N85/N170 were generated concurrently on
DevEnv `f550f3cb1f774374b6842aa648fda796` with the exact same frozen aarch64
binary (SHA256 recorded in `environment.txt`). The redundant N85/N170 leg of
the first host's sequential driver was stopped after N65/N130 completed; its
exit-143 line records that intentional cancellation, not a failure of any
accepted block. The accepted N85/N170 run has empty stderr and status `0,0`.

Every raw and derived artifact is covered by `checksums.sha256`. Exact scorer
commands are in `commands.txt`; machine-readable results are
`analysis/scalar_score.json` and `analysis/thermal_jet_score.json`.
