# P154 Phase E: topology versus ordinary bulk-energy coordinate

## Answer

The locked PR #273/#277 production archives do contain the ordinary
**integrated Bernoulli-energy** insertion, because every batch retains the
complete microcanonical threshold histograms.  But it is not a second
`J_bulk` direction.  It is exactly the already-scored topology-even birth
coordinate:

```text
J_bulk,integrated
  = sqrt(N/[p0(1-p0)]) P4[S_mode 1]
  = P4[S_historical']
  = (1/2) P4[E_top'].
```

After the existing P154 normalization this is exactly the scalar `U` used in
the frozen norm-4 score.  Across N=65,130,260,520 and N=85,170,340,680:

```text
max |J_bulk,integrated - J_top,even-birth|       = 1.04e-12
max |alias residual| in 16 explicit delete-one rows = 9.20e-12
max |reconstructed U - committed U|             = 4.44e-16
```

The full 16 by 16 same-stream covariance is archived in `latest.json` in the
order `(J_top,J_bulk)` at each size.  It is rank-duplicated by the exact
identity, not evidence for two correlated fields.

## Consequence for the proposed Phase-E comparison

The proposed `J_top` versus ordinary integrated `J_bulk` vote is
**structurally non-identifiable on these archives**: the two labels name the
same production coordinate.  The borderline PR #273 Jordan score therefore
already is the topology-even integrated-energy score; recomputing it under a
bulk label cannot absorb an additional residual or support a new bulk
primary.

This is useful model elimination.  It removes the apparent zero-cost second
direction without spending another sample.  A genuinely independent
`J_bulk` must be a separately declared microscopic local Potts/cluster
singlet, not another derivative or Krawtchouk view of the threshold curve.

## Archive audit

The authoritative blocks are:

- PR #273, commit `8b26a30`: N65/85/130/170/260/340 production;
- PR #277, commit `3e855ce`: N520/N680 generation-four blocks.

Their only production columns are marginal K1/K2 histograms and the integer
moments

```text
sum K1, sum K2, sum K1^2, sum K2^2, sum K1*K2,
sum gap, sum gap^2.
```

These are sufficient for every integrated Bernoulli score mode and the full
batch-estimator covariance.  They contain no local pair-connectivity,
cluster-energy, cluster-count, or event-times-local-field cross moment.

## Minimum future statistic, only if a local singlet is wanted

Freeze a fixed `p_ref`, a D4-scalar translation-averaged local singlet
`J_local`, its microscopic cutoff, and its centering.  Per batch and
orientation the minimal sufficient sums are

```text
samples,
sum J_local, sum J_local^2,
sum I_rank0, sum I_rank2,
sum I_rank0*J_local, sum I_rank2*J_local.
```

Those sums permit connected rank-sector couplings and their full batch
covariance.  No pilot was started here because the requested ordinary
integrated energy statistic was already present and exactly aliased; choosing
a new microscopic `J_local` would be a new observable contract, not a missing
adapter field.

## Scientific card

- **Mechanism space changed:** topology-even companion versus an allegedly
  distinct ordinary integrated bulk-energy direction.
- **Result:** the two are one exact coordinate; the distinct integrated-bulk
  alternative is eliminated as an identification, not statistically fitted.
- **Not proved:** an independent microscopic local singlet is neither present
  nor excluded.
- **Observer / sector / source / geometry:** normalized `P4[E_top']/2` /
  Alexander-even birth sector / Bernoulli occupation-score insertion /
  P154 norm-4 Gaussian lineages.
- **Dependency:** all rows are exact re-expressions of PR #273/#277 and must
  not be counted as new evidence.
- **Natural upweight:** only a predeclared local-singlet mark with event cross
  moments can test a genuinely new bulk direction.

## Reproduction

```bash
python3 scripts/score_p154_phase_e_bulk_alias.py \
  --output results/p154-phase-e-bulk-energy-alias/latest.json
python3 -m unittest discover -s tests \
  -p 'test_p154_phase_e_bulk_alias.py'
```

No new Monte Carlo sample was generated and no Huawei environment was used.
