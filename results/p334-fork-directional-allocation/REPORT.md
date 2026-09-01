# P334 fork directional allocation: nonclosure is not yet an H4 mechanism

**The within-checkpoint degree-dispersion excess has a resolved positive mean, but its directional contrast is unresolved in this design.** This zero-new-sample result projects the parent's complete checkpoint covariance; it does not rerun the already completed cooperative decomposition.

The exact real-checkpoint counterexamples remain decisive against the specified scalar state `(N, orientation, k0, H2, b2, age, ell)`. That is a different claim from identifying which coordinate carries an orientation response. No new significance threshold, exponent fit, Monte Carlo block, or scalar-state test is introduced here.

## Definitions and dependence

Let `d=N-k0`, `b1=d-H2`, `c_v` be the safe second-site count after safe first site v, and `b2=sum(c_v)/2`. Per checkpoint, `s1=b1/d`, `s2=2*b2/[d(d-1)]`, and the one-common-update/two-clone probability is `sum(c_v²)/[d(d-1)²]`. Capital `S1,S2` denote checkpoint means, `Q` the mean exact branching probability, and `B=E[s2²/s1]`.

`total_gap = Q-S2² = common_gate + between_checkpoints + within_checkpoint`, where the three parts are `(1/S1-1)S2²`, `B-S2²/S1`, and `Q-B=E[s1 Var(c_v/(d-1) | safe,C)]`. `Q` itself is separately reported and must not be confused with `total_gap`.

The source has 20k base permutations per size. A permutation, both orientations and all clone/exact-count rows form one checkpoint cluster. The input 22×22 covariance is propagated through the parent's differentiable aggregate definitions and then the linear `second-first` contrast. The two sizes remain separate; all output rows reuse their original dependency groups, not additional independent evidence.

## Direction-resolved outputs

All ± quantities below are parent checkpoint-cluster/delta-method SEs. Standardized contrasts are descriptive post-reveal summaries, not multiplicity-adjusted claims. Cos4 normalization is a geometric comparison only: with two directions it cannot distinguish H4 from other angular contributions.

### N325

Physical pair: `[17, 6]` → `[18, 1]`. Exact Δcos4 = `16128/21125`.

| observable | first ± SE | second ± SE | second−first ± SE | contrast / Δcos4 ± SE |
|---|---:|---:|---:|---:|
| branch_success_Q | 0.82078417 ± 0.0010279 | 0.81814486 ± 0.0010324 | -0.0026393091 ± 0.0014544 | -0.0034570563 ± 0.001905 |
| total_gap | 0.052667722 ± 0.00027354 | 0.053335486 ± 0.0002723 | 0.0006677647 ± 0.00038529 | 0.00087466079 ± 0.00050466 |
| common_gate | 0.051255415 ± 0.00026719 | 0.051912187 ± 0.00026669 | 0.00065677207 ± 0.00037691 | 0.00086026228 ± 0.00049369 |
| between_checkpoints | 0.0012934924 ± 1.9681e-05 | 0.0013008689 ± 1.8932e-05 | 7.3765566e-06 ± 2.7332e-05 | 9.6620634e-06 ± 3.5801e-05 |
| within_checkpoint | 0.00011881426 ± 2.1576e-06 | 0.00012243033 ± 2.4034e-06 | 3.6160695e-06 ± 3.2507e-06 | 4.7364501e-06 ± 4.2579e-06 |

The within-checkpoint contrast / SE is `1.112`. No component/total-directional-gap fractions are formed: signed, uncertain directional denominators would turn additive responses into unstable mechanism shares.

### N425

Physical pair: `[16, 13]` → `[19, 8]`. Exact Δcos4 = `32256/36125`.

| observable | first ± SE | second ± SE | second−first ± SE | contrast / Δcos4 ± SE |
|---|---:|---:|---:|---:|
| branch_success_Q | 0.84822319 ± 0.00090309 | 0.84956895 ± 0.00089641 | 0.0013457601 ± 0.0012665 | 0.0015071795 ± 0.0014184 |
| total_gap | 0.045569828 ± 0.00025065 | 0.045235406 ± 0.00024851 | -0.00033442221 ± 0.00035142 | -0.00037453505 ± 0.00039358 |
| common_gate | 0.04453738 ± 0.00024587 | 0.044202929 ± 0.00024416 | -0.00033445141 ± 0.00034496 | -0.00037456774 ± 0.00038633 |
| between_checkpoints | 0.00095547056 ± 1.4636e-05 | 0.00095703464 ± 1.3961e-05 | 1.5640862e-06 ± 2.0092e-05 | 1.7516931e-06 ± 2.2502e-05 |
| within_checkpoint | 7.697726e-05 ± 1.5008e-06 | 7.5442368e-05 ± 1.4499e-06 | -1.5348918e-06 ± 2.0863e-06 | -1.718997e-06 ± 2.3366e-06 |

The within-checkpoint contrast / SE is `-0.736`. No component/total-directional-gap fractions are formed: signed, uncertain directional denominators would turn additive responses into unstable mechanism shares.

## What changes next

1. **Keep the exact nonclosure result; do not resample its already saved inputs.** The safe-insertion graph's degree second moment (equivalently its 2-star overlap count) distinguishes real checkpoints with the same scalar tuple. It is a concrete microscopic coordinate, not an unexplained generic memory label.
2. **Do not call successor-H2 one-step closure a new mechanism experiment.** In the new archive, on `branch_common_safe=1`, `H2_after=(N-k0-1)-branch_q_after_safe_count`. Given this value, the independent one-site clones have exact success `q=1-H2_after/(N-k0-1)` and product expectation `q²`; this is calibration. Absorbed rows use q=0 and do not define a rank-one successor H2.
3. **Use existing states to target a genuinely richer question.** Frozen seed/counter/k0/period-matrix metadata reconstructs the current microscopic configuration. A deterministic replay can expose the safe-insertion degree distribution and ask what boundary organization carries its 2-star variation; a specifically chosen third-clone fan-out reads the cubic degree moment. These are possible follow-ups on existing states, not requests for another generic production block. The present direction contrast supplies no identified H4 carrier and does not justify enlarging the same experiment by default.

## Provenance, reproducibility and narrow check

Source branch `experiment/p334-cooperative-closure-pilot-20260830` at `6147e22f53902a94e5f133739f2c1d423691d0b8`; production result `e81dd59ff6be69056e504e0e81cfeccf73dc5e97`. Input score SHA256 `ba51551cecf4feb8e48c7c95f105c38a744476eb625c30504fe9574e176f8d11`. Parent exact witnesses: `notes/p334-real-checkpoint-scalar-nonclosure.md`. This is a retrospective derived analysis, not independent confirmation.

`score.json` preserves the parent 22×22 covariance, Jacobians, both-direction 10×10 covariance, 5×5 raw and normalized contrast covariances, exact rational geometry, environment versions and the additive check. Only one additive/parent-projection check is run; no test suite and no Monte Carlo are invoked.

```bash
/Users/lc/python-envs/research-py311/bin/python scripts/analyze_p334_fork_directional_allocation.py
```

Assessment: share with the stated finite-design, post-reveal and angular-identification caveats. No full-state temporal memory, scaling law, or continuum field identity follows from this projection.
