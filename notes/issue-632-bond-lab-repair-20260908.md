# #632 bond laboratory repair (2026-09-08)

Local exact repair. No new L, no MC.

## Old → withdrawn / corrected / kept

| Claim | Status |
|---|---|
| Bond `dual_fail=118133/262144` | **Withdrawn as physics.** Compound of forest sign `(ax+dx)` and dual occupancy `not occ[pi]`. |
| Geometric dual + corrected rank | **Corrected.** `dual_fail=0`; pairs `(0,2)/(1,1)/(2,0)=75460/111224/75460`; `M(1/2)=0`, `F(1/2)=1/2`. |
| `eval_F_float` site L=3 at 1/2 = 115.936… | **Withdrawn.** Coefficients are already counts; extra binomial removed. Now matches `43/128`. |
| Same-M / different-F via `P11→3/2 P11` | **Withdrawn as a simplex counterexample.** Mass `593/512`; after normalization M changes. Kept: equal M need not determine the full joint `(P11,P20,P02)`. |

Contractible plaquette `[(0,1,1,0),(1,4,0,1),(3,4,1,0),(0,3,0,1)]` has rank 0 after the sign fix (was 1).

Repro:

```bash
python3 scripts/probe_invariant_shape/exact_rank_census.py --bond --check-606
python3 tests/test_probe_invariant_shape.py
```

Does not edit `docs/STATUS.md`. Does not close #628/#631/#635.
