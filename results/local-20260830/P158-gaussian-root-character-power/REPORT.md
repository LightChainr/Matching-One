# P158 Gaussian commuting-square power gate

No target production was run.  Source-only forecast:

- `DeltaM` factor: `0.0132796927517053`;
- slope factor: `2.37137370566166`;
- root factor: `7/1250 = 0.0056`;
- N650 leading target: `DeltaM=+1.50913e-5`, root gap `-7.59549e-7`;
- N850 leading target: `DeltaM=+1.07152e-5`, root gap `-4.88096e-7`;
- per-lineage samples for two-sided `alpha=.01`, 80% power versus zero:
  `27.21B` and `55.18B`;
- 500M expected z: `0.470` and `0.337`;
- even perfect fresh parent-child CRN can reduce variance by at most `2.15%`
  and `2.34%`.

The integer-period runner already accepts both target matrix pairs.  The
existing root-amplitude scorer must be generalized beyond its hard-coded
`[65,85]` pair before a target pilot.

Machine-readable details: `power.json`.

Validation:

```text
python3 tests/test_gaussian_root_character_power.py -v       # 5 passed
python3 tests/test_gaussian_commuting_square_root.py -v      # 5 passed
threshold_rank_integer_period_mc --self-test                 # passed
100-replica custom N650/N850 CLI smoke                        # passed
```

