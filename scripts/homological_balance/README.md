# Homological-balance exact torus census

Ledger: `notes/homological-balance-root-ledger-20260906.md`
Verification: `notes/literature-officer-20260906-homological-balance-verify.md`
JSON: `results/homological-balance-exact-torus/latest.json`

```bash
python3 scripts/homological_balance/exact_torus_enum.py
python3 scripts/homological_balance/verify_independent.py
python3 -m unittest tests.test_homological_balance_exact_torus
```

`verify_independent.py` must be run with this directory on `sys.path` (it imports `exact_torus_enum`); running it as a file from this folder works.
