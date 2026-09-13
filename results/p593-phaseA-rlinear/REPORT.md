# p593 Phase A — r_linear definition pinned (raw: candidates.json)

- date: 2026-09-13
- issue: 593 (Phase A)
- question: which candidate definition of `r_linear` reproduces 10/26/72 at widths 4..6?
- answer: the joint Krylov dimension `dim span{G0^k f : k >= 0, f in D0 readouts}` under exact elimination mod 2147483647; candidates 2 (lumping) and 3 (pencil ranks) are falsified.
- new value: true `r_linear` at width 8 = 689 (repository had budget-truncated `>= 150`).

## Commands

```bash
cd scripts
python p593_phaseA_rlinear_candidates.py
```

## Environment

- Python 3.13 (managed env), numpy 2.5.2, scipy 1.18.1
- repository generator/codec code unchanged; all arithmetic exact over Z/2147483647
- no frozen configuration was substituted: seeds, readouts, dictionary nesting and the prime are the repository's own declarations

## Verdict per candidate

| candidate | definition | verdict |
|---|---|---|
| 1 | joint Krylov over D0 seeds, exact mod-p | **matches 10/26/72/218 exactly** |
| 2 | linear lumping (coarsest exact strong lumping block count) | gives 76/232/750 at widths 6..8 = `r_positive`, not `r_linear` |
| 3 | pencil ranks (G0, H, J, D mod p; row patterns; value tuples) | no match at any width |

Full Matching-One repository CI has not been run for this commit.
