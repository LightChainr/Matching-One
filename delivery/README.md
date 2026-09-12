# Matching One multi-advance handoff — 2026-09-13

## New completed results

Read `delivery/notes/multi-advance-review-20260913-zh.md` first. The six technical notes
cover (i) full spatial-source Hessians, (ii) invisible positive spatial marks,
(iii) annealed versus quenched roots, (iv) the rank-source zero-free strip and
extensive-source phases, (v) geometrically anchored conditional-odds integration,
and (vi) continuum/lattice limit order. There is no novelty claim or new production.

`delivery/` contains only NEW repository files. The patch is additive.

## Prerequisites

Python 3.10+ and mpmath. Exact enumeration/polynomial/sign components use only the
standard library; mpmath is used for explicitly labelled numerical roots/integrals.
The complete test set and independent comparison also require these existing files:

- `results/research-control-20260912/width4-rank-closure-certificate.json`
  from PR #708, Git blob `50b7297deefe7c50215aea2ed534ca5810461af3`.
- `results/research-control-20260912/width4-parametric-definition.json`
  from PR #710, SHA-256
  `5edc624377399ad0878c7a608e722ebb65a7be454a6f022a1b6382a0661e91b0`.

They are supplied unchanged in `source_inputs/` solely for reproduction, not duplicated
by the new patch. A tree based on #710 has both. Neither is assumed present on main.

## Reproduce after applying the patch

```sh
git apply --check matching-one-multi-advance-20260913.patch
git apply matching-one-multi-advance-20260913.patch
python -m unittest discover -s tests -p 'test_spatial*.py' -v
python -m unittest discover -s tests -p 'test_rank_source_zeros.py' -v
python scripts/torus_source_hessian.py --certificate results/research-control-20260912/width4-rank-closure-certificate.json --out /tmp/full-site-hessian-new.json
python scripts/verify_source_hessian.py --certificate results/research-control-20260912/width4-rank-closure-certificate.json --hessian results/research-control-20260912/full-site-hessian.json --out /tmp/hessian-check-new.json
python scripts/invisible_spatial_marks.py --certificate results/research-control-20260912/width4-rank-closure-certificate.json --hessian results/research-control-20260912/full-site-hessian.json --out /tmp/marks-new.json
python scripts/topological_source_zeros.py --definition results/research-control-20260912/width4-parametric-definition.json --out /tmp/zeros-new.json
python scripts/conditional_odds_integration.py --hessian results/research-control-20260912/full-site-hessian.json --out /tmp/integration-new.json
```

Scripts refuse to overwrite existing result paths. Use fresh output names on reruns.

## Validation boundaries

27 focused tests passed. The patch was applied in a clean MINIMAL Git tree containing
the two source JSONs; all newly applied bytes were compared; all five result generators
were rerun. This is not a checkout of the full repository and not full repository CI.
`VALIDATION.json` and `HIGH_PRECISION.json` record checks actually executed.
`REMOTE_ACTIONS.json` records two posted PR conversation comments only. No branch,
issue lifecycle, PR merge, or STATUS mutation was performed.
