# Matching One

Matching One is an open computational research program on the square-lattice site-percolation threshold and the exact relation to its matching lattice.

\[
p_c^{\mathrm{site}}(\mathbb Z^2)
+ p_c^{\mathrm{site}}(\mathrm{NN+NNN}) = 1.
\]

Numerical estimates place the square-site threshold near `0.59274605079`. The project does **not** treat a rounded estimate as a definition, and it does not claim a known closed form. Published high-precision estimates and finite-size sequences are tracked as method-specific results with provenance; canonical reconciliation is issue #4.

## Scientific direction

The project has moved beyond guessing combinations of familiar constants. Its main research direction is the finite-size structure of the exact matching observable, including orientation dependence, matching parity, threshold-rank reconstruction, and possible spin-4 correction sectors.

The current server/research archive contains finite-size evidence for:

- a nonzero same-area, same-shape orientation effect on primitive Gaussian tori;
- the sign expected from the tested `cos(4 theta)` design;
- independent-seed reproduction over five prescribed sizes;
- held-out support for a `cos(4 theta) N^(-13/8)` law over the current size range;
- local linear closure between the matching residual, its slope, and the finite-size root shift;
- two prospective Gaussian `1+i` lineages compatible with the no-fit doubling ratio `-2^(-13/8)`.

These observations do **not** establish a unique asymptotic exponent, unique H4 harmonic, an `x=21/4` logarithmic-CFT operator identification, an exact value of `p_c`, or a closed form. Clean-source replay, covariance hardening, canonical archive import, prospective harmonic/root tests, and exact parity controls remain active gates.

The repository-wide claim ledger is [`docs/STATUS.md`](docs/STATUS.md).

## Canonical repository policy

`main` is the reviewed integration line. It currently contains the governance/CI baseline and the foundational research methods integrated from PR #15 and PR #18.

A research branch, issue, server directory, or pull-request description is evidence and working history; it becomes canonical only through reviewable integration to `main` under [`GOVERNANCE.md`](GOVERNANCE.md) and [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).

The large Huawei server campaign is intentionally preserved in PR #21 / `server/huawei-analysis-20260828`. It is **not** to be bulk-merged into `main`. Issue #59 governs a curated import of production source/tests, archive manifests, and bounded immutable result families. PR #46 (cross-size covariance infrastructure) and PR #56 (post-doubling decision layer) are currently drafts pending their P0 dependencies.

GitHub Actions checks `main` and pull requests across Python 3.9/3.11/3.13 plus C++17 build/self-tests. Hosting-side protection of `main` is still tracked by issue #52.

## Exact and reference quantities

| Quantity | Status |
|---|---|
| Square-site `p_c` | numerical, method-dependent; no known closed form |
| NN+NNN matching-site `p_c` | exactly `1 - p_c` |
| Square-bond `p_c` | exactly `1/2` |
| Triangular-site `p_c` | exactly `1/2` |

`constants/pc.yaml` records reference values; until #4 is complete it must not be interpreted as an adjudication of disputed last digits. Broad PSLQ/constant search is P2 and blocked by #4 under issue #1.

## Layout

```text
constants/      reference values and provenance notes
data/           canonical source datasets when integrated
notes/          theory, derivations, negative results, and research decisions
scripts/        reproducible analysis and reference checks
experiments/    frozen protocols and computation queues
predictions/    preregistered prediction artifacts
results/        immutable result archives when curated onto main
tests/          regression and scientific-contract tests
docs/           governance, status, and roadmap
```

Not every directory present on a research/server branch is already canonical on `main`; consult `docs/STATUS.md` and issue #59 for integration status.

## Local checks

Install the declared Python dependencies when needed, then run the integrated unit/regression suite:

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

A minimal numerical smoke check is:

```bash
python3 scripts/compare_candidates.py
```

Production result families must additionally follow the clean-checkout, source/binary hash, RNG-domain, batch, checksum, and raw-sufficient-statistic requirements in `REPRODUCIBILITY.md`.

## Working principles

1. Prefer lattice-native identities, exact controls, and falsifiable predictions over numerology.
2. Freeze models, geometry/order conventions, RNG domains, and held-out data before scoring.
3. Preserve failed models, null results, raw sufficient statistics, and provenance.
4. Separate finite-size observations from asymptotic, universal, or operator-level interpretation.
5. Treat topology, RNG, covariance, and threshold-rank code as high-risk and require independent/exact checks.
6. Use expensive hardware only after a CPU oracle, power calculation, and information-per-wall-time benchmark.

See [`CONTRIBUTING.md`](CONTRIBUTING.md) before opening a research pull request.

## License

MIT. See [`LICENSE`](LICENSE).