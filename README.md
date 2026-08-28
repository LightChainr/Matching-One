# Matching One

Matching One is an open computational research program on the square-lattice site-percolation threshold and the exact relation to its matching lattice.

\[
p_c^{\mathrm{site}}(\mathbb Z^2)
+ p_c^{\mathrm{site}}(\mathrm{NN+NNN}) = 1.
\]

Numerical estimates place the square-site threshold near `0.59274605079`. The project does **not** treat a rounded estimate as a definition, and it does not claim a known closed form. Published high-precision estimates and finite-size sequences must be tracked as method-specific results with provenance; see #4.

## Scientific status

The project has moved beyond guessing combinations of familiar constants. Its main research direction is now the finite-size structure of the exact matching observable, including orientation dependence, matching parity, threshold-rank reconstruction, and possible spin-4 correction sectors.

Current evidence in the active research stack supports the following finite-size statements:

- a nonzero same-area, same-shape orientation effect on primitive Gaussian tori;
- the sign predicted by a `cos(4 theta)` harmonic;
- an independent-seed reproduction over five sizes;
- a held-out preference for a `cos(4 theta) N^(-13/8)` model over a zero effect in the current size range;
- local linear closure between the matching residual, its slope, and the finite-size root shift.

These observations do **not** yet establish a unique asymptotic exponent, an `x=21/4` logarithmic-CFT operator identification, an exact value of `p_c`, or a closed form. Clean-checkout replay and cross-size covariance are active P0 requirements in #39 and #46.

The repository-wide claim ledger is maintained in [`docs/STATUS.md`](docs/STATUS.md).

## Repository policy

`main` is the reviewed integration line. Research may proceed on short-lived branches, but code, protocols, and result archives become canonical only through a pull request to `main` with tests and provenance.

The existing research history is currently represented by a stacked sequence of pull requests:

1. #15 — executable research program and finite-size audits;
2. #18 — matching/orientation reference implementation and second-wave protocol;
3. #21 — server implementations, raw aggregates, and confirmation analyses;
4. #41 and #46 — post-confirmation decisions and covariance audit.

Large result archives are retained for auditability, but a committed result is not automatically an accepted scientific claim. Claim strength is governed by [`GOVERNANCE.md`](GOVERNANCE.md) and [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).

## Exact and reference quantities

| Quantity | Status |
|---|---|
| Square-site `p_c` | numerical, method-dependent; no known closed form |
| NN+NNN matching-site `p_c` | exactly `1 - p_c` |
| Square-bond `p_c` | exactly `1/2` |
| Triangular-site `p_c` | exactly `1/2` |

The machine-readable constants file records reference values and must not be used to turn rounded decimals into definitions.

## Layout

```text
constants/      reference values and provenance notes
notes/          theory, derivations, negative results, and research decisions
scripts/        reproducible analysis and reference checks
experiments/    frozen protocols and computation queues (research stack)
data/           source datasets with manifests and checksums (research stack)
results/        immutable result archives and reports (research stack)
tests/          regression and scientific-contract tests (research stack)
docs/           governance, status, and roadmap
```

## Quick check

The small script currently present on `main` can be run with:

```bash
python3 scripts/compare_candidates.py
```

Research branches containing the full toolchain should be run from a clean checkout using the dependency and command records committed with each result. Pull requests to `main` are checked by GitHub Actions across the supported Python versions.

## Working principles

1. Prefer lattice-native identities, exact controls, and falsifiable predictions over numerology.
2. Freeze models, seeds, geometry sets, and held-out data before scoring.
3. Preserve failed models, null results, raw sufficient statistics, and provenance.
4. Separate empirical finite-size observations from asymptotic or operator-level interpretation.
5. Use expensive hardware only after a CPU oracle, power calculation, and information-per-wall-time benchmark.

See [`CONTRIBUTING.md`](CONTRIBUTING.md) before opening a research pull request.

## License

MIT. See [`LICENSE`](LICENSE).
