# Matching One

Matching One is a computational research project on square-lattice site percolation, its matching-lattice identity, and the finite-size structure behind the threshold.

The exact structural anchor is

\[
p_c^{\mathrm{site}}(\mathbb Z^2)+p_c^{\mathrm{site}}(\mathrm{NN+NNN})=1.
\]

## Current empirical picture

The strongest numerical structure is a matching-odd orientation signal on primitive Gaussian tori. Independent same-`N` confirmation and prospective Gaussian lineages are compatible with

```text
DeltaM ~ DeltaCos4 * N^(-13/8).
```

The N=185/265 prospective new-geometry block gives

```text
x=21/4 H4-like: chi2 = 3.046 / 2
zero:            chi2 = 29.409 / 2
x=17/4:          chi2 = 30.246 / 2
```

The completed norm-5 N=325/425 prospective discriminator then resolves the old odd-harmonic alias:

```text
H4:  chi2 =  0.4163 / 2
H12: chi2 = 35.1931 / 2
H8:  chi2 = 16.0120 / 2
zero:chi2 =  1.7764 / 2
```

So the frozen H4 transfer strongly outperforms the tested H12/H8 aliases. The child block alone does **not** reject zero effect, so this is a transfer/harmonic discrimination result rather than a new standalone nonzero-effect detection.

The matching-even N=185/265 result is also compatible with its frozen amplitude after the exact channel conversion

```text
DeltaS_cross = -DeltaS_either,
```

with corrected `chi2 = 0.5700315 / 2` and no refit.

## Derivative/metric problem

The prospective intrinsic-center score remains

```text
P4[S]   ~ N^-1        survives
P4[D]   ~ N^-13/8     survives
P4[D']  ~ N^-5/8      survives
P4[S']  ~ N^-5/4      fails
```

On the norm-5 full-curve derived view, q=2 and rank-2/Jordan remain inconclusive:

```text
q2 analytic: chi2 = 10.648 / 6
Jordan/log:  chi2 =  9.020 / 6
```

Two simple scalar correction stories have now failed on post-reveal P57 analyses:

- `E[K_plus-K_minus] = A N^(5/8)+B` fails its frozen N325/N425 target with joint `chi2 = 155.22 / 2`;
- using the observed canonical rank-gap width as a single scalar rescaling does not collapse the higher Krawtchouk/Hermite thermal jet.

This shifts the live mechanism question toward low-rank transfer/operator mixing rather than another scalar correction exponent.

## Highest-value next work

1. **N=145 -> 290 full curve — #50.** This is now the single highest-information new production block. Score the frozen center-slope/root correction first, then reuse the same data for held-out low-rank transfer, Krawtchouk/Hermite, rank-gap and metric-free-ratio tests.
2. **Use the completed P57 block harder.** The raw 500M N325/N425 histograms and joint moments are canonical and support additional covariance-aware analyses without new simulation.
3. **Norm-4 dyadic closure — #154.** Ready if the q=2/Jordan or low-rank transfer ambiguity survives N290.
4. **Independent controls.** Self-matching tangent (#155), pivotal/four-arm, FK/Potts, square-bond duality, and modulus/Pell work can proceed in parallel whenever they add more information per cost.

N=1105 is lower current information per CPU, not prohibited.

## Analysis infrastructure

Threshold-rank production preserves reusable sufficient statistics rather than only final decimals. The same archive now supports:

- roots, center slopes and derivative channels;
- exact finite Russo/pivotal checks;
- activation/reliability signatures;
- Krawtchouk/Hermite full-curve response coordinates;
- paired rank-gap/neutral-window observables;
- low-rank matrix/semigroup transfer analysis;
- metric-free amplitude ratios;
- primary-only prequential evidence aggregation.

`scripts/run_norm5_analysis_bundle.py` is a thin orchestration entrypoint for six-size norm-5 blocks. It calls existing frozen/typed scorers independently and does not redefine their statistics.

## Research navigation

- `docs/STATUS.md` — authoritative current claim ledger.
- `docs/RESEARCH-MAP.md` — how scientific tracks and evidence fit together.
- `notes/SYNTHESIS-20260828.md` — execution-facing synthesis.
- `docs/ROADMAP.md` — information-gain priorities.
- `analysis/research_ledger.yaml` — compact machine-readable evidence/work state.
- `analysis/artifact_registry.yaml` — lightweight artifact/navigation index.
- `results/evidence-ledger/latest.md` — primary-only predictive evidence view.

Old reports, queues, closed PRs and negative results remain provenance, not competing current-status documents.

## Execution policy

Useful analysis, exact work and pilots are allowed by default. The repository keeps only three hard constraints:

1. do not rewrite frozen predictions or committed result history;
2. do not silently score incompatible observable semantics;
3. do not add correlated views of one raw random block as independent primary evidence.

Claim strength is controlled by evidence and chronology, not process ceremony.

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

No closed form for square-site `p_c` is claimed. Published numerical estimates remain method-specific and are tracked in `data/literature_threshold_sources.json`.

## License

MIT. See `LICENSE`.
