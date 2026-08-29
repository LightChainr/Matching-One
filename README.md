# Matching One

Matching One is a computational research project on square-lattice site percolation, its matching-lattice identity, and the finite-size structure behind the threshold.

The exact structural anchor is

\[
p_c^{\mathrm{site}}(\mathbb Z^2)+p_c^{\mathrm{site}}(\mathrm{NN+NNN})=1.
\]

## What we currently observe

The strongest numerical structure is an orientation-dependent matching signal on primitive Gaussian tori.

- Independent 100M same-`N` confirmation at `N=65,85,130,145,170` reproduces the prescribed odd-square-harmonic sign and is compatible with

  ```text
  DeltaM ~ DeltaCos4 * N^(-13/8).
  ```

- Three prospective `1+i` Gaussian lineages are compatible with the no-fit raw-contrast transformation

  ```text
  DeltaM(2N)/DeltaM(N) = -2^(-13/8).
  ```

- The genuinely new N=185/265 full-curve target block, 500M paired permutations per size, gives

  ```text
  DeltaM, x=21/4 H4-like: chi2 = 3.046 / 2
  DeltaM, zero:            chi2 = 29.409 / 2
  DeltaM, x=17/4:          chi2 = 30.246 / 2
  ```

  so the matching-odd `x=21/4` H4-like radial law survives the new-geometry test and strongly outperforms these two frozen alternatives.

### Important channel-map erratum

The original #108 report described the matching-even N=185/265 result as a prospective sign reversal. That interpretation was a protocol error, not a new physical effect.

The frozen matching-even source amplitude was fitted from P31 `either/even`, while the threshold-rank target reconstructs rank-2 `cross/even`. Complementary torus topology gives

```text
DeltaS_cross = -DeltaS_either.
```

After applying only this exact source-to-target channel map, with no refit,

```text
corrected DeltaS chi2 = 0.5700315 / 2
marginal residual z   = +0.67, -0.12
```

The matching-even prospective result is therefore compatible with the frozen amplitude after channel conversion. The raw #108 files remain preserved.

The durable lesson is narrow: claim-bearing scores must state their observable semantics and apply an exact registered map when source and target conventions differ.

### A real prospective failure remains

The normalized derivative law

```text
P4[S'] ~ N^-5/4
```

fails prospectively on N=185/265, while two predeclared corrections remain viable:

```text
pure N^-5/4:  chi2 = 52.716 / 2
rank-2/log:    chi2 =  1.204 / 2
analytic 1/N:  chi2 =  0.862 / 2
zero:          chi2 = 1278.555 / 2
```

The correction mechanism is not yet uniquely identified.

## Working interpretation

A compact empirical law for the central matching-odd sector remains

\[
\Delta M_N \approx A\,\Delta\cos(4\theta)\,N^{-13/8}.
\]

It has independent-seed, held-out, Gaussian-semigroup, root-closure, and prospective-new-geometry support. It does **not** yet prove that H4 is unique rather than H12/H20, nor uniquely identify an `x=21/4` LCFT operator.

The matching-even sector is compatible with its frozen `N^-1` amplitude once cross/either semantics are aligned. Full-curve data also resolve a small finite-size correction to the bare center-slope multiplier `2^(3/8)`.

The finite threshold-rank archive is now more than a source of roots and slopes: it supports reliability signatures, pivotal identities, Krawtchouk/Hermite response modes, paired rank-gap width, and low-rank transfer analysis without rerunning the simulation.

## Next experiments and analysis

The priority order is short, but it is a scheduling order rather than a permission system.

1. **Norm-5 N=325/425 — #57.** Highest-information new block. Keep the frozen H4/H12/H8/zero primary score, then reuse the same full threshold-rank data for q=2/Jordan, Krawtchouk/Hermite thermal-jet, paired rank-gap width, root/slope, derivative and multi-u analyses.
2. **Third full-curve lineage — #50.** Score `145 -> 290`, including the frozen finite-size slope correction and induced root prediction, then use the same data as a held-out score-mode/low-rank transfer test.
3. **Use existing full curves harder.** The prequential evidence ledger (#95) is complete. Active source-data routes include intrinsic quantile centers (#101), multi-u response (#119), joint operator mixing (#125), full standardized threshold profiles (#122), and pivotal/four-arm continuation (#100/#121).
4. **Choose later geometry by information gain — #102.** Maximin Gaussian design tooling is canonical; update cost/variance inputs as new runs arrive rather than defaulting to larger N.

Norm-4, self-matching tangent, Pell/modulus work, N=1105, exact algebra and other controls are valid parallel directions. Their current ranking reflects expected information per compute, not a blanket gate.

### One-command norm-5 analysis

When a six-size norm-5/full-curve data block is available, `scripts/run_norm5_analysis_bundle.py` provides a thin orchestration entrypoint. It infers exact shared counter groups from metadata and runs the existing scorers independently. A frozen scorer that is inapplicable to an exploratory input is recorded as a failed subanalysis without blocking the other views.

```text
--run 65:HIST:MOMENTS:METADATA
--run 85:HIST:MOMENTS:METADATA
--run 130:HIST:MOMENTS:METADATA
--run 170:HIST:MOMENTS:METADATA
--run 325:HIST:MOMENTS:METADATA
--run 425:HIST:MOMENTS:METADATA
```

This keeps one expensive data block reusable rather than turning every new observable into a new production campaign.

## Research navigation

- `docs/STATUS.md` — authoritative current claim ledger.
- `docs/RESEARCH-MAP.md` — how scientific tracks and evidence fit together.
- `notes/SYNTHESIS-20260828.md` — execution-facing scientific synthesis.
- `docs/ROADMAP.md` — information-gain priorities.
- `analysis/research_ledger.yaml` — machine-readable questions, evidence and work state.
- `analysis/artifact_registry.yaml` — lightweight artifact/navigation index; it does not block research integration.
- `results/evidence-ledger/latest.md` — primary-only predictive evidence view.

Old reports, wave notes, queues, closed PRs and negative results remain provenance, not competing current-status documents.

## Repository

```text
constants/      reference values and exact relations
data/           literature datasets and provenance
docs/           canonical status, research map, and roadmap
analysis/       evidence/work indexes and analysis manifests
notes/          synthesis, theory, derivations, negative results
scripts/        production analysis, scorers, exact checks, orchestration
experiments/    frozen and historical protocols
predictions/    preregistered/frozen predictions
results/        immutable raw and derived research archives
tests/          smoke, regression, exact-contract tests
src/            production C++ engines
```

This is an exploratory mathematics/computational-physics repository. Useful research assets should enter `main` quickly; claim strength is controlled by evidence and chronology, not by branch location or process ceremony.

The repository-wide hard constraints are intentionally minimal:

1. do not rewrite frozen predictions or committed result history;
2. do not silently score incompatible observable semantics;
3. do not add correlated views of one raw random block as independent primary evidence.

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

No closed form for square-site `p_c` is claimed. Published numerical estimates remain method-specific and are tracked in `data/literature_threshold_sources.json`.

## License

MIT. See `LICENSE`.
