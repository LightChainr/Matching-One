# Matching One

Matching One is a computational research project on square-lattice site percolation, its matching-lattice identity, and the finite-size/operator structure behind the threshold.

The exact structural anchor is

\[
p_c^{\mathrm{site}}(\mathbb Z^2)+p_c^{\mathrm{site}}(\mathrm{NN+NNN})=1.
\]

## Start here

- [`docs/NEXT-TARGETS.md`](docs/NEXT-TARGETS.md) — fast decision board: what can change the research picture next.
- [`docs/STATUS.md`](docs/STATUS.md) — authoritative current claim ledger.
- [`docs/RESEARCH-MAP.md`](docs/RESEARCH-MAP.md) — how the scientific tracks fit together.
- [`docs/ROADMAP.md`](docs/ROADMAP.md) — default research-attention priorities.
- [`analysis/research_ledger.yaml`](analysis/research_ledger.yaml) — machine-readable evidence, hypotheses and priorities.
- [`results/evidence-ledger/latest.md`](results/evidence-ledger/latest.md) — primary-only predictive evidence view.

Priorities are advisory. They do not lock tasks, prevent parallel work, require branch consolidation, or imply that a lower-ranked mechanism is rejected.

## Current empirical picture

The strongest numerical structure is a matching-odd orientation signal on primitive Gaussian tori. Independent primary blocks strongly disfavor global zero and remain compatible with a leading H4-like transfer

```text
DeltaM ~ DeltaCos4 * N^(-13/8).
```

The completed norm-5 N=325/425 discriminator resolves the old angular alias in favor of H4 over the tested H12/H8 alternatives:

```text
H4:   chi2 =  0.4163 / 2
H12:  chi2 = 35.1931 / 2
H8:   chi2 = 16.0120 / 2
zero: chi2 =  1.7764 / 2
```

The child block alone does **not** reject zero, so the value of this result is transfer/harmonic discrimination rather than a new standalone detection.

The N145 -> N290 full-curve block adds a different fact: the finite-size response cannot be compressed into one scalar multiplier. The corrected center/slope/root structure survives, while a resolved shape direction fails the one-multiplier description.

## The main identification problem

The project should no longer treat “spin 4” as a unique operator label.

Two exact continuum candidates now matter directly:

```text
V_(2,+/-2):  x=17/4, |spin|=4     four-leg primary candidate
Q4 epsilon:   x=21/4, |spin|=4     thermal descendant / Jordan candidate
```

The Q4 construction is exact at the representation level and has precise modular/Jordan fingerprints, but the lattice overlap is not proved. The lower `x=17/4` field is a serious competitor until its physical `Q=1` representation/multiplicity and global matching overlap are derived.

Therefore the highest-information question is increasingly a **selection-rule / representation question**, not another radial exponent fit.

## Derivative and local-state problem

Prospective/full-curve analyses support the broad pattern

```text
P4[S]   ~ N^-1        survives
P4[D]   ~ N^-13/8     survives
P4[D']  transfers cleanly
P4[S']  requires nontrivial finite-size mixing
```

Simple scalar-width and one-coordinate shell descriptions do not close the full state. Rank-2/Jordan remains viable, but ordinary diagonalizable low-rank mixing and cover/topological memory remain live alternatives.

The current local pivotal/two-cutoff results should therefore be read as evidence that the chosen scalar coordinate is insufficient, not as a blanket rejection of local or logarithmic mechanisms.

## Three questions organize the project

```text
Q1. Which continuum sector produces the global matching H4 signal?
Q2. Why does the lattice matching observable select that sector?
Q3. Can this structure explain or constrain the microscopic value of p_c itself?
```

The first two form the operator/selection program. The third is a distinct threshold-origin program: even a complete CFT operator identification does not automatically determine the nonuniversal microscopic critical probability.

## Highest-value current attention

The default order is documented in [`docs/NEXT-TARGETS.md`](docs/NEXT-TARGETS.md). In compact form:

1. **Resolve the lower spin-4 competitor and its global selection rule** — representation/projector/multiplicity, charged positive controls, and typed `Q`-parameter fingerprints.
2. **Infer the smallest predictive state before naming fields** — use existing norm-2/norm-5/N290/local/score-mode data to compare ordinary rank-2, Jordan rank-2, cover-enriched state, or no stable low-rank closure.
3. **Calibrate the machinery on independent controls** — triangular known logarithmic pair, square-bond/self-matching controls, and tunable exactly-critical anisotropy.
4. **Buy orthogonal prospective evidence** — same-N coalescence, norm-4 composition, or modulus/shape, chosen by expected model-space reduction per compute.
5. **Build the matching-observable bridge** — RG tangent, FK/Potts sector derivative, defect/intertwiner or transfer-matrix selection rule.
6. **Keep a direct threshold-origin line visible** — correlated-hyperedge/self-dual embeddings and exact finite structural work aimed at the value of `p_c` itself.

These activities can run in parallel; the list is a default attention order, not a sequence of approvals.

## Analysis infrastructure

Threshold-rank production preserves reusable sufficient statistics rather than only final decimals. The archive supports:

- roots, centers, slopes and derivative channels;
- exact Russo/pivotal checks;
- Krawtchouk/Hermite and full-curve coordinates;
- rank-gap and intrinsic-quantile observables;
- low-rank/minimal-realization and semigroup tests;
- Gaussian cover/Smith/deck arithmetic;
- metric-free amplitude ratios;
- primary-only prequential evidence aggregation.

Different analyses of the same histogram are coordinates on one raw block, not automatically independent evidence.

## Working philosophy

The repository should stay permissive about exploration while keeping claims conservative.

Useful exact work, reanalysis, pilots and side programs may proceed whenever they are cheap or open a new information axis. Priority changes are preferred over task locking or premature consolidation.

The scientific objective is to make live mechanisms predict different things as quickly as possible—not to maximize the number of explanations, and not to prematurely collapse them into one story.

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

No closed form for square-site `p_c` is claimed. Published numerical estimates remain method-specific and are tracked in `data/literature_threshold_sources.json`.

## License

MIT. See `LICENSE`.
