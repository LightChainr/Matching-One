# Matching One

Matching One studies square-lattice site percolation through its matching-lattice identity, finite topology, operator sectors and the microscopic origin of the threshold. The exact structural anchor is

\[
p_c^{\mathrm{site}}(\mathbb Z^2)+p_c^{\mathrm{site}}(\mathrm{NN+NNN})=1.
\]

The repository is organized to expose the next mechanism-changing observation, not to turn priorities into permissions. Exact work, reanalysis, pilots and independent theory lines may proceed in parallel. A lower priority is not a rejection; no task is locked by the overview.

## Start here

- [`docs/NEXT-TARGETS.md`](docs/NEXT-TARGETS.md) — seven highest-information moves and their exact decision outputs.
- [`docs/STATUS.md`](docs/STATUS.md) — authoritative claim ledger, including the boundary between `main` and unmerged frontier work.
- [`docs/RESEARCH-MAP.md`](docs/RESEARCH-MAP.md) — the state/source/observer/geometry/acquisition atlas.
- [`docs/ROADMAP.md`](docs/ROADMAP.md) — default attention order and parallel research portfolio.
- [`analysis/research_ledger.yaml`](analysis/research_ledger.yaml) — machine-readable nodes, sectors, sources, experiments and dependency groups.
- [`results/evidence-ledger/latest.md`](results/evidence-ledger/latest.md) — primary-only predictive evidence view.

## The current scientific picture

### One exact finite object, two activations

On honest periodic square-cell tori, the integrated digital Alexander theorem gives, configuration by configuration,

```text
r_black + r_white = 2,
q = r_black - 1 = 1 - r_white.
```

Along a uniformly ordered site permutation, the ambient black homology rank is nondecreasing and activates at most twice. With `K1=K_minus` and `K2=K_plus`,

```text
q_n = -1 + 1{n>=K1} + 1{n>=K2},
M_N(p) = -1 + E[H_K1(p)] + E[H_K2(p)].
```

Thus the matching curve is exactly an equal mixture of first-wrap and second-direction-completion distributions, not generically one latent threshold. The paired midpoint `C=(K1+K2)/2` and gap `G=K2-K1` retain information discarded by the scalar mixture.

The completed Draft reanalysis scores all ten archived sizes without new Monte Carlo production or a fitted exponent. `K1` is the larger linearized H4 root-shift component at every size. `K2` reinforces it at N=65,85,130,145,170,185,290 and has an opposite-sign point estimate at N=265,325,425. Those three cancelling `K2` terms are individually below `|z|=2`, so the sign change is a geometry/sector clue rather than a confirmed transition. The full estimates, aligned jackknife covariance and dependency groups are in [`results/two-activation-h4/latest.md`](results/two-activation-h4/latest.md).

Each activation contributes a strictly positive beta density for `0<p<1`, so `M_N'(p)>0`; together with the endpoint signs this gives one unique, simple physical root. This is a finite-volume theorem, not a claim that the root has a closed form.

### The global ordinary and charged sectors disagree usefully

Independent square-site primary blocks strongly disfavor global zero and remain compatible with an H4-like ordinary, unmarked response. The completed norm-5 ordinary transfer favors H4 over the tested H12/H8 aliases, although the child block alone does not reject zero.

A different, deck-charged N325 observable instead gives a frozen likelihood ordering of approximately

```text
H8 / H12 / H4 = 71 / 21 / 8.
```

This is not a discovered spin-8 field: H12 and H4 are not both rejected at a strict 5% gate. It is evidence that the charged and ordinary global observables should not be compressed into one scalar angular amplitude.

### “Spin 4” is not a field identity

The main continuum candidates remain

```text
V_(2,+/-2): x=17/4, |spin|=4, four-leg [2] sector
Q4 epsilon:  x=21/4, |spin|=4, thermal singlet descendant/Jordan candidate
```

An unmerged exact global-endpoint calculation finds an ordinary linear selection zero for the four-leg `[2]` carrier under the explicit hypothesis that the matching observable is a regular unlabeled zero-leg generic-Q endpoint. A separate unmerged width-four Q=4 transfer calculation reproduces the singlet-to-`[2,2]` matrix-element zero while showing that charged/twisted traces can still see that block. This narrows the global selection question; it does not prove absence from every trace, singular Q derivative, defect sector or marked observable.

### The finite-size state is compact but not scalar

The N145->290 full-curve result rejects a one-multiplier description through a resolved shape direction. Norm-4 subsequently rejects the frozen analytic q2 scalar and common thermal-jet generator, while the corresponding Jordan laws survive at 5% with visible tension. Post-reveal work suggests “Jordan plus one conjugation-even mode” as an economical next transfer model, not as an identified unique matrix.

In the annulus line, the simplest matching-odd scalar shell law fails, yet a two-state radial recurrence determined from N325/N425 passes a held-out N365 third geometry with joint `p=.965`. The correct update is therefore not “radial structure failed,” but “one scalar coordinate failed while a two-state recurrence survived.”

### Geometry and acquisition carry different kinds of memory

The exact CRT join says an unmarked final quotient does not remember the order of its factors. An exact path-filtration witness nevertheless shows that intermediate homology ranks can record when rank two first appears. Chronological/morphism memory therefore belongs to an intermediate filtration, dynamic lineage or marked acquisition—not to the final endpoint alone.

The triangular-site energy/log-pair control has also moved beyond implementation readiness: an open-PR cross-cutoff production score passes its frozen joint model and resolves a nonzero top-partner cutoff shear. Its `kappa_proxy` is gauge-dependent; a top field with two normalized macroscopic spin insertions is the cheapest proposed gauge-invariant continuation.

## Five coordinates for every research node

The scientific atlas uses five independent coordinates:

| coordinate | examples | question it prevents us from conflating |
|---|---|---|
| **state** | finite topological state, continuum module, low-dimensional dynamics, microscopic threshold structure | What object is evolving or being identified? |
| **source** | Bernoulli thermal, Potts-Q measure, projector derivative, relative cluster fugacity, deck charge, boundary source | What perturbation creates the response? |
| **observer** | global unmarked, topology-typed, local pivotal, matching-odd, deck-charged, boundary/cross-microscopic | Which sector is actually being measured? |
| **geometry** | Gaussian cover, annulus, modulus/Hecke child, boundary cross-ratio, triangular/square controls | What change supplies the discriminator? |
| **acquisition** | existing-data reuse, exact construction, small pilot, new production | What new cost and independence does the result have? |

The same exponent or harmonic in two different observer sectors is not automatically the same operator. Several views derived from one histogram are coordinates on one random block, not independent evidence.

## Default attention

The current default order is:

1. explain why the completed `K1/K2` decomposition is first-activation dominated while the `K2` point estimate changes from reinforcement to cancellation on three geometries;
2. measure the complex C3 response across the `rho` geometries to distinguish scalar-Q4 zero, weight-4 KdV and weight-8 conjugate characters;
3. form a common-source/readout Gaussian/annulus context rectangle to test one generator versus a morphism-enriched state;
4. turn the triangular log-pair control into a gauge-invariant cross-microscopic H4/E6/top-field statistic;
5. calibrate the full Potts-Q derivative—measure, confluent projector and explicit field—on a boundary or tiny VJS control;
6. reuse the N325 charged engine for the two primitive Z5 cubic fusion channels;
7. move the relative-source theory from its rigid three-sector scalar algebra to a connectivity/defect radical.

These are parallel attention lanes, not approval gates.

## Claim and provenance discipline

The overview uses four lifecycle labels:

- `main_integrated`: present on the current shared `main` line;
- `open_pr`: inspectable in an open PR but not yet integrated;
- `branch_only`: inspectable remote branch, neither promoted to `main` nor treated as canonical;
- `hypothesis`: a proposed mechanism or future discriminator.

Frozen chronology is preserved. Claim-bearing comparisons require identical observable semantics or an exact registered map, and correlated views of one random block are not counted as independent primary evidence. Those are the only evidence constraints; they do not serialize research.

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

No closed form for square-site `p_c` is claimed. Published numerical estimates remain method-specific and are tracked in `data/literature_threshold_sources.json`.

## License

MIT. See `LICENSE`.
