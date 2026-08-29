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

Thus the matching curve is exactly an equal mixture of first-wrap and second-direction-completion distributions, not generically one latent threshold. The raw midpoint `C_raw=(K1+K2)/2` and raw gap `G_raw=K2-K1` retain information discarded by the scalar mixture; the normalized clock/lifetime coordinates used below are `C=C_raw/(N+1)` and `W=G_raw/(N+1)`.

The main line now contains three complementary finite checks. PR #282 identifies `K_minus/K_plus` with the first and second essential-H1 births in the exact filtration oracle. PR #283 reconstructs `P0=1-F1`, `P1=F1-F2`, `P2=F2` and the mean rank-one lifetime from an existing rank archive. PRs #284/#285/#286 extend the short-period quotient frontier through index 8: 55 HNF quotients and 654,678 filtrations have no birth, reflection, rank-sum, reconstruction or plateau-line failures. PR #286 also verifies cached versus direct semantics on all 24 axis-L2 paths and archives 5,980 subset states plus 1,967,984 rank-one plateau steps; maximum observed saturation index remains one. This is a strong finite frontier, not an unrestricted theorem for every degenerate quotient. The older production rank archive still lacks the plateau line `ell`, Smith saturation index `iota` and local birth mark, which is precisely why the next stream should retain them and reuse the exact cached oracle rather than reimplement topology.

The completed Draft reanalysis scores all ten archived sizes without new Monte Carlo production or a fitted exponent. `K1` is the larger linearized H4 root-shift point estimate at every size, while the magnitude classification is `K1`-dominant at nine sizes and shared at one. `K2` reinforces it at N=65,85,130,145,170,185,290 and has an opposite-sign point estimate at N=265,325,425. The full-curve Bernstein reuse resolves the apparent reversal: the integrated `K2` area has one positive sign at all ten sizes, while each negative root-point value belongs to a local lobe with a delete-one-stable nearby zero branch at approximately `.59462`, `.59635` or `.59694`. The node-minus-`p_bar` offsets are smaller than their position errors, so no significant ordering is claimed. See [`results/two-activation-h4/latest.md`](results/two-activation-h4/latest.md) and [`results/activation-curve-nodes/latest.md`](results/activation-curve-nodes/latest.md); the isolated 16-worker ARM64 replay of the root split is recorded in [`results/two-activation-h4/server-run.json`](results/two-activation-h4/server-run.json).

An additional Draft reuse of the branch-only P205 quotient prism closes the cheapest character question. A six-coordinate GLS selects `H4/H4` for the two activations (`chi2=2.585155/4`, `p=.629455`): both components reinforce, with fitted amplitudes `A1=.508581` and `A2=.295266`; `K2` supplies 36.732% of the signed fitted H4 amplitude. This retrospective small-N result identifies neither a continuum field nor an asymptotic law, but it shows that the ordinary P205 H4 signal is not purely a first-birth effect. See [`results/two-activation-prism/latest.md`](results/two-activation-prism/latest.md).

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

The unmerged double-projector staircase sharpens the ordinary rank-plane question: the vacuum/KdV spin-4 direction lies on the Alexander-even line and has `A_top=0`, while the regular unlabeled `[2]` four-leg carrier also has zero ordinary overlap. Thermal `Q4 epsilon`/Jordan is the first listed candidate allowed by both selectors. “Allowed” is still not “coupled”: the missing discriminator is a nonzero lattice-to-module coupling together with its modulus/Jordan fingerprint. A new complex-C3 acquisition should therefore name a typed non-`A_top` module or charged observer, rather than repeat the already-resolved ordinary H4/H8/H12 character competition.

Main also contains an exact Q=1 velocity fingerprint that makes the future Q-tangent test sharper. The four-leg `V_(2,+/-2)` row has `dx/dQ=-5 sqrt(3)/(16 pi)`, while thermal `Q4 epsilon` has `-9 sqrt(3)/(16 pi)`; their velocity gap is `sqrt(3)/(4 pi)`. The generic-loop spin-8 and spin-12 rows in that oracle are controls, not assignments of the experiment-design H8/H12 aliases. The analytic velocities therefore do not need to be rederived; the missing work is the lattice measure + projector + explicit-field interface and its overlap.

### The finite-size state is compact but not scalar

The N145->290 full-curve result rejects a one-multiplier description through a resolved shape direction. Norm-4 subsequently rejects the frozen analytic q2 scalar and common thermal-jet generator. The post-reveal “Jordan plus one conjugation-even mode” recurrence then passes its frozen N520/N680 generation-four pilot at `lambda=1/2` (`1.314/2`, `p=.518` for scalar U; `9.298/10`, `p=.504` for the jet), removing the visible residual tension. Nearly identical `lambda=0` and `lambda=1` scores mean that pure Jordan, the analytic even-mode choice and persistent curvature remain unresolved; the pilot does not establish a nonzero extra mode or a unique transfer matrix. This raises the value of a geometry/operator coordinate chosen to maximize the frozen models' Mahalanobis separation rather than another free lambda fit.

In the annulus line, PR #247 supplies a useful split inside one correlated raw block: the ordinary/local plus-shell control is compatible with one common per-log amplitude (`p=.3564`, amplitude about 3.3 source standard errors from zero), the matching-odd minus-shell scalar law fails (`p=.00382`), and matched-cutoff equality is marginally tense (`p=.04544`). A two-state radial recurrence determined from N325/N425 nevertheless passes a held-out N365 third geometry with joint `p=.965`. The correct update is therefore not “radial structure failed,” but “the ordinary control works, one odd scalar coordinate fails, and a two-state recurrence survives.”

### Geometry and acquisition carry different kinds of memory

The exact CRT join says an unmarked final quotient does not remember the order of its factors. An exact path-filtration witness nevertheless shows that intermediate homology ranks can record when rank two first appears. Chronological/morphism memory therefore belongs to an intermediate filtration, dynamic lineage or marked acquisition—not to the final endpoint alone.

The triangular-site energy/log-pair control has also moved beyond implementation readiness: an open-PR cross-cutoff production score passes its frozen joint model and resolves a nonzero top-partner cutoff shear. Its `kappa_proxy` is gauge-dependent; a top field with two normalized macroscopic spin insertions is the cheapest proposed gauge-invariant continuation.

The branch-only rank-birth algebra now supplies the next common-stream interface. The exact birth gates `I01` and `I12` give an even rank-increment row `S=I01+I12=M'` and an odd lifetime row `D=I12-I01=-partial_p P1`; every nonzero `D` insertion carries a canonical plateau line `ell`. Complement symmetry also separates the normalized midpoint `C` as an odd clock translation from the gap `W` as an even rank-one lifetime. Tiny exact/Gaussian controls give a nonzero finite `A_top` coupling to the line-resolved `D`-H4 source, removing a symmetry-zero obstruction without yet establishing survival scaling or a Q4 identity.

## Five coordinates for every research node

The scientific atlas uses five independent coordinates:

| coordinate | examples | question it prevents us from conflating |
|---|---|---|
| **state** | finite topological state, continuum module, low-dimensional dynamics, microscopic threshold structure | What object is evolving or being identified? |
| **source** | Bernoulli thermal, rank-birth geometric tilt, Potts-Q measure, projector and explicit-field derivatives, relative cluster fugacity, deck charge, boundary source | What perturbation creates the response? |
| **observer** | global unmarked/rank projector, topology-typed, rank-birth line/landing, matching-even/odd, local pivotal, deck-charged, boundary/cross-microscopic | Which sector is actually being measured? |
| **geometry** | Gaussian cover, annulus, modulus/Hecke child, boundary cross-ratio, triangular/square controls | What change supplies the discriminator? |
| **acquisition** | existing-data reuse, exact construction, small pilot, new production | What new cost and independence does the result have? |

The same exponent or harmonic in two different observer sectors is not automatically the same operator. Several views derived from one histogram are coordinates on one random block, not independent evidence.

## Default attention

The current default order is:

1. bridge the established global ordinary H4 response to rank-birth mechanism on one common stream retaining `q`, `I01/I12`, `ell`, `iota`, line-H4, landing-H4 and the required cross-moments;
2. measure a typed non-`A_top` module or charged complex C3 response across the `rho` geometries, since the ordinary P205 H4 character is already resolved and `A_top` kills the vacuum/KdV line;
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

No closed form for square-site `p_c` is claimed. The exact threshold-origin chain now fixes block-event semantics, an O(`s^2`) evaluator and a familywise plan with 400 trials and acceptance cutoff 373, but contains no fresh IID result. Its next empirical step is acquisition under that frozen plan. Separately, the exact microscopic H4 stencil gate shows that a nonnegative axis-only family cannot cancel its fourth-moment proxy; an admissible search needs a negative-phase oblique orbit (axis/integer-diagonal cancellation is `4:1` in the declared normalization). Neither statement is a theorem about the renormalized H4 coupling. Published numerical estimates remain method-specific and are tracked in `data/literature_threshold_sources.json`.

## License

MIT. See `LICENSE`.
