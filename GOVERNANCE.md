# Governance

Matching One is an exploratory computational-mathematics project. Governance exists to increase discovery speed while keeping the evidence chain legible. The default is **do the useful work, integrate it quickly, and control claim strength afterward**.

## Current phase — user-directed contraction, 2026-08-31

The user has now asked this workstream to stop automatically adding archive observables, prefix features and generic certificates. The two prespecified P0 experiments, prospective P154 transmission and P334 independent intervention, have completed and triggered their fixed downgrade rules. Their broader issues remain P1; no new production is automatically queued. A new P0 analysis must state, before target generation, which outcome ends a specified mechanism line; freeze its predictions, score and fixed budget first. Archive analyses remain C2 even when individual calculations have valid covariance or cross-validation, because hypothesis generation has repeatedly used that archive. Failed candidates are downgraded before any post-result replacement model is proposed.

This dated instruction governs the present workstream and supersedes the automatic-expansion and merge defaults below. It does not change frozen history or erase exploratory assets. The user also requests Draft delivery without merging. Current scientific status and execution order remain in `docs/STATUS.md` and `docs/NEXT-TARGETS.md`; the older general policy below is retained for provenance.

## 1. Default mode: run and integrate

`main` is the shared research line, not a publication-only branch.

Useful scripts, exact calculations, source-data reanalyses, pilots, theory notes, frozen predictions, negative results, and result archives should normally enter `main` as soon as they are understandable and minimally checked. External approval is not required for ordinary exploratory work in a solo-maintainer repository.

A registry or documentation conflict must not block a scientifically useful analysis asset. Integrate the asset first and repair navigation in a follow-up commit when that is faster.

Branches and PRs are coordination tools, not permission gates. Close or bypass duplicate entrypoints aggressively once the useful content is canonical.

## 2. The three hard constraints

Only three repository-wide constraints should routinely block a claim-bearing score or evidence aggregation.

### A. Do not rewrite frozen chronology

A prediction or scoring contract frozen before target reveal stays immutable. Committed raw results and historical reports are preserved. Corrections are additive: keep the old artifact and add the corrected interpretation or replacement result.

This does **not** prevent new models, new analyses, new observables, or post-reveal discovery. It only prevents presenting a later change as if it had been frozen earlier.

### B. Do not silently compare different observable semantics

For a claim-bearing scorer, channel, primal/matching or even/odd combination, probability coordinate, orientation order, normalization, and quantity must either match or use a named exact map. Unsupported mappings fail closed.

Exploratory work may inspect alternative conventions freely, provided it is labeled as exploratory and does not silently enter a frozen score.

### C. Do not count one random block as several independent evidence blocks

Multiple roots, slopes, derivatives, score modes, quantiles, or other views derived from the same histograms may all be useful. They simply must not be added as independent primary evidence unless their joint covariance/evidence construction justifies it.

Everything else in this document is guidance, priority, or claim-labeling policy rather than a permission gate.

## 3. Scientific claim levels

Every important conclusion should fit one of these levels.

| Level | Meaning |
|---|---|
| C0 | hypothesis, conjecture, design, or theory candidate |
| C1 | method/control validated by exact identity, oracle, or deterministic regression |
| C2 | exploratory numerical signal; analysis may be adaptive |
| C3 | reproduced/frozen finite-size numerical result, e.g. independent seed or prospective/held-out test |
| C4 | asymptotic/mechanistic interpretation supported by multiple discriminating tests |
| C5 | rigorous result or independently checkable certificate/proof |

A result can be on `main` at any level. Merging is not a claim upgrade and lack of preregistration is not a reason to discard useful C2 evidence.

The current project-wide summary lives in `docs/STATUS.md`; the current analysis sequence lives in `docs/NEXT-TARGETS.md`. Dated `notes/SYNTHESIS-*.md` files and `docs/ROADMAP.md` preserve their historical context.

## 4. Research execution policy

### Existing-data analysis

Default: **run it**.

If the required sufficient statistics already exist, analysis should not wait for a separate issue, approval, or roadmap promotion. Add exact/semantic checks where cheap, record whether the result is retrospective or prospective, and integrate useful outputs.

### Exact calculations and controls

Default: **run them**.

Bounded enumeration, algebraic identities, tiny-system oracles, symbolic calculations, and synthetic controls are low-friction ways to kill weak ideas early. They do not need to wait behind the primary compute queue.

### Pilots and method experiments

Default: **run small pilots**.

A pilot may be used to estimate variance, runtime, numerical stability, or information gain. Pilot outcomes should not be relabeled as prospective target evidence, but they may immediately change engineering choices and future priorities.

### New production

There is no general production gate. Run new production when the expected information gain justifies the cost and the observable is sufficiently defined to interpret the result.

For an expensive confirmatory question, freezing the target/model/sign/score before reading the target is strongly preferred because it can support C3 evidence. If that is not done, the run remains useful exploratory C2 evidence rather than being blocked or discarded.

GPU, large CPU, Pell, N=1105, norm-4, norm-5, modulus scans, or other campaigns are therefore **priority decisions, not permission classes**. The roadmap may say “later”, “low leverage”, or “needs a better observable”; it should not imply that useful exploratory computation is forbidden.

### Sequential stopping

Predeclared e-process/confidence-sequence stopping is preferred when a run is intended as confirmatory evidence. Exploratory monitoring is allowed; it simply cannot inherit the same optional-stopping guarantee after the fact.

## 5. High-risk numerical machinery

Topology, homology, RNG, threshold-rank reconstruction, covariance propagation, and exact polynomial machinery can contaminate many downstream analyses if wrong. These deserve exact or independent regressions where practical.

The check is proportional to reuse and consequence. It should not freeze unrelated exploratory work while a stronger oracle is being developed.

## 6. PR and branch practice

Use the smallest workflow that keeps work legible.

Good patterns include direct integration of focused research assets, one small PR for a reusable method, one combined code+result PR for a tightly coupled experiment, or one archival PR for a coherent compute campaign.

Prefer canonicalization over stacked coordination debt:

- if a PR contains a useful standalone asset, integrate it even if its registry edit conflicts;
- if a later PR supersedes an earlier theory/protocol entrypoint, keep the history but close the duplicate;
- do not maintain deep stacks merely to preserve process order;
- merge commits are useful when provenance ancestry matters; squash/rebase are optional.

## 7. Results and corrections

Negative, null, failed, underpowered, and contradictory results are first-class research assets. Keep them.

When an error or interpretation change is found:

1. preserve the old artifact;
2. add a correction or replacement result;
3. state what changed and why;
4. update `docs/STATUS.md` only if the current claim boundary changes.

Do not build a more elaborate correction workflow than the scientific risk requires.

## 8. Statistical discipline without bureaucracy

For strong quantitative claims, prefer full covariance, effect sizes, held-out/prospective data, and parameter-free tests. For discovery, flexible fits and approximate diagnostics are allowed when labeled accordingly.

A statistical refinement should block a strong claim only when it could realistically change that claim. It should not block qualitative exploration, alternative coordinates, mechanism discovery, or generation of the next sharper test.

## 9. Roadmap semantics

Roadmap labels are scheduling hints:

- **active** — high expected information now;
- **ready** — can be run/analyzed whenever resources are available;
- **later** — lower information per cost at present;
- **dependency** — another result would make interpretation sharper, but exploratory work may proceed;
- **historical** — retained for provenance.

Avoid using `gated`, `blocked`, or `do not start` except when one of the three hard constraints would be violated.

## 10. Operating principle

Scientific language stays conservative; research execution stays permissive.

In practice:

- integrate useful analysis quickly;
- use existing sufficient statistics harder before assuming new data are required;
- run cheap exact controls and pilots early;
- choose expensive work by information gain rather than ceremony;
- keep chronology and observable semantics explicit;
- aggregate correlated evidence once, not repeatedly;
- let failed tests redirect the program immediately.

## 11. Releases

A paper-oriented or archival release should contain a claim-ledger snapshot, source/result hashes, major limitations, and enough information to reconstruct the reported tables/figures. A release tag does not itself upgrade a claim.
