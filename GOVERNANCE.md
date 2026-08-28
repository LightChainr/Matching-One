# Governance

Matching One is an exploratory computational-mathematics project. Governance exists to keep evidence understandable and reversible, not to make research wait for product-style process.

## 1. Default mode: fast research integration

`main` is the shared research line.

The maintainer may merge exploratory code, notes, protocols, and result archives once they are understandable and the relevant tests/checks run. External approval is **not** required for ordinary exploratory work in a solo-maintainer repository.

Large result archives may be merged in one PR when splitting them would mostly create administrative work. Code and data may live in the same PR when their relationship is easier to audit that way.

The main distinction is not “merged versus unmerged.” It is **how strong a scientific claim the evidence supports**.

Use branches when they help isolate work, not because every small research step needs ceremony. Avoid force-rewriting published numerical history; ordinary follow-up commits and explicit corrections are preferred.

## 2. Scientific claim levels

Every important conclusion should fit one of these levels.

| Level | Meaning |
|---|---|
| C0 | hypothesis, conjecture, design, or theory candidate |
| C1 | method/control validated by exact identity, oracle, or deterministic regression |
| C2 | exploratory numerical signal; analysis may still be adaptive |
| C3 | reproduced/frozen finite-size numerical result, e.g. independent seed or prospective/held-out test |
| C4 | asymptotic/mechanistic interpretation supported by multiple discriminating tests |
| C5 | rigorous result or independently checkable certificate/proof |

A result can be on `main` at C0, C1, or C2. Merging exploratory work is not a claim upgrade.

The current project-wide summary lives in `docs/STATUS.md`; the execution-facing synthesis lives in the latest `notes/SYNTHESIS-*.md`.

## 3. What actually needs review rigor

### Exploratory notes and scripts

Usually enough:

- purpose is clear;
- script runs or compiles;
- output is labeled exploratory;
- no existing frozen result is silently overwritten.

A smoke test is preferred for numerical scripts. It does not need exhaustive coverage before the script can be useful.

### Expensive numerical runs

Preserve enough information that the run can be understood later:

- source commit or source hash;
- command/configuration;
- RNG seed/counter convention when stochastic;
- sample/batch count;
- raw sufficient statistics or the most reusable aggregates;
- a short report saying what passed and what failed.

For expensive or decisive tests, preregistering the sign/model/geometry before looking at the target is strongly preferred.

### High-risk numerical machinery

Topology, homology, RNG, threshold-rank reconstruction, covariance propagation, and exact polynomial code deserve stronger checks because a bug can contaminate many downstream experiments.

Aim for one independent/exact reference or a deterministic regression where practical. This is a priority, not a reason to freeze all downstream exploration until every edge case is formalized.

### C4/C5 or paper-facing claims

This is where independent review matters most.

Before calling an exponent/operator/universality mechanism established, seek at least one of:

- independent implementation;
- genuinely new prospective geometry/size;
- exact control model;
- independent collaborator review;
- analytic derivation/certificate.

External review is recommended for paper-facing C4/C5 claims, but it is not a merge prerequisite for ordinary research progress.

## 4. PR and branch practice

Use the smallest workflow that keeps the work legible.

Good patterns include:

- one small PR for a focused analysis change;
- one combined code+result PR for a tightly coupled experiment;
- one large archival PR for a coherent compute campaign;
- direct follow-up commits on a research branch when several experiments share the same engine.

Avoid maintaining deep stacks of PRs after their base work is already integrated. Retarget or merge them onto `main` and close superseded coordination PRs.

Merge commits are appropriate when provenance ancestry matters. Squash is fine for documentation/governance cleanup. Rebase is optional, not a policy goal.

## 5. Results and corrections

Negative results, failed models, and underpowered experiments are useful evidence. Keep them.

When an error or interpretation change is found:

1. preserve the old artifact;
2. add a correction or replacement result;
3. state what changed and why;
4. downgrade a claim in `docs/STATUS.md` if the old conclusion no longer holds.

Do not spend time engineering a correction workflow more elaborate than the scientific risk requires.

## 6. Statistical discipline without bureaucracy

For confirmatory experiments, prefer:

- frozen target/sign/model before the target run;
- held-out or prospective data when available;
- full covariance when it materially changes the conclusion;
- reporting effect sizes and uncertainties, not only p-values;
- keeping flexible/free-exponent models secondary to parameter-free tests.

For exploratory work, approximate diagnostics are acceptable if clearly labeled. A covariance or finite-sample refinement should block a strong quantitative claim only when it could realistically change that claim; it should not block qualitative exploration that is robust to the issue.

## 7. Current project operating principle

Scientific language should remain conservative; engineering integration should be fast.

In practice:

- put useful evidence on `main`;
- mark whether it is exploratory, reproduced, or prospective;
- keep the strongest three next discriminators visible in the synthesis note;
- avoid spawning a new branch/note/issue unless it creates a sharper test;
- spend process effort in proportion to the scientific consequence of being wrong.

## 8. Releases

A paper-oriented or archival release should contain a claim-ledger snapshot, source/result hashes, major limitations, and enough information to reconstruct the reported tables/figures.

A release tag does not itself upgrade a claim.
