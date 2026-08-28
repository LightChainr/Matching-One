# Governance

This document defines how Matching One turns exploratory work into canonical code, data, and scientific claims.

## 1. Canonical repository state

`main` is the reviewed integration line. A branch, issue, notebook, server directory, or pull-request description is evidence and working history, but it is not canonical until merged to `main`.

Direct commits and force-pushes to `main` are prohibited by policy. Changes should arrive through pull requests and pass the repository checks. Repository rules should enforce this policy when the hosting integration permits it.

Research branches should be short-lived and named by purpose:

- `research/<topic>` for theory, experiments, and analysis;
- `fix/<topic>` for defects and corrections;
- `governance/<topic>` for repository policy and maintenance;
- `archive/<topic>` only for immutable historical imports.

A stacked research branch must identify its parent PR. Once the stack is integrated, obsolete PRs should be closed as superseded rather than left indefinitely ambiguous.

## 2. Roles

### Maintainer

The repository owner is the final maintainer for operational decisions. The maintainer may merge, close, revert, or archive work, but scientific claim upgrades remain subject to the evidence requirements below.

### Contributor

A contributor proposes code, data, theory, experiments, or documentation and is responsible for the provenance and tests of that contribution.

### Independent reviewer

For high-impact numerical or theoretical claims, an independent reviewer should verify the reasoning, implementation, or reconstruction without relying solely on the author's execution path. Independence may mean a separate implementation, derivation, seed/counter range, machine, or data transcription.

## 3. Scientific claim levels

Every result report and claim-bearing PR should state one of these levels.

| Level | Name | Minimum meaning |
|---|---|---|
| C0 | Proposal | Hypothesis, design, or conjecture; no empirical support required |
| C1 | Validated method/control | Exact identity, regression vector, oracle agreement, or implementation contract has passed declared controls |
| C2 | Exploratory signal | Effect observed in a non-final or discovery analysis; model and sample choices may still be adaptive |
| C3 | Confirmed finite-size result | Independent seed or implementation plus a frozen or held-out test supports the declared finite-size statement |
| C4 | Asymptotic/mechanistic interpretation | Multiple sizes/geometries and discriminating alternatives support an asymptotic exponent, universality class, or operator mechanism |
| C5 | Rigorous result | A proof or independently checkable certificate establishes the statement under explicit assumptions |

A higher level must not be inferred from sample size alone. In particular, C3 evidence does not automatically establish a C4 exponent or operator identification.

The current project-wide status is recorded in `docs/STATUS.md`. A claim-level upgrade requires a dedicated PR that updates the ledger and links the supporting artifacts.

## 4. Change classes and review gates

### Documentation-only changes

Required:

- no change to numerical artifacts or scientific meaning unless declared;
- links and terminology checked;
- CI passes.

### Code changes

Required:

- regression tests for the changed contract;
- deterministic behavior where promised;
- no silent change to frozen conventions, geometry order, units, or RNG domains;
- performance claims supported by end-to-end measurements;
- C++ production paths checked against a Python or exact oracle where feasible.

Topology, homology, RNG, threshold-rank, and covariance code are high-risk. A new implementation of one of these components should be checked against a structurally independent reference, not only against itself.

### Experiment protocols

A production experiment must declare before evaluation:

- hypothesis and competing models;
- geometry/size set and orientation order;
- training, pilot, evaluation, and held-out partitions;
- RNG domain, seed/counter policy, and batch structure;
- primary statistic and covariance model;
- power or sensitivity target;
- stopping rule and acceptance/falsification criteria;
- required raw sufficient statistics and metadata.

A pilot may select sample size or frozen variance-reduction weights. It must not select them from whether the point estimate looks favorable.

### Result archives

Required:

- immutable raw sufficient statistics or exact source data;
- source commit and source-file hash;
- executable hash for production binaries;
- environment, compiler/interpreter, dependency versions, and commands;
- seed/counter ranges and batching;
- checksums;
- analysis code and a concise `REPORT.md`;
- positive, negative, and failed-model results.

Previously published result files must not be overwritten. Corrections use a new directory or append-only correction record that links the superseded artifact.

### Scientific conclusions

A conclusion PR must separate:

1. direct observations;
2. model-dependent deductions;
3. asymptotic or theoretical interpretation;
4. known failure modes and unresolved alternatives.

The PR must not use words such as “proved,” “exact,” “universal,” or “confirmed” beyond the applicable claim level.

## 5. Pull-request policy

Prefer one logical change per PR. Code/protocol changes and bulk result imports should normally be separate so that executable logic can be reviewed without thousands of generated lines.

A large result PR must provide a compact manifest and must identify which files are source, generated, raw, derived, and narrative. Generated plots should be reproducible from committed data.

Merging uses the method appropriate to the history:

- merge commits for an intentional stacked branch whose ancestry must be preserved;
- squash merge for a self-contained governance or maintenance PR;
- rebase only when commit identity is not part of provenance.

## 6. Decisions, corrections, and reversals

Scientific disagreement is resolved by adding discriminating tests, not by deleting inconvenient outputs. Negative results and failed gates remain first-class artifacts.

When an error is found:

1. open or update an issue describing the affected scope;
2. preserve the original artifact;
3. add a corrected artifact and compatibility analysis;
4. downgrade the claim in `docs/STATUS.md` if required;
5. identify whether the cause was code, data, covariance, provenance, interpretation, or documentation.

A revert is preferred to an opaque history rewrite.

## 7. Releases

A research release should contain:

- a claim ledger snapshot;
- a data and result manifest with hashes;
- tested source and dependency information;
- known limitations;
- a citation record.

A release tag does not itself raise a scientific claim level.

## 8. Governance changes

Changes to this document require a pull request. Material changes to claim levels, reproducibility requirements, or merge policy should explain their effect on existing work.
