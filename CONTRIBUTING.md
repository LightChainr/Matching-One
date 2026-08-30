# Contributing

Matching One accepts contributions to numerical methods, exact checks, statistical analysis, theory, data provenance, documentation, and research infrastructure.

Read `GOVERNANCE.md` and `REPRODUCIBILITY.md` before proposing production computation or a claim-bearing result.

## Start with an issue

Open an issue before work that changes a scientific protocol, consumes substantial compute, imports a dataset, or introduces a new interpretation. The issue should identify:

- the question and why it matters;
- the current evidence and literature boundary;
- the proposed discriminator or deliverable;
- dependencies and resource requirements;
- acceptance and failure criteria.

Small bug fixes and documentation corrections may proceed directly to a pull request.

## Branches

Use one focused branch:

```text
research/<topic>
fix/<topic>
governance/<topic>
```

Do not build a long hidden stack. If a branch depends on another open PR, state the dependency in both PR descriptions.

## Development environment

The baseline code supports Python 3.9 and later. Research branches may add dependencies in `requirements.txt`; production work should also record exact versions in the result metadata or an environment file.

Typical checks are:

```bash
python3 -m compileall -q scripts tests
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

For C++ code, use C++17 or later only when required, record the compiler and flags, and retain deterministic self-tests. GitHub Actions compiles source files and runs declared self-tests where available.

## Code requirements

### Python

- Prefer standard-library code for reference/oracle implementations.
- Use type hints on public functions and data contracts.
- Validate inputs and fail explicitly on missing sizes, singular covariance, invalid rank conventions, or incomplete metadata.
- Keep parsing, estimation, model selection, and reporting separable.
- Avoid binary floating-point when exact decimal transcription or arbitrary precision is part of the scientific contract.

### C++

- Make RNG and reductions independent of thread scheduling where promised.
- Use fixed-width integer types for counters and histograms; document overflow bounds.
- Keep a slow reference path or exact regression vectors.
- Do not enable `fast-math` when it can alter Bernoulli decisions, rank ordering, or reproducibility.
- Record source and executable hashes for production runs.

### Tests

A test should protect a scientific or software contract, not merely repeat the implementation. Important examples include:

- exact polynomial or enumeration vectors;
- independent union-find/BFS topology agreement;
- basis and geometry invariance;
- RNG test vectors and batch-partition invariance;
- held-out leakage checks;
- covariance normalization and propagation;
- deliberately corrupted metadata or certificates being rejected.

## Data and results

Prefer reviewable text formats: CSV, JSON, YAML, Markdown, and exact integer coefficient files. Large binary files require a documented reason and storage policy.

Every imported source dataset needs:

- source citation and exact table/equation/location;
- estimator and geometry definition;
- decimal-preserving transcription;
- row-count and sentinel-value tests;
- SHA-256 checksums.

Every stochastic result needs the metadata described in `REPRODUCIBILITY.md`. Commit raw sufficient statistics, not only final roots or fitted coefficients.

Do not overwrite a committed result. Add a correction directory or a new version and link the superseded artifact.

## Experiment design

Before production, freeze:

- the primary hypothesis and alternatives;
- all sizes/geometries and signed orientation order;
- training and held-out partitions;
- seed/counter domains and batching;
- primary estimand and covariance treatment;
- power target and sample-count rule;
- success, failure, and stopping conditions.

Pilot data may be used for variance and power. Evaluation choices must not depend on whether the pilot point estimate favors the hypothesis.

Report signed effects, standard errors, covariance/correlation, condition numbers, and sensitivity or power. Preserve negative results and rejected models.

## Pull requests

Keep source/protocol changes reviewable. Bulk raw results should normally be a separate PR or clearly isolated directory with a manifest.

A PR description should include:

- linked issue and dependency PRs;
- change class and proposed claim level;
- exact commands/tests run;
- provenance and generated-file boundaries;
- compatibility or migration implications;
- known limitations;
- explicit statements of what the change does **not** establish.

Use the repository PR template. Do not mark a PR ready for review while required result files, tests, or metadata are still expected.

## Scientific writing

Use precise language:

- “observed in these sizes” rather than “asymptotic” without an asymptotic test;
- “compatible with” rather than “equal to” for a numerical candidate;
- “candidate operator” rather than “identified operator” until competing sectors are excluded;
- “exact” only for a proved identity, exact arithmetic result, or certified computation.

A result may be valuable because it falsifies a promising route. Do not optimize the narrative by removing failed attempts.

## Review priorities

Reviewers should examine, in order:

1. definition and sign/unit conventions;
2. provenance and frozen design;
3. exact controls and independent implementation checks;
4. covariance and held-out integrity;
5. numerical conditioning and finite-size alternatives;
6. performance and presentation.

Correctness and auditability take precedence over throughput.
