# Reproducibility Standard

> **When this applies.** At publication time, and for expensive runs nobody wants to
> repeat. It is not a gate on exploratory work — see `GOVERNANCE.md` §0 and §2.
>
> **While exploring, the reproducibility minimum is two lines:** record the command
> and the seed, and commit the raw sufficient statistics rather than only the fitted
> numbers. Everything below is what you fill in later, for the runs that turned out
> to matter.
>
> The reason the full form exists is that some of it genuinely cannot be
> reconstructed afterwards — a seed, a compiler flag, a dirty tree. Where a field is
> in that class it is marked **irrecoverable** below and is worth recording at the
> time even while exploring. The rest can wait, and should.

This standard applies to computation, imported data, statistical analysis, and claim-bearing reports in Matching One.

## 1. Reproduction target

A result is reproducible when another person can identify the exact inputs, source, executable or interpreter environment, stochastic domain, commands, and analysis needed to regenerate the reported sufficient statistics and conclusions.

A screenshot, final decimal, or prose description alone is not a reproducible artifact.

## 2. Minimum production metadata

For a published or expensive result, record, as applicable — the **irrecoverable**
ones first, since those are the only fields worth stopping for at the time:

```text
full_git_commit        irrecoverable
dirty_tree             irrecoverable
command_line           irrecoverable
seed / counter domain  irrecoverable  (see section 3)
compiler_flags         irrecoverable in practice
```

The remainder can be reconstructed later from the commit, and should be filled in
when the result is written up rather than when it is produced:

```text
result_id
created_utc
repository
full_git_commit
source_paths
source_sha256
executable_sha256
dirty_tree
os_and_kernel
machine_or_instance
cpu_gpu_and_memory
compiler_or_interpreter
compiler_flags
dependency_versions
command_line
working_directory
input_paths_and_hashes
output_paths_and_hashes
```

For a confirmation run, `dirty_tree` must be false. If historical work was run from a dirty tree, preserve it, identify the source-file hash, and label the result as provenance-limited until clean replay.

## 3. Stochastic computation

Record:

```text
rng_algorithm
rng_version_or_implementation_hash
domain_separation_scheme
seed
counter_first
counter_last_exclusive
replica_or_permutation_definition
batch_count
samples_per_batch
thread_count
scheduling_policy
```

The domain key must explicitly state which of model, geometry, orientation, size, observable, and engine tag are mixed into the stream. Cross-size or cross-geometry coupling may be intentional, but it must be declared and its covariance retained.

Required checks where applicable:

- deterministic test vectors;
- one-thread versus multi-thread equality;
- batch concatenation/partition invariance;
- non-overlapping confirmation domains;
- integer overflow analysis;
- exact or oracle agreement on tiny systems.

## 4. Sufficient statistics

Store enough raw aggregate data to recompute the reported analysis without rerunning the simulation.

For threshold-rank work this normally includes per-batch integer `K_minus` and `K_plus` histograms, sample counts, aligned batch identifiers, joint moments, and the cross-orientation information required for covariance-aware differences and roots.

For fixed-parameter work this normally includes per-batch indicator or score sums for every retained orientation/sector/channel, not only pooled means.

Do not store only normalized floating-point probabilities when exact integer counts are available.

## 5. Imported data

Every source table or sequence must have:

- full citation;
- exact page, table, equation, repository, or supplementary-file location;
- estimator, topology, boundary condition, and unit definition;
- decimal strings preserved exactly as published;
- transcription date and method;
- row count, sentinel values, and independent check;
- file hash.

Conflicting published estimates remain separate method-specific records. They must not be collapsed into one “accepted” decimal without a documented decision.

## 6. Blind, frozen, and held-out analyses

The repository must retain the artifact that freezes:

- model family;
- training/pilot/evaluation split;
- target sizes or geometries;
- sign/orientation order;
- primary statistic;
- covariance model;
- sample-count or stopping rule;
- acceptance threshold.

The freeze artifact must precede target evaluation. A commit timestamp alone is not sufficient when the target data were already visible to the analyst; the report should state what was known at freeze time.

Changing held-out values must not alter model selection or fitted training objects. This property should be regression-tested when practical.

## 7. Covariance and uncertainty

Use paired or joint batch statistics whenever observables share randomness. Report whether covariance is measured, analytically known, approximated, or ignored.

A fit on correlated training data and a held-out score sharing the same random stream must propagate:

- training covariance;
- held-out covariance;
- training/held-out cross-covariance;
- parameter-estimation uncertainty.

Do not convert a deterministic fit residual into a confidence interval for the infinite-size limit. Separate arithmetic precision, Monte Carlo uncertainty, finite-size model uncertainty, and data/provenance uncertainty.

## 8. Exact and certified work

Exact claims must specify the arithmetic domain and normalization. Examples include integer coefficient arrays, rational Virasoro calculations, exhaustive enumeration, interval arithmetic, or machine-checkable certificates.

A certified computation should provide a separate verifier that does not share the discovery/optimization path. The verifier must reject deliberately corrupted inputs.

## 9. Result directory contract

Preferred layout:

```text
results/<campaign>/<task>/
  REPORT.md
  metadata.json
  commands.txt
  environment.txt
  checksums.sha256
  raw/
  derived/
  logs/
```

Small tasks may use a flatter layout, but the distinction between raw and derived files must remain clear.

`REPORT.md` should contain:

1. question and frozen protocol;
2. source and environment;
3. exact controls;
4. primary result with covariance;
5. held-out or confirmation score;
6. negative results and failed gates;
7. scientific claim level;
8. limitations and next discriminator;
9. reproduction commands.

## 10. Corrections

Never overwrite a committed result to make history appear clean.

A correction must include:

- the affected artifact and claim;
- the defect category;
- whether central values, uncertainties, or interpretation changed;
- corrected files under a new path or version;
- compatibility comparison;
- updated checksums and claim ledger.

## 11. Independent replication

The strongest replication changes more than the random seed. For high-risk or C4-level claims, seek at least one of:

- independent topology or estimator implementation;
- independent derivation;
- separate machine/compiler;
- separate data transcription;
- different geometry or microscopic control with a predicted relation.

Agreement between two analysis scripts consuming the same flawed sufficient statistic is not independent replication.

## 12. Archival boundary

Raw outputs may be retained even when later found underpowered, biased, or provenance-limited. Their reports must make the limitation explicit. Canonical scientific conclusions are governed separately by `docs/STATUS.md`.
