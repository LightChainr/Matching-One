## Purpose

<!-- What question or maintenance need does this PR address? Link the issue. -->

Closes/advances: #

Dependency PRs: none

## Change class

- [ ] Documentation/governance only
- [ ] Reference or production code
- [ ] Dataset/provenance import
- [ ] Frozen experiment protocol
- [ ] Raw/derived result archive
- [ ] Scientific conclusion or claim-ledger change
- [ ] Correction/retraction

Proposed claim level: `C0 / C1 / C2 / C3 / C4 / C5 / not applicable`

## What changed

<!-- Separate source, protocol, raw data, derived data, and narrative files. -->

## Frozen design and anti-leakage

<!-- For experiments: geometry/size set, orientation order, training/held-out split, seed/counter domains, sample-count rule, primary statistic, covariance model, and acceptance criteria. State what was known before target evaluation. -->

Not applicable / details:

## Provenance

- Full source commit:
- Dirty tree: `false / true / not applicable`
- Source-file hash:
- Executable hash:
- Environment/compiler/interpreter:
- Command record:
- Input and output checksums:
- RNG/domain/batch record:

## Validation

Commands run:

```text

```

- [ ] Unit/regression tests pass
- [ ] Exact or independent reference checks pass where required
- [ ] One-thread/multi-thread and batch-partition invariance checked where promised
- [ ] Invalid/corrupted input rejection tested where relevant
- [ ] Generated files can be reproduced from committed inputs

## Statistical review

- [ ] Signed effects and uncertainties are reported
- [ ] Shared-randomness covariance is retained and propagated
- [ ] Conditioning and model instability are reported
- [ ] Held-out data did not affect model selection
- [ ] Null/negative results and failed gates are preserved
- [ ] Power or sensitivity is stated for unresolved effects

Not applicable / details:

## Scientific boundary

Direct observations:

Model-dependent deductions:

This PR does **not** establish:

## Compatibility and migration

<!-- Frozen conventions, units, signs, rank definitions, file formats, or downstream scripts affected. -->

## Reviewer checklist

- [ ] Definitions, units, and sign/order conventions are explicit
- [ ] Source/protocol/result boundaries are reviewable
- [ ] Claim language matches `GOVERNANCE.md`
- [ ] `docs/STATUS.md` is updated if claim strength changes
- [ ] Historical artifacts are preserved rather than overwritten
- [ ] CI passes
