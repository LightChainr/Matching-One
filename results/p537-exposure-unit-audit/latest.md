# P537 exposure-unit correction sidecar

## Decision

The N25 scorer stores the selected positive exposure as `P`; the N65 scorer
stores it as `P/N`.  The historical scale diagnostic divided these coordinates
directly.  This sidecar puts both sizes into one declared convention before
forming conditional signed density.  It does not overwrite the frozen result.

## Entry / double-contact cell

| convention | N25 exposure | N65 exposure | exposure decay | density decay | signed decay |
|:--|--:|--:|--:|--:|--:|
| historical mixed `P` vs `P/N` | 0.0326940900527 | 0.000272689693714 | 5.009479622 | -1.387909616 | 3.621570006 |
| common unweighted `P` | 0.0326940900527 | 0.0177248300914 | 0.640733134 | 2.980836872 | 3.621570006 |
| common source-weighted `P/N` | 0.00130776360211 | 0.000272689693714 | 1.640733134 | 1.980836872 | 3.621570006 |

For common unweighted `P`, the conditional signed densities are
`-8.98415542361e-05` at N25 and
`-5.20606276782e-06` at N65.  The repaired
finite-pair decomposition is therefore approximately
`0.640733 + 2.980837 = 3.621570`.

## What remains frozen

- N65 transmission decision: `CONTACT_FUSION_COMPLETION_TRANSMITS`.
- N65 determinant: `-8.68821605512176e-14`.
- Full `J65`: `-0.00162250988939 +/- 0.000185530082422`.
- The signed N25/N65 matrices and the post-reveal shape score are numerically unchanged.

The repair changes the exposure/density attribution, not signed mass.  The old
`[3, 29/8, 3, 3]` score (`Q=0.636435678`, nominal
`p=0.958930060`) remains a post-reveal shape fingerprint.
Its additional `5/8` is no longer attributed to conditional signed strength and
does not identify a field or occurrence-frequency exponent.

## Mechanism language

- A nonzero stage-by-contact determinant establishes a nonseparable signed table.
  Calling it a physical commutator still requires independently defined operations
  `F`, `B`, `FB`, and `BF`; the present construction is an algebraic encoding candidate.
- The selected cells account for `2.550516%` of
  the full `J65` point estimate.  The rest is named only
  `complement_of_selected`, not a spatially nonlocal mechanism or causal share.
- Exact joint uncertainty for that share awaits the shared full/selected delete-one
  factors; none is invented here.

## Provenance and boundary

- Manifest: `analysis/p537_exposure_unit_audit_manifest.yaml` (`b1fd74547cda38b9bea348789f74352a673f7c797d1b714975a1e130897233e6`).
- Every input is read from a pinned Git blob and its Git blob SHA-1 is verified.
- Existing fixed data only; no raw TSV replay, new random samples, GPU, cloud job,
  full test suite, or historical-result overwrite.
- These are two-size finite coordinates, not universal exponents.  Priority is
  attention allocation, not permission or a task lock.
