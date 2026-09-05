# Publication checklist

Everything here is **publication-time work**. None of it gates exploratory work — see
`GOVERNANCE.md` §0 and §2, which is the whole rule set for that.

Reach for this file when there is a specific paper and a specific claim. Then do the
items that the *claim* needs, not all of them: a methodology paper and a numerical
result need different halves of this list.

## Why it is all here and not in the workflow

Verification is cheaper after the fact. By publication time you know which numbers
survived, which is usually one in ten. Verifying during exploration means verifying
all ten, and the nine you abandon take their verification with them.

The second cost is worse: a check written into the workflow runs forever, and every
subsequent contributor reads it as the standard. That is how a research repository
turns into a compliance system — not by one bad decision, but by a hundred
individually reasonable ones.

So: do this work once, late, for the things that are actually being claimed.

## 1. The numbers

- [ ] Every number in the paper is generated, not hand-typed. Tables render from a
      committed artifact and a regression test fails if they drift.
- [ ] Every claimed-exact result was computed in exact arithmetic, with no binary
      floating point anywhere in the claim.
- [ ] Each headline number has one independent confirmation: a second method, a
      closed form, a limiting case, or a second implementation.
- [ ] Transcribed literature values were compared digit by digit against the primary
      source, and that comparison is recorded with its date, coverage and method.
      A checksum is not this check.

## 2. Provenance

- [ ] Source citation precise enough to find the table again: table/equation number,
      version, DOI or arXiv id.
- [ ] Estimator, geometry and normalization defined for every imported quantity.
- [ ] SHA-256 for every canonical machine-readable input, and a check that the pins
      agree with each other and with the files.
- [ ] Frozen inputs and generated outputs clearly separated, with the command that
      regenerates each output recorded and known to work.
- [ ] Environment recorded for anything expensive to rerun: interpreter, compiler,
      flags, versions.

## 3. Design and chronology

Only for claims that rest on it — a C3 or above.

- [ ] Hypothesis, alternatives, sizes, geometries and signed orientation order frozen
      before the target was read, with the freeze visible in history.
- [ ] Training and held-out partitions declared, and model selection demonstrably
      independent of held-out data.
- [ ] Seed and counter domains, batching, and the sample-count rule recorded.
- [ ] Stopping rule declared in advance if optional stopping is claimed.
- [ ] Where a design was *not* frozen, the paper says so and claims C2.

## 4. Statistics

- [ ] Signed effects with uncertainties, not just p-values.
- [ ] Full covariance where views share randomness; correlated views counted once.
- [ ] Condition numbers and conditioning failures reported.
- [ ] Power or sensitivity stated for every unresolved effect — a null with unknown
      sensitivity is not a result.
- [ ] Negative results, failed gates and rejected models retained and reported.

## 5. Reproducibility

- [ ] A reader with the repository can regenerate every figure and table.
- [ ] Regeneration commands committed and verified to run from a clean checkout.
- [ ] Raw sufficient statistics committed, not only fitted coefficients.
- [ ] Randomized results reproduce from the recorded seeds, or the paper says why not.

## 6. Claim language

- [ ] Every claim carries its level (`GOVERNANCE.md` §5) and the language matches:
      "compatible with" not "equal to", "candidate" not "identified", "observed at
      these sizes" not "asymptotic".
- [ ] An explicit statement of what the work does **not** establish. This is usually
      the most-read paragraph and the easiest to get wrong.
- [ ] `docs/STATUS.md` updated if the claim boundary moved.
- [ ] Limitations that a referee will find are named before they find them.

## 7. External inputs

- [ ] Any result taken from an external model, tool or unchecked paper is labelled as
      a theory input, and either verified independently or reported as unverified.
      The ordering is: answer -> prediction -> score against existing data -> claim.
      An external answer is never itself the evidence.

## Release

A paper-oriented release contains a claim-ledger snapshot, source and result hashes,
the major limitations, and enough to reconstruct the reported tables and figures.

A release tag does not upgrade a claim.
