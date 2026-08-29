# Content-addressed pre-registration for the finite-size confidence gate

Status: protocol-only audit slice for Issue 112. The checked-in records are synthetic validator
fixtures, not percolation observations.

## Frozen statistical plan

The plan imports the exact constants already verified by the confidence-gate oracle:

- familywise error probability `1/1,000,000`;
- upper and lower sides;
- at most three pre-indexed attempts per side;
- 400 trials per attempt;
- null block-event probability `8457/10000`;
- Bonferroni per-run alpha `1/6,000,000`;
- acceptance at 373 or more successes, using the exact rational binomial tail.

The full plan is serialized as canonical JSON and addressed by SHA-256. Every final record carries
that digest. Changing any plan constant, forbidden exploration digest, stopping rule, or record
contract invalidates previously prepared records.

## Fail-closed ledger rules

Each record declares its side, graph, attempt index, tested parameter, trial and success counts,
stream domain, data digest, and an independence attestation. The validator rejects:

- a plan-digest mismatch;
- anything except phase `final`;
- a square/matching graph inconsistent with the upper/lower side;
- an attempt outside `1..3`, or a repeated side/attempt pair;
- a trial count other than 400 or a success count outside `0..400`;
- repeated record IDs, stream domains, or data digests;
- a final digest listed as exploration data;
- disagreement between the exact binomial tail and the frozen cutoff.

Stopping after an accepted pre-indexed attempt is compatible with the six-test Bonferroni budget.
Inventing a seventh test, changing `N`, or choosing a test after seeing final outcomes is outside the
frozen protocol.

## What the digest cannot prove

Content addressing detects plan and record substitution; it does not prove that occupation fields
were generated independently or from genuine Bernoulli randomness. The required attestation is a
provenance statement whose truth must be established outside this validator. Domain separation is
also an engineering safeguard, not a mathematical randomness proof.

Consequently this artifact proves no event probability and no bound on square-site `p_c`. It prepares
an auditable envelope into which future, genuinely fresh final trials could be placed.
