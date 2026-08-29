# Pre-registration and trial-ledger audit gate

This artifact freezes a content-addressed statistical plan. It contains only synthetic
validator fixtures and makes no empirical percolation claim.

## Frozen plan

- plan SHA-256: `8beeb3977457f391993c19012b14bb264f76744052c616fdb9eaf76b09d7e29a`;
- familywise alpha: `1/1000000`;
- two sides, at most `3` attempts per side;
- `400` trials per attempt; per-run alpha `1/6000000`;
- exact acceptance cutoff: `373` successes.

## Validator controls

The synthetic six-record fixture has `2` accepted records and exercises both sides. Tests
also require rejection of plan-digest tampering, duplicate side/attempt pairs, duplicate
domains or data digests, exploration-data reuse, graph/side mismatches, and invalid counts.

## Boundary

A nonempty independence attestation is provenance, not proof. This validator checks metadata
and exact arithmetic only; genuine randomness and IID Bernoulli sampling remain external
conditions for any future confidence statement.
