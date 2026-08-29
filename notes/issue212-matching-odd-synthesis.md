# Issue #212: matching-odd signal-existence synthesis

## Boundary

This is a post-reveal synthesis of two evidence-ledger rows that were already
registered as primary evidence. It is not a new preregistration and must not be
appended to the ledger as another primary row.

The scorer reads exactly `issue43_n185_n265_deltaM` and `issue57_norm5`. It
requires both to be scored `primary` rows in the `matching_odd -> matching_odd`
channel and freezes their distinct raw-data groups. Within each block it reads
only the registered `zero_effect` score and that block's fixed H4 prediction.
Roots, derivatives, Krawtchouk modes, P49, P50, and other correlated views do
not enter.

## Result

Under the Issue #212 block-diagonal independence contract, the zero-effect
scores add to

```text
chi2 = 31.18573555150965, dof = 4, p = 2.805595267905808e-6.
```

The corresponding fixed H4 predictions add to

```text
chi2 = 3.4622795373044296, dof = 4, p = 0.48363695393249573.
```

The fixed-H4 and zero-effect predictive NLPDs are respectively
`-35.7946059274312` and `-21.988187137702105`, hence
`Delta NLPD(H4-zero) = -13.806418789729096` nats. Lower NLPD favors fixed H4.
NLPD is the appropriate direct predictive comparison here; a raw chi-square
difference is not treated as a likelihood-ratio statistic because the frozen
predictive covariances need not agree.

Thus the norm-5 children alone remain compatible with zero, while the two
independent matching-odd target blocks jointly reject a global zero-effect
description and remain compatible with their fixed H4 predictions. This moves
the scientific bottleneck from signal existence to mechanism and transfer.

## Reproduction

```bash
python3 scripts/score_matching_odd_synthesis.py \
  results/evidence-ledger/latest.json \
  --output results/evidence-ledger/issue212-matching-odd-synthesis.json
python3 tests/test_matching_odd_synthesis.py
```

The JSON records the source ledger SHA-256, the exact selection contract, the
two block-level inputs, the additive scores, and explicit non-primary
governance flags.
