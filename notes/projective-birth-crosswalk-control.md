# Exact projective-birth crosswalk control

This control addresses one bounded part of Issue 439: reconstructing the
rank-sector observables from typed, same-stream projective birth rows. It does
not import or score production data.

Each row carries exact birth times `tau1 <= tau2`. A row is either:

- `direct_rank2`, with `tau1 == tau2` and no rank-one plateau line; or
- `plateau`, with `tau1 < tau2` and a nonempty projective line label.

Births are inclusive: a birth at `tau` is present at threshold `p` when
`tau <= p`. The rank is therefore zero before `tau1`, one on
`[tau1, tau2)`, and two from `tau2` onward.

For every threshold, the exact integer counts obey

```text
P0 + P1 + P2 = n
F1 = P1 + P2
F2 = P2
P1 = sum(line plateau counts)
M = P2/n - P0/n = F1/n + F2/n - 1.
```

The committed certificate uses six deterministic synthetic rows, including
both direct and plateau cases. Its five threshold snapshots are:

| p | (P0,P1,P2) | (F1,F2) | active plateau lines | M |
|---|---|---|---|---|
| 0 | (5,1,0) | (1,0) | L0:1 | -5/6 |
| 1/4 | (4,2,0) | (2,0) | L0:2 | -2/3 |
| 1/2 | (2,2,2) | (4,2) | L0:1, L1:1 | 0 |
| 3/4 | (0,2,4) | (6,4) | L1:1, L2:1 | 2/3 |
| 1 | (0,0,6) | (6,6) | none | 1 |

The parser rejects floats, incomplete or extra fields, invalid birth order,
out-of-range times, inconsistent direct/plateau tags, and empty inputs.

## Boundary

This certificate establishes only the exact row-level crosswalk and closure
identities. It does not establish availability or correctness of raw archives,
aligned covariance or jackknife estimates, wedge/common-ray/transfer fits,
the direct-02 decomposition, N=1360 forecasts, or any production/physical
claim. Issue 439 must remain open for those tasks.
