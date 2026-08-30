# P429 shared-prefix branching continuation

## Acquisition

At the fixed rank-one checkpoint, one ordinary permutation update is shared.
If it remains rank one, the exact successor is cloned and two independently
tagged counter streams select one suffix site each.  Common absorption scores
both suffixes zero.  The same prefix, common update and suffix counters are
shared across the paired orientations.

The 20k gate passed (`z=40.45` at N325, `z=34.84` at N425), triggering the
predeclared fresh 100k block.  Production used 100 batches per size:

| size | host | rank-one risk rate | common-safe defined rate |
|---|---|---:|---:|
| N325 first / second | HZsCM6 | 0.4549 / 0.4569 | 0.9370 / 0.9369 |
| N425 first / second | TgFr7R | 0.4476 / 0.4461 | 0.9479 / 0.9482 |

Both remote self-tests passed and both stderr logs are empty.

## Frozen primary

In `(N325 first, N325 second, N425 first, N425 second)` order,

- `b2_survival_estimate=E[(y1+y2)/2]` is
  `(0.876094, 0.875164, 0.896875, 0.897138)`;
- `clone_dependence_gap=E[y1*y2]-E[y1]E[y2]` is
  `(0.053146, 0.053115, 0.045441, 0.044560)`.

The covariance-weighted common gaps are

- N325: `0.053131 +/- 0.000622`, `z=85.47`;
- N425: `0.044955 +/- 0.000611`, `z=73.53`.

The complete eight-coordinate delete-one-batch covariance is in
`production_score.json`.

## Mechanism decomposition

The large primary gap is mostly a structural shared-gate effect: if the common
update is absorbed, both clone responses are simultaneously zero.  After
conditioning on a common-safe successor, the residual heterogeneity is much
smaller but remains positive:

| size | common-safe successor gap | SE | z |
|---|---:|---:|---:|
| N325 | 0.001610 | 0.000224 | 7.18 |
| N425 | 0.000928 | 0.000211 | 4.41 |

Across individual environments, the successor component accounts for only
`1.4%--2.9%` of the unconditional gap.  Thus `97%--99%` of the large headline
gap is the common absorption gate and must not be called memory.  The small
remainder is the production-relevant branching-sensitive successor
heterogeneity.

## Scientific boundary

The exact N16 result in #429 already proves that a complete unbranched survival
signature need not be an autonomous Markov state.  This production block shows
that a small common-safe branching remainder persists at N325 and N425 under a
frozen same-prefix experiment.  It does **not** condition on the full microscopic
survival vector per configuration, so it supports the exact mechanism rather
than independently proving scaled non-Markov closure.  It identifies neither a
continuum memory field nor a scale exponent.
