# P250 radius-five morphism result

## Result

The fresh radius-five shell resolves the ambiguity left by the 80k stream under
the frozen alpha `0.01` decision rule.  Both hand-specific quadratic
annihilators extend to the new degree-five rows.  Identity plus conjugation
fails, while exactly one of the four frozen Alexander conventions survives:
`Alexander_R2_conjugation`.

| gate or map | chi-square / df | finite-batch p | frozen decision |
|---|---:|---:|---|
| plus independent extension | `0.00638 / 16` | `1.0000` | survives |
| minus independent extension | `0.00274 / 16` | `1.0000` | survives |
| identity + conjugation | `32.0406 / 10` | `0.000728` | rejected |
| Alexander R0 + conjugation | `35.8500 / 10` | `0.000192` | rejected |
| Alexander R1 + conjugation | `30.8953 / 10` | `0.001080` | rejected |
| Alexander R2 + conjugation | `23.2988 / 10` | `0.013379` | survives |
| Alexander R3 + conjugation | `24.3110 / 10` | `0.009686` | rejected |

Thus the frozen scorer returns
`Alexander_family_selected_over_identity`.  This is not best-p selection: all
five maps were judged against the same preregistered threshold, and the
survivor list happens to contain only R2.

## What changed

The earlier radius-four archive admitted identity plus conjugation and all four
Alexander reflection conventions.  The 20-point radius-five shell adds the
minimum degree-five moments that make the truncated morphism identifiable.
The new result says the two hand sectors are not related by the naive identity
intertwiner.  They remain compatible with independent five-state quotients,
and their unique surviving parameter-free cross-hand line map is the
Alexander-reflected, twice-rotated, conjugated convention.

This is direct evidence that the previous raw shared-rank failure was a wrong
cross-hand identification, rather than failure of the within-hand five-state
description.  It also selects a concrete convention for the next joint
`Tx,Ty` or higher-degree flat-extension calculation.

## Boundary

R2 has `p=0.01338`, while R3 is immediately below the gate at `p=0.00969`.
The frozen rule makes the decision unambiguous, but the numerical margin is
thin.  The result identifies a degree-two annihilator-line morphism on the
available radius-five moment domain; it does not prove a complete transfer
algebra, a finite-quotient graph isomorphism, or a continuum operator identity.

## Execution and chronology

Huawei-CodeBuddy-XPk2PZ ran the exact frozen `1,200,000` replicas in counters
`[0,1200000)`, seed `25050510120261250`, 400 batches and 16 workers.  It exited
zero after 869 seconds with empty stderr.  Raw hashes were committed before
scoring in `ad162d3`.

The first scorer invocation stopped before computing any statistic because an
empty fresh-row append had NumPy shape `(0,)` instead of `(0,6)`.  Commit
`dfcfc53` fixes only that shape and adds a regression test; candidates,
statistics, alpha, and decision logic are unchanged.  The fixed invocation is
the single completed scientific score.

## Scientific card

- Mechanism space changed: identity cross-hand identification is removed;
  Alexander R2 plus conjugation is selected within the frozen parameter-free
  family.
- Not proved: a full two-generator state model, exact graph morphism, or
  continuum OPE statement.
- Observer/sector/source/geometry: fixed-p projective-leg Z5 charged two-point
  moments; plus/minus Hecke hands; charges 1 and 2; norm-505 radius-five shell.
- Dependency group: the 80k N505 bivariate archive plus this independent 1.2M
  degree-five stream.
- Next lift: use the selected R2 convention as a fixed gauge in the minimal
  higher-degree flat-extension or shared `Tx,Ty` realization; do not reopen the
  five-map vote on the same moments.
