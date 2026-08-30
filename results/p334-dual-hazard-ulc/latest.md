# Direct aggregate-current route to fixed-line ULC

## Exact conditional lemma

Write `A_k=|F_k|`, `beta_k=B_(k-1)/(k A_k)` for the normalized rank-zero birth boundary, and `xi_k=X_k/((N-k) A_k)` for the normalized rank-two exit boundary. Internal-edge double counting gives

`q_(k+1)/q_k = (1-xi_k)/(1-beta_(k+1)).`

Complement duality swaps birth and exit and reverses the layer: `beta_P(k)=xi_M(N-k)` and `xi_P(k)=beta_M(N-k)`. Therefore nondecreasing `xi` on both primal and matching carriers is sufficient for ULC. This is an exact conditional lemma, not yet a general proof of exit-hazard monotonicity.

## What rank monotonicity proves, and the precise missing term

Along every internal fixed-line edge, exit-pivotal sites nest upward and birth-pivotal sites nest downward. Thus the edge-weighted normalized exit hazard is nondecreasing. Passing from edge-weighting to uniform layer-weighting introduces exactly

`xi_(k+1)-xi_k = Delta_edge + Cov_k(u,h_x)/E_k[u] - Cov_(k+1)(d,h_x)/E_(k+1)[d]`,

where `u,d` are internal up/down degrees and `h_x` is local exit hazard. `Delta_edge>=0` is proved pointwise. The covariance correction has no fixed sign, so the missing aggregate lemma is exactly that the edge slack dominates its negative part.

## Bounded exact audit

Across 83 HNFs through N=12, complement line/degree duality passes on 59922 primal rank-one states with zero failures.

| carrier | fixed lines | adjacent pairs | exit hazard monotone | birth hazard monotone | strict ULC sequences | negative bias corrections |
|---|---:|---:|---:|---:|---:|---:|
| primal | 240 | 492 | 240 | 240 | 240 | 208 |
| matching | 240 | 492 | 240 | 240 | 240 | 276 |

The most negative exit degree-bias correction is -1/12 at `[[2, 0], [0, 3]]`, carrier `matching`, line `[1, 0]`, lower layer 2; edge slack 5/12 leaves the positive uniform increment 1/3.

## Status

- **Proved:** current identity, ratio identity, complement birth/exit swap, and edge-coupled pivotal nesting.
- **Exact finite evidence:** uniform exit hazard is nondecreasing on all 480 primal/matching line sequences through N=12, hence every audited q sequence is ULC.
- **Not proved:** the covariance domination inequality for arbitrary quotient size.
- **Revised proof target:** an aggregate two-step path or boundary double count for `Delta_edge + degree_bias >= 0`; statewise matching is neither required nor possible.
