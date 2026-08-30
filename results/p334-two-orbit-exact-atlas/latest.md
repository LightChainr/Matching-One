# Bounded exact two-orbit HNF atlas

Status: `bounded_exact_two_orbit_atlas`.

HNF index 4..12; four distinct face corners; connected primal graph; no physical quarter-turn lattice symmetry; exactly two projective line orbits under the exact D4 stabilizer; distinct nonzero effective orbit characters. Only geometry, line support and chi4 enter.

Scanned 119 HNFs and included 6. The frozen close threshold is `0.0409499231130266717212867090127`.
The revealed atlas contains only 1024 subset states and 3840 directed boundary edges in total.

| HNF | N | orbit roots | separation | close | Gram | topology below/between/above |
|---|---:|---|---:|---|---:|---|
| [[7, 2], [0, 1]] | 7 | 0.592783237894885157837638571036, 0.571428571428571428571428571429 | 0.0213546664663137292662099996079 | True | 527/625 | reinforce/cancel/reinforce |
| [[7, 3], [0, 1]] | 7 | 0.571428571428571428571428571429, 0.592783237894885157837638571036 | 0.0213546664663137292662099996079 | True | 527/625 | reinforce/cancel/reinforce |
| [[7, 4], [0, 1]] | 7 | 0.592783237894885157837638571036, 0.571428571428571428571428571429 | 0.0213546664663137292662099996079 | True | 527/625 | reinforce/cancel/reinforce |
| [[7, 5], [0, 1]] | 7 | 0.592783237894885157837638571036, 0.571428571428571428571428571429 | 0.0213546664663137292662099996079 | True | 527/625 | reinforce/cancel/reinforce |
| [[8, 3], [0, 1]] | 8 | 0.571428571428571428571428571429, 0.622035526990772772785483463766 | 0.0506069555622013442140548923372 | False | -7/25 | cancel/reinforce/cancel |
| [[8, 5], [0, 1]] | 8 | 0.622035526990772772785483463766, 0.571428571428571428571428571429 | 0.0506069555622013442140548923372 | False | -7/25 | cancel/reinforce/cancel |

## Atlas answer

One simple balance root per orbit is universal in the bounded atlas, but closeness under the inherited N13/N17 envelope is not: the N8 class is the minimal counterexample. Cooperation topology stratifies exactly by character-Gram sign across every full p curve.

- one simple root per orbit: 6/6;
- close under frozen envelope: 4/6;
- minimal closeness counterexample: `hnf-8-3-0-1`;
- exact mechanism signatures after HNF symmetry copies: 2;
- all coefficientwise compression and cross-term identities: `True`.

Exact theorem:

```text
Re[(chi1 J1) conjugate(chi2 J2)] = Gram(chi1,chi2) J1 J2
```

Positive Gram produces reinforce/cancel/reinforce across the two simple zeros; negative Gram produces cancel/reinforce/cancel. This holds on every included full curve.

Next prediction: Beyond N12, root multiplicity and separation are dynamical questions, but any two-orbit quotient passing the same compression gate must retain Gram-sign cooperation topology over every root interval.

## Boundary

- The atlas exhausts the declared HNF geometry gate only through index 12.
- HNF variants sharing one exact signature are symmetry copies, not independent evidence.
- Root closeness is judged by the frozen parent envelope and is not retuned on this atlas.
- No Monte Carlo sample, Huawei production, new PR, or merge is used.
