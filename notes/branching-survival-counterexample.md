# Trace-equivalent survival states can branch differently

This exact bounded certificate uses the honest `4 x 4` square torus. Sites are specified by physical coordinates modulo four; row-major masks are included only as reproducibility checks.

The two rank-one states A and B both have `k=8`, primitive line `(1,0)`, and the same complete future rank-one survival counts

```text
b = (1, 7, 18, 20, 8, 0, 0, 0, 0).
```

Therefore their entire unbranched future rank trajectory has the same law. In particular, both have two-step survival probability `9/14`.

They differ after a state-reading operation:

1. choose one common vacant insertion uniformly;
2. if still rank one, clone the successor;
3. choose one independent vacant insertion in each clone;
4. require both clones to remain rank one.

The exact successor distributions and branching scores are

| state | absorbed | safe H2=1 | safe H2=2 | safe H2=3 | successful branches | probability |
|---|---:|---:|---:|---:|---:|---:|
| A | 1 | 3 | 2 | 2 | 190/392 | 95/196 |
| B | 1 | 1 | 6 | 0 | 186/392 | 93/196 |

The gap is exactly `1/98`. Direct enumeration of all `8 * 7 * 7 = 392` branch choices agrees with the calculation from the successor H2 distribution. Cloning immediately at the original checkpoint gives `(7/8)^2=49/64` for both states and does not distinguish them; the common update is essential.

## Boundary

This result certifies one coordinate-defined N=16 counterexample. It does not reproduce the selected-prefix history table, the six-quotient strong-Markov refinement census, predictive-algebra construction, approximate bisimulation, scaled continuation, production acquisition, or a continuum/field interpretation. Issue #429 remains open.
