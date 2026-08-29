# Post-reveal rho-child residual handoff

The frozen E4+E6 relation fails coherently, not at one bad child.

| child | residual (re, im) | chi-square contribution |
|---|---:|---:|
| 2omega | ['-0.00050167546103909842897', '0.000016583537475307391022'] | 1.8797962130414171458 |
| omega_over_2 | ['0.00025356196773282770399', '0.00040246736459129475727'] | 1.8476881818063176458 |
| omega_plus_1_over_2 | ['0.00026759278601981573154', '-0.00051770657594355970997'] | 2.3370451865515898102 |

The residual DFT is

- r=0: `['6.4930975711816688533e-6', '-0.000032885224625652520559']`, magnitude `0.00003352011805988320573`
- r=1: `['0.000011547056829398409834', '0.000028784729407951586566']`, magnitude `0.000031014434841059885941`
- r=2: `['-0.00051971561543967850766', '0.000020684032693008325016']`, magnitude `0.00052012705192124859757`

The r2 magnitude is `15.516862171906708319` times the largest other residual row.
Thus the frozen two-character function misses one coherent r2-shaped complex relation; this is distinct from the raw-data statement that standalone r2 remains unresolved.

## Next recognition experiment

Use the three degree-2 children of the opposite-side Pell N30 parent:

- `[[6,6],[0,10]]`
- `[[12,3],[0,5]]`
- `[[12,9],[0,5]]`

Within each triple construct the exact complex annihilator `w` of E4 and E6, then compare the phase of `(w dot delta_H4)/(w dot E4_squared)` between N60 and N112. This removes overall normalization. E4-squared completion predicts phase preservation; signed Pell leakage predicts reversal or loss of coherence. Production is not authorized by this post-reveal note.
