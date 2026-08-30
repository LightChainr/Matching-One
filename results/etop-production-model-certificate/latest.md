# Production E_top model-elimination certificate

The exact state basis is `A_top=P2-P0`, `E_top=P0+P2=1-P1`.
Each fixed model uses a separate 99% familywise Gaussian-Bonferroni outer confidence set over eight high-statistics production datasets; no cross-dataset p-values are pooled.

| model | equation | decision | incompatible production rows |
|---|---|---|---|
| M_ETOP_ZERO | E_top=0 | eliminated | P49-N130, P50-N145, P49-N170, P43-N185, P43-N265 |
| M_F2_ZERO | E_top=-1 A_top | eliminated | P50-N145 |
| M_F1_ZERO | E_top=1 A_top | eliminated | P49-N130, P50-N145, P49-N170, P43-N185, P43-N265, P50-N290 |
| M_COMMON_STATE_LINE | E_top=r A_top, free r | not_eliminated | common r=[-0.9305631428471134, -0.4551779794216242] |

## Interpretation

Production data eliminate a pure Alexander-odd state response (`E_top=0`) and both exact endpoint-cancellation lines (`F1=0` and `F2=0`) under their declared outer confidence sets. A one-dimensional state line with a free common mixing ratio is not eliminated, so the result establishes a required even component without yet proving a two-dimensional continuum module.

The `F2=0` exclusion is the narrowest of the three and should be read from its saved simultaneous interval, not as a generic absence of finite-size cancellation.

## Boundary

This is a machine-checkable certificate relative to the declared Gaussian first-order outer confidence sets. It is not an exact finite-replica theorem, a continuum field identity, an exponent fit, or an independence-weighted meta-analysis.
