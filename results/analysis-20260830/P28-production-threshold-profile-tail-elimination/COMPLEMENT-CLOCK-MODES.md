# P28 exact-complement clock-mode diagnostic

This final post-reveal compression uses no new Monte Carlo and fits no new tail family.

## Exact transformation

The repository's digital-Alexander filtration identity is

```text
K1^G(pi) + K2^Ghat(reverse(pi)) = N+1,
K2^G(pi) + K1^Ghat(reverse(pi)) = N+1.
```

For continuous priority thresholds this implies reflected density maps.  After each component is
centered and scaled by its own reflected mean and standard deviation,

```text
K1^G right(z) = K2^Ghat left(z),
K1^G left(z)  = K2^Ghat right(z).
```

Thus K1-right can be compared with K2-left as matching versus primal realizations of the same
second-birth tail; K1-left pairs with K2-right similarly.  The identity exchanges `G` and
`Ghat`: it does not itself assert that the two profiles are equal.

For each mapped tail and curvature residual, define

```text
E = (matching + primal)/sqrt(2),
O = (matching - primal)/sqrt(2).
```

The inverse transformation reconstructs all component residuals with maximum error `2.8e-17`.

## Result

| coordinate | global chi-square / df | interpretation |
|---|---:|---|
| complement even | `28,632,665.24 / 48` | dominant common curvature |
| complement odd | `1,892,387.37 / 48` | smaller but resolved mismatch |
| amplitude-free pair closure | `26,071.06 / 32` | common shape rejected |

The closure survival probability has `log10 p approximately -5612`.  Nevertheless the raw
curvature vectors are extremely close to collinear: every paired cosine lies between `0.9992` and
`0.9999`.  Allowing a free amplitude leaves relative shape mismatches of only `1.5–2.1%` on the
mapped-left pair and `2.6–3.9%` on the mapped-right pair.  Production precision resolves these
small deviations decisively.

The lowest smooth curvature mode dominates both sectors:

| coordinate | mode-2 marginal chi-square | mode-3 | mode-4 |
|---|---:|---:|---:|
| even | `5,651,958.77` | `18,404.00` | `3,353.20` |
| odd | `497,064.09` | `19,262.32` | `586.58` |

Marginal mode scores are correlated and do not sum to the full coordinate score.

## Size drift

| N | odd/even chi-square ratio | mapped-left matching/primal amplitude | mapped-right amplitude |
|---:|---:|---:|---:|
| 265 | `0.0736` | `1.84–1.88` | `0.38–0.46` |
| 290 | `0.0685` | `1.81–1.84` | `0.38–0.47` |
| 325 | `0.0733` | `1.78–1.79` | `0.39–0.42` |
| 425 | `0.0487` | `1.67` | `0.45–0.49` |

The mapped-left amplitude and odd/even ratio decrease at the largest size, but no zero-odd closure
is visible.  The `N=290` mapped-left pair inherits the previously recorded K1-right count-gate
boundary and is diagnostic only; the same structure is independently present at `N=265,325,425`.

## Mechanism compression

The complement-related tails are two faces of a **shared dominant curvature direction**, but not
the exact same finite-size mode.  A smaller complement-odd shape component survives even after
amplitude removal.  The data therefore support:

```text
dominant complement-even Alexander/clock large-deviation mode
  + resolved complement-odd subleading correction,
```

rather than either an exact single mode or two independent leading channels.  This statement is
about the measured finite-size threshold profiles and does not identify a continuum operator.
