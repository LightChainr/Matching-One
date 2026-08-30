# P28 pure-tail rejection diagnostic

This is a post-reveal diagnostic of the frozen production result.  It does not alter the original
three-model decision and launches no simulation.

## Attribution

The principal conclusion is structural: every size, orientation, and side already had its own
intercept and decay constant in the frozen test.  Cross-size drift, orientation amplitudes, and
left/right amplitude asymmetry were nuisance-saturated before the lack-of-fit contrasts were
formed.  The rejection is therefore caused by **single-tail internal curvature**.

| diagnostic partition | covariance-aware result |
|---|---:|
| left-tail marginal chi-square | `2,758,737.89` |
| right-tail marginal chi-square | `3,525,748.88` |
| first-orientation marginal chi-square | `2,996,561.11` |
| second-orientation marginal chi-square | `3,158,304.45` |
| lowest curvature-mode marginal chi-square | `1,439,918.77` |
| next curvature mode | `13,613.10` |
| highest curvature mode | `4,260.50` |

Marginal partitions under correlated GLS do not sum to the global chi-square.  The machine artifact
also gives an exactly additive but signed point attribution.  It alternates across the interior
tail bins (`2.75`, `3.0`, `3.25`); the outer `3.5` bin is not the driver.

## Effective-exponent map

The local diagnostic uses adjacent log-density slopes and retains its full delete-one covariance.

| N | left beta-effective range | right beta-effective range |
|---:|---:|---:|
| 265 | `1.658–1.701` | `1.945–2.099` |
| 290 | `1.639–1.685` | `1.922–2.057` |
| 325 | `1.626–1.652` | `1.911–2.010` |
| 425 | `1.609–1.639` | `1.832–1.894` |

There is a clear side difference and a downward finite-size drift.  Nevertheless every one of the
16 individual orientation/side curves rejects a constant effective exponent over the window; the
largest survival probability is below `2e-21`.  Thus asymmetry and drift describe the curvature
but do not explain it away.

## Minimal post-reveal correction

Adding one lower-power term to the original target,

```text
log rho(z) = a - c z^(4/3) + d z^(2/3),
```

reduces the global score from `6,280,338.87 / 48` to `49,893.28 / 32`.
This is a large descriptive reduction but still an overwhelming rejection
(`log10 p approximately -10780`).  It is not a renewed model contest and was not used to change
the frozen result.  A successful continuation would need a separately frozen finite-size or
two-birth correction family, still using the existing archives before any new Monte Carlo.
