# Issue #43 matching-even channel-map erratum

## Summary

The prospective N=185/N=265 data in #108 do **not** support the published interpretation that the matching-even sector reverses sign relative to the P31 frozen amplitude.

The apparent sign reversal comes from comparing different torus wrapping channels:

- the frozen source amplitude was fitted from P31 `channel=either, sector=even`;
- the threshold-rank engine stores rank-2 **cross** thresholds;
- the Issue #43 full-curve scorer therefore observes `cross/even`.

For the complementary primal/matching construction on the torus, cross and either obey the exact topological exchange

`S_either = 1 - S_cross`,

so orientation contrasts satisfy

`DeltaS_cross = -DeltaS_either`.

No physical amplitude changes sign under this conversion; the sign is part of the channel definition.

## Source evidence available before the N185/N265 targets

The frozen artifact states explicitly:

```text
matching_even:
  source: P31
  channel: either
  sector: even
  A_I: +0.010603216462677735
```

The P31 confirmation CSV itself contains, for the same geometries and same seed, opposite cross/either even amplitudes. For example:

```text
N=65  cross/even  N*P4 = -0.010259341579861...
N=65  either/even N*P4 = +0.010259341579861...

N=85  cross/even  N*P4 = -0.013005642849392...
N=85  either/even N*P4 = +0.013005642849392...
```

`tests/test_issue43_cross_either_correction.py` checks this relation across every common P31 size in the committed source table.

## Corrected no-refit score

The original frozen either-channel means were

```text
N=185  +0.00006752163745881449
N=265  +0.00006891944697034459
```

The correct cross-channel predictions are therefore

```text
N=185  -0.00006752163745881449
N=265  -0.00006891944697034459
```

The #108 observations are

```text
N=185  -0.000060815376233385734 +/- 0.000008089565565582085
N=265  -0.00007024950784522366  +/- 0.000009385620077608982
```

Using the same fully correlated frozen source-amplitude uncertainty as the original scorer gives

```text
residual z:  +0.6672, -0.1189
chi-square:   0.5700315436 / 2
```

This uses zero target-fit parameters. The exponent, amplitude magnitude and source uncertainty remain exactly the preregistered values.

## Scientific consequence

The matching-even result should be reclassified from

> prospective sign reversal / failure of the x=4 even sector

to

> prospective **compatibility** with the frozen matching-even amplitude after the exact either-to-cross channel map.

The matching-odd `DeltaM` score is unaffected.

Therefore #108 remains a valid prospective target reveal, but its reported two-sector *interpretation* needs an erratum. The original files should remain in history; a corrected score/report should be added rather than silently rewriting the raw result.

## Governance consequence

Future prediction artifacts and scorers must include explicit fields for:

- source wrapping channel;
- target wrapping channel;
- exact channel map, if any;
- raw versus normalized angular convention.

A hash match alone is insufficient when the source artifact and target sufficient statistic live in different but exactly related topological channels.
