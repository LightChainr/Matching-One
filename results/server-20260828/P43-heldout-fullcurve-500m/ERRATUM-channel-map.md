# Erratum: Issue #43 matching-even channel map

The original prospective target values remain unchanged. The correction concerns only which frozen source prediction is comparable to the target statistic.

The frozen matching-even source amplitude came from P31 `either/even`. The threshold-rank Issue #43 target is rank-2 `cross/even`. Complementary torus topology gives

```text
DeltaS_cross = -DeltaS_either.
```

Therefore the frozen positive means must be negated before scoring the cross target.

| N | observed cross DeltaS | corrected frozen cross mean | residual / combined SE |
|---:|---:|---:|---:|
| 185 | -6.08153762334e-5 | -6.75216374588e-5 | +0.6672 |
| 265 | -7.02495078452e-5 | -6.89194469703e-5 | -0.1189 |

Using the same fully correlated frozen source-amplitude uncertainty as the original scorer:

```text
corrected DeltaS chi-square = 0.5700315436 / 2 df
survival probability (df=2) = 0.7520023938
```

No target-fit parameter is introduced. The original source amplitude, exponent, target observations and uncertainty model are unchanged.

Consequences:

- withdraw the interpretation that N185/N265 show a physical matching-even sign reversal;
- the prospective matching-even amplitude is instead compatible with its frozen source law after exact channel conversion;
- the matching-odd DeltaM score is unaffected;
- the original #108 artifacts remain historical provenance and should not be silently overwritten.

Reproduce with:

```bash
python3 scripts/score_issue43_cross_either_correction.py \
  --primary results/server-20260828/P43-heldout-fullcurve-500m/analysis/primary_score.json \
  --prediction predictions/two_spin4_heldout_20260828.yaml \
  --json /tmp/channel_map_corrected.json

python3 -m unittest tests.test_issue43_cross_either_correction -v
```
