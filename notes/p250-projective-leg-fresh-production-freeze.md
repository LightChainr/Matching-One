# P250 fresh projective-leg production freeze

The 2k smoke is used only for denominator and eight-real support power.  The
frozen grid is `5k,10k,20k,40k`, with production thresholds

```text
weakest d1/d2 charged pair |z| >= 5,
minimum d1/d2 support power >= 0.90 at alpha=0.01.
```

The noncentrality estimate is the explicit exploratory alternative
`max(chi2-df,0)`, scaled linearly with sample count.  At 5k, the weakest d2
pair is projected at only `4.585 sigma` and d1 support power is `0.716`.  At
10k they become `6.484 sigma` and `0.986`, while d2 support power is already
effectively one.  Therefore 10k is the first qualifying grid point; 20k and
40k are deliberately not selected.

The fresh stream uses seed `25033433720260930`, counters `[0,10000)`, 50
batches and the unchanged projective-leg operator.  The support scorer first
requires, separately at d1 and d2:

1. all four pair denominators at `|z|>=5`;
2. a nondegenerate eight-real cubic covariance;
3. rejection of the zero cubic vector at `p<0.01`.

Only if all six conditions pass does the code enter the phase branch and form
the joint four-real d1/d2 cross-product closure.  A failed support stage emits
a locked record and never computes a phase statistic.

The run is assigned to `Huawei-CodeBuddy-XPk2PZ`.  Results are revealed once,
after completion and checksum transfer.
