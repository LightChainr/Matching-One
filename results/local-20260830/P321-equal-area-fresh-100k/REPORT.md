# P321 fresh 100k multiscale score

Status: **fresh preregistered decision scored once after all three scales
completed**.

## Decision

The frozen scale law passes:

```text
fixed N^-2/N^-3 scale fit: chi2 = 2.0002264951440183 / 4
p = 0.735717220889469
frozen alpha=.01 critical value = 13.276704135987622
```

The frozen conditional ordinary-thermal-Q4 E4 shape is **not rejected**:

```text
E4 residual score: chi2 = 4.027463739637316 / 3
p = 0.2585137726342605
frozen alpha=.01 critical value = 11.34486673014437
```

Thus the scale gate is open and the primary E4 statistic remains well below
its rejection boundary.  The exploratory 20k tension (`7.531750/3`) did not
strengthen on the fresh streams.  The fresh score is also far below the
pilot-planning alternative's expected score `25.65875`; this is a comparison
to the frozen power calculation, not a separately preregistered alternative
hypothesis test.

## Acquisition and provenance

- freeze/source commit: `82060874d54c832c9bc000d8f16f52104e78445e`;
- remote: `Huawei-CodeBuddy-TgFr7R`;
- remote path: `/workspace/p321-fresh-100k`;
- scales: `N=144,576,1296`;
- samples: 100,000 per shape and scale;
- covariance: 50 aligned batches, full `5x5` within-scale covariance;
- RNG: common random numbers across shapes within each N, independently
  domain-separated across N;
- scoring: `scripts/score_p321_equal_area_rectangles.py` from the frozen
  commit, launched once after every `campaign.json` existed.

The 42 remote acquisition/log files were copied without deleting or changing
the remote source.  The remote and local SHA-256 lists matched file by file
before scoring.  `remote_source_checksums.sha256` preserves that source list;
`checksums.sha256` covers the complete local result package except itself.

The scorer's JSON field remains `status: variance_smoke` because the frozen
implementation labels every campaign with at most 100,000 samples that way.
That legacy schema label was not rewritten after seeing the data; the
scientific lifecycle is the fresh preregistered score frozen in `8206087`.

## Roots

Parentheses contain aligned delete-one-batch standard errors.

| N | rho=1 | rho=16/9 | rho=9/4 | rho=4 | rho=9 diagnostic |
|---:|---:|---:|---:|---:|---:|
| 144 | 0.5927336312 (1.329e-4) | 0.5926751647 (1.583e-4) | 0.5927246342 (1.373e-4) | 0.5927054330 (2.033e-4) | 0.5915776222 (4.568e-4) |
| 576 | 0.5926451222 (8.091e-5) | 0.5928946672 (1.042e-4) | 0.5928160961 (8.714e-5) | 0.5925926606 (1.174e-4) | 0.5928212794 (2.417e-4) |
| 1296 | 0.5927706448 (7.119e-5) | 0.5927006550 (7.804e-5) | 0.5928003201 (6.909e-5) | 0.5926344193 (1.066e-4) | 0.5927271031 (2.662e-4) |

All repeated-square histogram and moment rows are byte-identical at every N.

## Frozen width-normalized shape score

The common fitted root limit is

```text
pc = 0.5927471379931379.
```

The fitted width-normalized leading coefficients in the frozen order
`(1,16/9,9/4,4,9)` are

```text
(-39.221008860498884,
  17.099036141737532,
   7.155067025676801,
  -4.846364092610365,
   0.5340006443586955).
```

Against the unchanged E4 ratios from `65b3830`, the three primary residuals
are

```text
(44.13204229285056,
 34.101648142954026,
 22.095529524318422).
```

Their full covariance is retained in `multiscale_score.json`; no diagonal
approximation, E4 amplitude fit, exponent fit or curve refit was introduced.
`rho=9` remains diagnostic and does not enter the decision.

## Interpretation boundary

This result says the fresh data do not reject the frozen conditional E4 shape
under the frozen scale model.  It does not prove the Virasoro-transparent
projector assumption, identify a field, or turn the pilot residual into a new
curve.  The clean scientific update is narrower and useful: the pilot's
near-threshold E4 tension did not reproduce at the preregistered precision.
