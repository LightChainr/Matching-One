# Issue #225 norm-5 multi-radius pivotal production

## Production and integrity

The frozen 200,000-replica, 200-batch job completed on the 16-vCPU Huawei
ARM64 environment at commit `89dd3606ac14f1d7630321e6326c8660bbb93c8d`.
The engine used the fresh counter interval
`[15000000000,15000200000)`, the four frozen norm-5 geometries and the
Euclidean radii `2,4,7,8`.

Wall time was 4.02 seconds, aggregate CPU utilization was about 1401%, peak
resident memory was 6552 KiB, and the process exited zero.  SHA256 values for
every transferred raw file match the values computed on the remote host.

The local scorer was rerun after transfer and now records all three declared
covariance-aware structural tests in `analysis.json`.

## Same-R UV contrasts

Every row is the frozen first-minus-second orientation contrast.  Plus and
minus are correlated views of the same configurations.

| N | R | Delta A_plus | Delta A_minus |
|---:|---:|---:|---:|
| 325 | 2 | -0.00711 +/- 0.00737 | +0.00647 +/- 0.00740 |
| 325 | 4 | -0.00851 +/- 0.00890 | +0.00025 +/- 0.01015 |
| 325 | 8 | -0.03606 +/- 0.00906 | -0.01288 +/- 0.00871 |
| 425 | 2 | +0.00022 +/- 0.00773 | -0.00050 +/- 0.00750 |
| 425 | 4 | -0.00436 +/- 0.00953 | +0.03254 +/- 0.01039 |
| 425 | 8 | -0.03708 +/- 0.01046 | -0.01921 +/- 0.01074 |

These are UV views because fixed `R` means `delta=R/sqrt(N)` changes with
size.  They are not six independent detections.

## Dyadic shell result

For each size and channel the scorer forms

```text
[Delta A(2R)-Delta A(R)] / log(2)
```

inside every aligned delete-one replicate.  A full-covariance GLS fit asks
whether the four shell coordinates share one constant amplitude.

| channel | common amplitude | SE | chi-square / df | survival p |
|---|---:|---:|---:|---:|
| plus | -0.020953 | 0.006372 | 3.2375 / 3 | 0.356442 |
| minus | -0.011584 | 0.005978 | 13.4140 / 3 | 0.0038218 |

The plus channel is compatible with one common per-log shell amplitude.  This
is a consistency pass, not proof of a logarithmic field; its fitted amplitude
is nevertheless resolved from zero at about 3.3 source standard errors.

The minus channel decisively fails the common-shell model.  Its four shell
values change sign and direction strongly enough that a single logarithmic
generator is rejected at the current radii and sizes.  This is the opposite
of the simplest motivating picture in which the matching-odd local mark alone
would carry a stable mesoscopic Jordan flow.

## Matched-delta result

The separately frozen comparison uses

```text
N325, R7: delta=0.3882901374
N425, R8: delta=0.3880570001
```

whose relative cutoff mismatch is only `6.006e-4`.  The joint plus/minus
N325-minus-N425 contrast is

```text
(+0.0301546, +0.0190732)
```

with its full 2x2 covariance retained in `analysis.json`.  Equality at fixed
delta gives

```text
chi-square = 6.18262 / 2, p = 0.0454424.
```

This is a marginal failure of fixed-delta collapse, not a discovery threshold.
Taken together with the minus-shell rejection, it says the present local odd
observable still carries measurable size/UV/topology dependence after the
second cutoff is matched.

## Scientific decision

The production establishes that the two-cutoff experiment is statistically
usable.  It does not support the simplest hypothesis that `A_minus` is one
stable logarithmic shell direction.  The cleaner constant-shell behavior is
instead in `A_plus`, the ordinary/local control channel.

The next high-information step should not fit a free radius exponent to this
block.  It should either:

1. repeat one larger matched-delta pair to learn whether the p=0.045 boundary
   moves toward collapse; or
2. enlarge the local observable basis so a rotating two-component shell flow
   can be tested directly rather than forced into one scalar `A_minus`.

All shell, same-R and matched-delta summaries originate from one raw random
block and count as one correlated result.
