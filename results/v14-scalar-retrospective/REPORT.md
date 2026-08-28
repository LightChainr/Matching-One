# V_<1,4> scalar retrospective audit

Status: retrospective design/power analysis only. No claim upgrade.

## Why this was tested

The corrected critical-Potts convention

```text
h_(r,s)=((2r-3s)^2-1)/24
```

puts the percolation thermal field at `V_<1,2>` with `h=5/8` and gives

```text
V_<1,4>: h=hbar=33/8, x=33/4, spin 0.
```

If the conditional matching/interchiral parity recursion makes `V_<1,4>` odd,
then its central matching contribution is

```text
M_L(pc) ~ L^-25/4 = N^-25/8,
```

which is relative `q=3` to the leading `L^-13/4` H4 term and produces an
`L^-7` root correction after the leading term is annihilated.

## Retrospective fixed-p two-angle H4-null projector

For each existing P31 same-N pair, define

```text
M0 = [c1 M2 - c2 M1] / [c1-c2],
c_i = cos(4 theta_i).
```

This cancels a *pure H4* contribution. Using the independent 100M confirmation
seed at `p_ref=0.592746050790` gives:

| N | M0 | SE | z | N^(25/8) M0 |
|---:|---:|---:|---:|---:|
| 65 | +1.4893e-5 | 7.4674e-5 | +0.20 | +6.89 |
| 85 | -3.6192e-5 | 6.3340e-5 | -0.57 | -38.73 |
| 130 | +2.8757e-5 | 6.4806e-5 | +0.44 | +116.09 |
| 145 | -1.42397e-4 | 6.4450e-5 | -2.21 | -808.68 |
| 170 | -4.7068e-5 | 7.3934e-5 | -0.64 | -439.42 |

Four of five sizes are below `0.7 sigma`; N=145 is a lone `-2.21 sigma` point.
This is **not** evidence for or against V14.

The important design result is that the scaled error grows extremely quickly,
as expected for an `N^-25/8` signal. Moving the same two-angle strategy to
larger N is therefore a poor use of compute.

## Two-angle leakage audit

The H4-null weights do not isolate H0 in the presence of higher harmonics. The
response coefficients to H8 and H12 are:

| N | H8 leakage | H12 leakage |
|---:|---:|---:|
| 65 | -0.1484 | +0.6716 |
| 85 | +0.2225 | -0.7628 |
| 130 | -0.1484 | -0.6716 |
| 145 | +0.8384 | -0.1004 |
| 170 | +0.2225 | +0.7628 |

In particular the apparently largest N=145 scalar excursion occurs where the
H8 leakage is about `0.84`. It cannot be interpreted as an H0 observation.
This confirms the caveat already motivating the exact four-angle N=1105 design,
but it also shows why N=1105 should not be run merely to chase this rapidly
decaying scalar before cheaper tests succeed.

## Retrospective root point-value diagnostic

Using the clean P45/P49 full-curve point roots only (without recomputing their
joint jackknife here), the same H4-null root projector gives

```text
N=65   p_scalar=0.5927507191723271
N=85   p_scalar=0.5927529686036729
N=130  p_scalar=0.5927428381640026
N=170  p_scalar=0.5927465557129569
```

For the fixed V14 root exponent `beta_N=7/2`, the two no-external-pc doubling
reconstructions are

```text
65 -> 130: pc_hat = 0.5927420740345603
85 -> 170: pc_hat = 0.5927459339297528
```

They differ by about `3.86e-6`. Without the full propagated root covariance
this is a **point-value diagnostic only**, not a score. The new
`score_v14_scalar_root_projector.py` performs the correct delete-one-batch
calculation when run against the raw histograms.

## Main conclusion

The most efficient V14 test is not a larger-N angular scalar measurement.
A scalar H0 term cancels exactly in same-N orientation differences, so the
successful Gaussian orientation program is structurally blind to it at linear
order. The right observable is the historical **unprojected adjacent-size
Mertens-Ziff annihilator**, where the leading H4 amplitude is removed but an H0
`q=3` term survives.

The branch therefore adds:

- `src/threshold_rank_axis_mc.cpp`: missing axis-torus threshold-rank engine;
- `src/threshold_rank_axis_pair_mc.cpp`: exact nested-permutation CRN coupling
  for adjacent L and L-1;
- `scripts/score_axis_pair_annihilator.py`: direct fixed-p correction-shape and
  accelerated-root scorer;
- `experiments/v14_axis_annihilator_20260828.yaml`: pilot/heldout protocol.

The primary new statistic is

```text
F_L(p_ref)=L^(13/4)M_L(p_ref)-(L-1)^(13/4)M_(L-1)(p_ref),
```

with model

```text
F_L = C_q [L^-q-(L-1)^-q]
    + T   [L^4-(L-1)^4] + ...,
```

so `q={2,3,4,6}` can be compared while `T` absorbs the unknown microscopic
mistuning `p_ref-pc`. This is the recommended next V14 computation.
