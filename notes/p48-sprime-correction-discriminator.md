# Prospective extension of the P48 `P4[S']` correction test

## Status

This note extends **already-frozen** P48 correction models to the new prospective `N=185,265` full-curve geometries. No coefficient is refit here.

The chronological source is now canonical on `main` in `scripts/score_p49_fullcurve_doubling.py` at commit `204a3eb42862a23e7a3520695be44f5fda97bd65`.

The fresh 100M `N=130,170` P49 counters reproduce the P48 anomaly: the pure `N^-5/4` derivative law fails strongly, while both frozen correction models survive. Those fresh counters are replication evidence only and are **not** used to update the coefficients below.

## Observable

For the derivative channel define

\[
y_N=N^{5/4}P_4[S'].
\]

The original parity model predicts a constant `y_N`. P48 showed a systematic drift, motivating two pre-existing competitors.

## Frozen model 1: q=2 ordinary correction

The first surviving alternative in the preregistered P49 chronology is

\[
P_4[S'] = A N^{-5/4}+B N^{-9/4},
\]

equivalently

\[
y_N=A+B/N.
\]

The previously frozen parameters are

```text
A = 3.203310807356976
B = -90.59560328584558
```

Prospective extension:

```text
N=185: y=2.7136048436497027 +/- 0.19668994042701207
N=265: y=2.8614406062783133 +/- 0.23233037114485428
```

These uncertainties are source-parameter uncertainty only.

## Frozen model 2: rank-2 Jordan/log correction

The second frozen alternative is

\[
P_4[S']=N^{-5/4}\bigl(A+B\log N\bigr),
\]

with previously frozen parameters

```text
A = -2.422594685734799
B =  1.016646899281392
```

Prospective extension:

```text
N=185: y=2.8846638769766324 +/- 0.2344225574126884
N=265: y=3.2500203406819943 +/- 0.32328609208512216
```

A pass would justify a dedicated logarithmic/Jordan study; it would not identify a specific LCFT partner.

## Scoring order

Preserve the existing chronology:

1. original pure `N^-5/4` law;
2. q=2 ordinary correction;
3. rank-2 Jordan/log correction;
4. zero effect;
5. only then free exponents or extra terms.

Target sampling covariance must be added to the source-prediction covariance in `predictions/p48_sprime_correction_competitors_20260828.yaml`.

## Boundary

This is secondary to Issue #43's original central-amplitude predictions and to the norm-5 H4/H12 discriminator. It must not be used to retune those primary tests.
