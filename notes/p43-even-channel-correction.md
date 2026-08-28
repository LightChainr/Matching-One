# P43 even-sector channel correction

## Executive correction

The prospective Issue #43 `DeltaS` conjunction failed as frozen, but the failure must **not** be interpreted as a sign reversal of the physical cross-wrapping matching-even sector.

The preregistration accidentally froze the P31 `either/even` amplitude while the P48/P49/P43 threshold-rank engine reconstructs rank-2 **cross** wrapping only.

The matching-odd difference is configuration-identical across the supported wrapping channels, so this bookkeeping error does not affect the Issue #43 `DeltaM` result. The matching-even sum is channel-dependent, so it does affect `DeltaS`.

## 1. The two source observables are not the same

The P31 fixed-p analyzer defines, channel by channel,

\[
S_x=\frac{R_G^x+\widehat R^x}{2},\qquad
D_x=\frac{R_G^x-\widehat R^x}{2}.
\]

Issue #43 froze the common-amplitude fit from

```text
channel = either
sector  = even
```

with

\[
A_I^{\rm either}=+0.010603216462677735\pm0.000936687018246324.
\]

The later threshold-rank engine does not retain the `either` event. Its ranks are explicitly

```text
K_plus  = first black occupation count with a rank-2 primal component
K_minus = first black count where complementary white rank-2 matching wrapping is lost
```

so its reconstructed `S` is the **cross/even** observable.

## 2. The P31 source already contains the correct cross/even sign

The same committed P31 100M table gives the scaled cross/even amplitudes

```text
N=65   -0.01025934158 +/- 0.00144330999
N=85   -0.01300564285 +/- 0.00180075293
N=130  -0.01244966667 +/- 0.00307983055
N=145  -0.00718869102 +/- 0.00239125942
N=170  -0.00816274479 +/- 0.00375246345
```

Their inverse-variance common amplitude is

\[
\boxed{
A_I^{\rm cross}=-0.010603216462677733
\pm0.000936687018246330,
}
\]

with

\[
\chi^2=4.6580/4.
\]

Within the committed table, `cross/even` and `either/even` have equal magnitude and opposite sign at every listed size. This is why selecting the wrong channel preserved an apparently excellent source fit while flipping the prospective sign.

## 3. Post-hoc methodological closure on N=185,265

The original Issue #43 frozen means were

```text
N185 +6.7521637459e-5
N265 +6.8919446970e-5
```

for `either/even`.

The already-existing P31 cross/even source instead predicts the same magnitudes with negative sign:

```text
N185 -6.7521637459e-5
N265 -6.8919446970e-5.
```

The actual threshold-rank cross/even observations are

```text
N185 -6.0815376233e-5
N265 -7.0249507845e-5.
```

Using the exact covariance already committed by the Issue #43 primary scorer, the post-hoc corrected residual score is

\[
\boxed{\chi^2=0.5700/2},
\]

with marginal residuals about `+0.67 sigma` and `-0.12 sigma`.

This is excellent physical closure of the cross/even `N^-1` law, but it is **not** a preregistered pass because the wrapping channel was corrected after the target was visible.

## 4. Correct evidence ledger

The repository should retain both statements:

1. **Governance/statistical:** the original Issue #43 two-sector conjunction failed exactly as frozen because its `DeltaS` artifact used the wrong wrapping channel. It must not be rewritten as a preregistered success.
2. **Scientific/mechanistic:** the target cross/even data are highly consistent with the pre-existing P31 cross/even `N^-1` amplitude. Therefore the failed frozen score is not evidence against the matching-even dimension-4 sector and should not be described as an observed physical sign reversal.

The matching-odd `DeltaM` score is unaffected because the matching difference is configuration-identical across the wrapping channels used here.

## 5. Required future practice

Every prediction artifact involving `S`, `R_G`, or `R_hat` must include an explicit wrapping-event field such as

```text
wrapping_channel: cross | either | both | direction_0 | direction_1
```

and every scorer must fail closed if the target engine reconstructs a different channel.

For matching differences `M=R_G-R_hat`, a channel-invariance certificate may be referenced instead, but the channel must still be recorded for provenance.

The next genuinely prospective even-sector test should use `cross/even` and the intrinsic/full-curve threshold-rank definition from the start rather than translating an `either` fixed-p coefficient after the fact.
