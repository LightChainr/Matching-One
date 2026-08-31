# Common-next-label Euler-invisible tangent: complete paired thermal response

## Outcome: dominant same-parity channels, with a cancellation-hidden E shape

The common plus-mark perturbation has a resolved orientation-average response,
while its paired H4-difference response remains weak. The common minus-mark
perturbation strongly excites the paired A contrast at both sizes. At N425 its
paired E contrast exchanges sign across the thermal window and nearly cancels
after integration. These are distinct entries of the same finite
source/observer response matrix.

| Channel / descriptive principal shape | N325 | N425 |
| --- | ---: | ---: |
| plus->D A largest lobe: amplitude /pointwise SE | +8.1321e-5 /4.7946e-5 (1.70SE) | -6.9189e-5 /4.2431e-5 (1.63SE) |
| plus->D E largest lobe: amplitude /pointwise SE | +7.1519e-5 /3.5387e-5 (2.02SE) | -4.5647e-5 /3.1017e-5 (1.47SE) |
| minus->D A positive peak p | .63312566 | .62934444 |
| minus->D A positive peak | +3.05999e-4 +/-7.17445e-5 | +3.14221e-4 +/-4.54777e-5 |
| minus->D A integral | +2.73217e-5 +/-6.95976e-6 | +2.48726e-5 +/-3.49218e-6 |

No principal zero crossing is assigned to plus->D: at N325 the positive E
feature has only2.02 pointwise SE and its preceding negative feature has0.47;
at N425 even the main extrema stay below1.7SE. The source's immediate-law
invariance does not by itself establish a nonzero common-plus H4 tangent.
Small mathematical tail roots are retained as numerical descriptors only.

### N425 minus->D E: birth-component exchange with almost zero integral

| Feature | p | H1 | H2 | E=H2-H1 |
| --- | ---: | ---: | ---: | ---: |
| Earlier negative E extremum | .60628578 | +1.49270e-4 | +8.60802e-5 | -6.31893e-5 +/-2.09007e-5 |
| Later positive E extremum | .66784413 | +5.58103e-5 | +1.03467e-4 | +4.76565e-5 +/-2.02398e-5 |

The two normalized paired birth responses stay positive at these extrema;
their relative dominance switches. H1 is larger early and H2 later, while
their sum A retains its strong positive main lobe. These are signed paired
contrasts, so their signs are not a claim that both individual orientation
clocks accelerate or delay.

The numerical E crossing is `.63792755` with local-delta SE`.00871020`.
Its adjacent extrema have3.02 and2.35 pointwise SE, respectively. The negative
and positive lobe areas are
`-3.08404e-6 +/-1.20128e-6` and
`+2.60091e-6 +/-1.27919e-6`; the small late negative tail adds about-4.28e-8.
The full integral is only `-5.25797e-7 +/-2.21120e-6`. Thus the near-zero
integral conceals opposing thermal contributions rather than showing an
absence of finite-p response. All areas and the integral share the saved
twenty-batch covariance; they are not separate confirmations.

N325 has a compatible-looking later positive E feature (2.25 pointwise SE),
but its early negative feature is only1.34SE. Its numerical crossing near.609
is consequently **not** promoted to a corresponding two-lobe mechanism or a
size-dependent root law.

### Common plus->S: a distinct even response

The plus->S A response is negative across its principal thermal lobe, with
integrals `-1.06188e-5 +/-1.69614e-6` and
`-1.19278e-5 +/-1.90375e-6`. Its E response has an early positive and later
negative lobe at both sizes:

| N | Earlier E peak: p, amplitude +/-SE | Later E peak: p, amplitude +/-SE | Numerical exchange p +/-local SE |
| --- | --- | --- | --- |
| 325 | .59610, +1.90772e-5 +/-6.57845e-6 | .66440, -3.38614e-5 +/-1.23497e-5 | .62409 +/-.00979 |
| 425 | .59976, +2.06325e-5 +/-7.66037e-6 | .66169, -3.74039e-5 +/-1.40264e-5 | .62491 +/-.00953 |

At both extrema H1,H2 are negative; early first-birth depletion dominates,
then later completion depletion dominates. This is the common-label
orientation-average response, not the earlier orientation-specific policy
or the new plus->D H4 response. The two root point estimates being close
does not establish a universal shared root.

The cross-orientation antisymmetric reciprocal diagnostic is unresolved at
p_ref and in the full integral at both sizes. The score preserves it and
minus->S rather than imposing an exact reciprocity/selection rule. The
dominant plus->S and minus->D channels therefore suggest a useful finite
parity organization, not exact diagonal closure or continuum field identity.

## Source and linear maps

The source is the **same common next-label perturbation applied to both
orientations**, unlike the earlier equal mixture of two orientation-specific
perturbations at `7c60b8a7`. It keeps the jointly safe degree-pair class masses
fixed. Within each class it tilts by `exp(t*pi_class*g_plus/minus)`, where the
marks are half the sum/difference of the two R0-only loop counts. Immediate
joint rank and Euler-increment law is preserved by construction.

This readout consumes the already extracted complete signed integer birth
histograms at `4db356e1b026853468f94d59d938895a2367ceb7`, not the old histogram,
raw suffix files, or any new simulation. Source axes are
`batch,cell,mark,orientation,birth,k`. The five source cells are00,01,02,10,20;
the full original20k denominator remains in force. The producer's64000
divisor already includes both half-mark conventions.

## Full-p and integral maps

For each batch and mark, first accumulate integer birth coefficients

```
C_j(n)=sum_(k<=n) h_j(k),
H_j(p)=sum_n C_j(n) binomial(N,n) p^n (1-p)^(N-n) /64000.
```

Sum the five cell contributions, then form `S=(first+second)/2` and
`D=(first-second)/delta_cos4`, using the original normalization stored in the
source. The tangent constants vanish because every signed birth histogram
has zero total mass. Consequently, for each of the four mark/output channels,

```
H_A=H_1+H_2,             H_E=H_2-H_1,
integral H_j=-d<E K_j>/(N+1),
integral H_A=-(d<E K1>+d<E K2>)/(N+1),
integral H_E=(d<E K1>-d<E K2>)/(N+1).
```

The signed Kj first moments and all curve integrals are the same histogram
information in different coordinates, not independent confirmation.

## What the reciprocal diagnostic actually means

Write T_ij for the response in orientation j to the unhalved R0-loop mark
of orientation i, under the same common class policy. Algebra gives

```
(delta_cos4/2)*(g_plus -> D) - (g_minus -> S) = (T21-T12)/2,
(delta_cos4/2)*(g_plus -> D) + (g_minus -> S) = (T11-T22)/2.
```

Both diagnostics are retained with the original twenty-batch covariance.
There is no assumed Onsager law setting T12=T21, nor a direct comparison of
unadjusted S and D amplitudes with incompatible normalizations.

## Numerical scope

The fixed full grid has1001 equally spaced p values plus p_ref. Integer
cumulatives are formed before floating normalization; exact zero endpoint
factors are removed with the existing `7c60b8a7` routine before sign-bracket
root readout. No old "main root must exist" assumption is carried over.

All A/E sign lobes and their numerical extrema are retained. A crossing with
either adjacent extremum below2 pointwise batch SE is explicitly labeled
weak and receives no physical-phase interpretation. This is a descriptive
annotation after the curve readout, not a simultaneous significance gate.
Peak/root positions are not selection-adjusted confidence statements.

`thermal_curve.csv` holds every channel's F1,F2,A,E and pointwise SE. The
immutable signed histograms reconstruct the full covariance across any p
values; the score also retains joint twenty-batch p_ref/integral vectors,
their covariance, and the per-lobe batch amplitudes/areas. No matrix inversion,
new fit, raw-hist extraction or multipanel plotting is undertaken.

Reproduce this readout with the existing single-thread research Python:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
/Users/lc/python-envs/research-py311/bin/python scripts/p334_paired_euler_thermal_response.py
```

The source commit and histogram hash are recorded in `score.json`; no
histograms from the old orientation-specific calculation were reused as data.
