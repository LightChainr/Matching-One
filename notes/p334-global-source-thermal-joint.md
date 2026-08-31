# P334: first-birth debt and the first width-sensitive thermal moment

This is a zero-MC, zero-DP join of the **same** twenty original batches per N,
not an additional dataset or independent confirmation. It completes the
full-observable result in `notes/p334-full-global-conditional-clock.md`.

## Complete microscopic source, rather than F2 alone

Use the prefix-safe global gate already defined in `3edc785a`. On each accepted
R1 orientation, let the eventual second-birth source be D (an original H2 gate)
or G (collective). For its probability `pi_s` and source-resolved second-birth
kernel `F2_s`, the full observable is

```
A_s = F2_s - pi_s*(1-F1).
A_D + A_G + A_remainder = safe_global_A.
```

`A_remainder` retains all full-A terms outside accepted R1 orientations,
including R0-containing pairs, whole-pair solver fallbacks and R2 terms. It is
not a selectively discarded error term. D/G are operational path labels, not
identified continuum fields. The first-birth debt is nonnegative before the
orientation contrast; its H4 projection need not be.

The full integral D/G point signs are opposite to their F2-only signs in both
sizes. Their shared uncertainty does **not** resolve these signs:

| Full integral H4 source | N325 mean +/- SE | N425 mean +/- SE |
|---|---:|---:|
| original H2 direct A | -0.000062209 +/- 0.000696279 | 0.000241509 +/- 0.000649903 |
| collective A | 0.000144162 +/- 0.000509421 | -0.000440505 +/- 0.000491393 |
| remainder A | 0.000026719 +/- 0.000821785 | 0.001142030 +/- 0.000682870 |

The new mechanism accounting is nevertheless exact: the earlier positive-F2
source loading cannot be transferred to the full topology observer without
the same source's first-birth debt. All source, debt and total uncertainties
come from their shared batch vectors, including the large positive covariance
that can make the difference more precise than either term.

## J1 exposes a width direction that J0 integrates away

For the baseline full path, define `C=(K1+K2)/2`, `W=K2-K1` and
`D_N=(N+1)(N+2)`. The exact thermal identities are

```
J0 = integral A(p) dp = 1-2*C/(N+1)
J1 = integral p*A(p) dp = 1/2-(C*C+C+W*W/4)/D_N
J1_center = 1/2-(C*C+C)/D_N
J1_width = -W*W/(4*D_N).
```

These are baseline path readouts; this join does not invent a conditional
replacement for J1. The same-batch H4 contrasts are:

| Component | N325 mean +/- SE | N425 mean +/- SE |
|---|---:|---:|
| J1 center | 0.000063761 +/- 0.000455161 | 0.000448864 +/- 0.000511207 |
| J1 width | 0.000006624 +/- 0.000017950 | -0.000036429 +/- 0.000012149 |
| J1 total | 0.000070386 +/- 0.000450460 | 0.000412435 +/- 0.000510452 |

The N425 width direction is approximately three batch SE below zero, while
the total remains much less precise because the center component dominates
its uncertainty. This exploratory same-source readout gives a concrete
width-sensitive observable, without treating a marginal three-SE direction
as an independent discovery or inferring a scaling law from these two sizes.
Its next discriminant would be the same named W² loading on a new independent
production block; no such production is started here.

## Sources and covariance handoff

- Full-global observer and original covariance: `3edc785a`.
- Marked full-A sources, positive F2 and first-birth debt: `2dd865f0b26a4d5d43f52b300293016e6ffd19b8`.
- Six orientation-resolved J0/J1-center/J1-width columns: `e64febe4ff10ca9cfb2f094c1b8ee8f733177fe1`.
- Shared underlying full births: `9c495ab13e65f2bc93dc0849ee3b73f88724c4b1`;
  conditional clocks: `0d1e586dafbade5e7d1f9bfc598170d0c881e337`.

`results/p334-global-source-thermal-joint/score.json` preserves the original
orientation source and thermal columns, named H4 contrasts, full shared LOO
covariance and immutable input hashes. Its covariance rank is at most 19;
there is no covariance inversion or new omnibus test. Every batch add-back
and source alignment agrees within `9.7e-16`.

Reproduce only this thin committed-vector join with
`python3 scripts/p334_global_source_thermal_joint.py`.
