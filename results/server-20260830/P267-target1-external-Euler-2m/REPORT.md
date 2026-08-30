# P267 Target 1: external Euler/rank-birth 2M reveal

## Decision

The external observer passes its bulk/contact gate. The connected complex
`P4 Cov(O_ext,J_D4)` is nonzero by an enormous margin at both production
sizes, and almost all of it remains after removing the frozen radius-2 local
Euler nuisance. This is not the old q-only contact identity and it is not a
root-neighbourhood artifact.

The source identity is not yet resolved. JD and JS retain the same complex
transfer phase, so the production result establishes a bulk Euler-to-rank-
birth response but does not identify that response with thermal Q4 epsilon.

## Frozen production

The preregistration is commit `5ccab6f`; the runner source is `01d8a4e`.
Both ARM hosts independently produced the identical binary SHA-256
`30e27d1c9483c77dc63d434b320a1b76d6ae99fb099ce43b13e30188de645f67`.

| size | geometry pair | samples / batches | seed | counter interval | wall time |
|---|---|---:|---:|---:|---:|
| N325 | `(17,6)` / `(18,1)` | 2M / 100 | 20260830325 | `[11000000000,11002000000)` | 238.892 s |
| N425 | `(16,13)` / `(19,8)` | 2M / 100 | 20260830425 | `[13000000000,13002000000)` | 307.775 s |

N325's first launch only discovered that `/usr/bin/time` was absent and did
not start the runner. Its 72-byte failure log is retained. The frozen seed and
counter were then used exactly once by the successful job.

## Primary and contact decomposition

Complex component standard errors come from the frozen within-size delete-one
covariance. Mahalanobis values use the full real/imaginary 2x2 block.

| channel | N325 | N425 | chi2 N325 / N425 |
|---|---:|---:|---:|
| external D | `-22.11585(400)+22.04220(449)i` | `-27.12570(520)-18.20644(472)i` | 928486 / 914809 |
| far D | `-21.91038(472)+21.94033(501)i` | `-27.18853(612)-18.11880(552)i` | 639684 / 590521 |
| near D | `-0.20547(323)+0.10187(272)i` | `0.06283(289)-0.08765(309)i` | 91.51 / 9.43 |
| external S control | `2.20014(501)-2.26699(548)i` | `2.44088(525)+1.62521(634)i` | 7472.8 / 4174.6 |
| far S control | `0.74058(584)-0.71007(642)i` | `0.72527(578)+0.51878(706)i` | 652.0 / 331.6 |

The near-D amplitude falls from `0.2293` to `0.1078`; at N425 its 2D chi-square
is only `9.43`. In contrast the far-D amplitudes are `31.007` and `32.673`.
Thus the root-local R2 window neither creates nor explains the production
signal.

## Frozen transfer and null

| channel | `C425/C325` | amplitude | phase (rad) |
|---|---:|---:|---:|
| external D | `0.20370+1.02625i` | `1.04627 +/- 0.00160` | `1.37486 +/- 0.00245` |
| far D | `0.20613+1.03336i` | `1.05372 +/- 0.00195` | `1.37391 +/- 0.00280` |
| external S | `0.16893+0.91275i` | `0.92825 +/- 0.01799` | `1.38778 +/- 0.03104` |
| far S | `0.16031+0.85421i` | `0.86912 +/- 0.05879` | `1.38528 +/- 0.11428` |

The preregistered D/S transfer-phase contrast is

```text
arg[(T_D)/(T_S)] = -0.01293 +/- 0.03120 rad, z=-0.414.
```

The contact-phase-lock null is therefore not rejected. JD is nearly
antiparallel to JS at each size: `JD/JS=-9.883-0.164i` at N325 and
`-11.141-0.041i` at N425. The far-only ratios are likewise nearly real and
negative. This is a one-phase source-plane signal, even though its D and S
amplitudes are plainly different.

## Exact and provenance gates

- all 400 production complement-audit rows have zero endpoint, site, line,
  local-mark and saturation-index failures;
- every path row passed the exact active/inactive Horvitz S/D identities;
- Russo `B-A'` is exactly zero at N325 and at most `4.28e-50` at N425;
- the compiled self-test passed arbitrary periods, basis invariance, physical
  `P ell` spin-4, saturation gcd, direct `0->2`, Euler and Gram gates;
- the main-integrated `83e98fc` alias gate is respected: Euler is scalar and
  the H4 leg remains typed complex internal data;
- every locally retained aggregate matches the remote SHA file. The 484/497MB
  sparse tables remain identified by SHA
  `427a0bc1...6101a9a79` and `28c750fc...bb864c0`; their row-wise Python replay
  was interrupted when the released DevEnvs returned to Ready. No sparse
  replay result is claimed by the local score JSON's empty placeholder block.

## Lifecycle tuple

```text
state       = full pre-insertion configuration, not q alone
source      = typed rank-birth chi4(P ell) JD4, with JS4 Gram/control plane
observer    = Euler O_ext, split exactly into R2 O_near and O_far
geometry    = norm-5 Gaussian pairs N325 and N425
acquisition = common-permutation marked Newman-Ziff, 2M/100 per size,
              independent seed/counter groups across sizes
lifecycle   = exact gate -> 100k pilot -> 5ccab6f prereg -> 2M reveal ->
              bulk passed / source phase unresolved
```

## Next mechanism decision

Do not buy another larger-N repeat of the same scalar row. The recorded Gram
roots already define the next zero-cost test: construct, orientation by
orientation and inside every delete-one batch,

```text
J_D_perp = J_D - beta J_S,
beta = <J_D, J_S> / <|J_S|^2>,
```

then score `Cov(O_far,J_D_perp)` and its N325-to-N425 complex transfer. If that
orthogonal residual remains nonzero and coherent, the Euler bridge has a
source component independent of the projective JS lane. If it collapses, the
present result is a bulk response carried by one shared source direction, not
evidence for Q4 epsilon. Only after this stored-data projection should a new
charged-seam or macroscopically separated observer be acquired. qJ remains a
contact control and cannot rescue either outcome.
