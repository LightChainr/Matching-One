# Marked-birth production pilot

## Run identity

- Huawei DevEnv: `DevEnvC_HZsCM6` (`033945d8bf8b47a7acf475c595169e07`), ARM64, 16 vCPU / 32 GiB.
- Prereveal/runner commit: `6899b119db5b16e9918db53abf5280d990eb6653`.
- Binary SHA256: `46d8a2690b9a3b1899b3fe61e9a2c16019cb39487d493998c477ca302eaa1223`.
- Counter stream: seed `202608290215`, replica offset `9100000000`, 20 aligned delete-one batches.
- Frozen sizes: q2 `N=65` (20,000), q2 child `N=130` (20,000), max-leverage P50 `N=145` (10,000).
- The run used one thread while the host's other cores were occupied; counter RNG makes this identical to a multithreaded replay.

Every raw and scored file is covered by `SHA256SUMS`.

## Exact and production audits

All three runs have zero endpoint, site, primitive-line, local-mark, and
saturation-index complement failures. Sparse-schema failures are zero. The
canonical Russo residual `B-A'` is zero or below `2.2e-50` at 50-digit scoring
precision. Direct `0->2` rows are present and retain `ell=null`, `iota=0`,
`S=2`, `D=0`.

The path stores both halves of the source. At every microcanonical size,

```text
S_full = (S_active + S_inactive)/2,
D_full = (D_active - D_inactive)/2.
```

Tiny exhaustive controls certify this equality coefficientwise. The same
controls separately certify physical `P*ell` rather than period-coordinate
angles and raw-winding gcd accumulation for `iota`.

## Pilot scores

Values are same-size `P4` projections; parentheses give delete-one standard
errors for real and imaginary parts.

| design | intrinsic p | P4 mean J_D4 | P4 Cov(A_top,J_D4) | P4 gamma_D4 |
|---|---:|---:|---:|---:|
| q2 N65 | 0.593120626 | 0.03869+0.02907i (0.01138,0.01550) | 2.84779-0.85665i (0.01707,0.01315) | 0.33904-0.10099i (0.00191,0.00165) |
| q2 N130 | 0.593023248 | -0.01256-0.00993i (0.01155,0.01721) | 3.60813+1.05851i (0.02856,0.02332) | 0.33188+0.09841i (0.00256,0.00211) |
| P50 N145 | 0.592585217 | 0.01656+0.02369i (0.02765,0.00833) | 3.77238+0.18889i (0.03332,0.00718) | 0.33455+0.01697i (0.00307,0.00059) |

## New mechanism split

The q2 child changes the registered `Delta cos(4 theta)` exactly from
`+1152/845` to `-1152/845`. Two statistically distinct transfers appear:

1. The **mean line-odd source** has
   `J_D4(N130)/J_D4(N65)=-0.33074-0.00823i`, magnitude `0.33085`. Its effective
   area exponent is `1.596`, close to the conditional `x=21/4` prediction
   `13/8=1.625` and factor `2^(-13/8)=0.32421`. The N130 mean is individually
   noisy, so this is a high-value preregistered lead, not a confirmation.
2. The **connected response normalized by birth mass** instead obeys
   `gamma(N130)/conj(gamma(N65))=0.97854-0.00120i`. Its magnitude exponent is
   only `0.031`, while the raw connected response grows by `1.2644`, near the
   thermal mass factor `2^(3/8)=1.2968`.

Thus the production stream separates a possible high-dimension signal in the
mean `J_D4` from a leading thermal/gate-overlap contribution in
`Cov(A_top,J_D4)`. The original normalized connected proxy is not a clean
`Q4 epsilon` projector by itself. The next innovative move is not simply more
replicas: orthogonalize `J_D4` against the gate-level thermal source using the
same stored common-field covariance, then repeat the q2 transfer on the
residual. The near-conjugate q2 phase supplies an unusually sharp geometric
control for that subtraction.

## Scientific card

1. MECHANISM SPACE: splits a candidate high-dimension mean source from a thermal-dominated connected response in one common random field.
2. NOT PROVED: the N65/N130 mean-source ratio does not yet identify `Q4 epsilon`; N130 mean significance is low and no third radial point exists.
3. OBSERVER-SECTOR-SOURCE-GEOMETRY: `A_top` | Alexander odd | lifted-line `J_D4` | exact q2 Gaussian parent/child plus P50 leverage.
4. DEPENDENCY GROUP: all pilot channels and their covariances reuse one counter-coupled stream per size; they form one discovery block.
5. UPWEIGHT OBSERVATION: q2 transfer of the thermally orthogonalized `J_D4` residual, with the frozen conjugation/phase law and `2^(-13/8)` magnitude target.

## Orthogonalization follow-up

The requested held-out batch-Gram score is in `ORTHOGONALIZATION.md`. It
changes the interpretation: N65 gives a weak, noisy estimator-level
`alpha=0.00626-0.02099i`, and applying it to N130 leaves the connected q2
growth unchanged (`1.26442 -> 1.26447`). More importantly, exact gate algebra
closes the connected response in terms of source means. Therefore the
connected `A_top` response is a contact channel, not an independent matrix
element. The mean-`J_D4` radial lead remains separate and alive.
