# Project Status and Claim Ledger

**Status date:** 2026-08-28

`main` contains the governance/reference layer and the Huawei numerical research archive. Claim strength is determined by evidence, not by branch location. Claim levels are defined in `GOVERNANCE.md`; the execution-facing summary is `notes/SYNTHESIS-20260828.md`.

## Exact/background facts

| Statement | Level | Status |
|---|---:|---|
| Square-site `p_c` has no known closed form | background | Current project/literature position |
| Square-site and NN+NNN matching-site thresholds satisfy `p_c + p_c_hat = 1` | C5 | Exact structural constraint |
| Square-bond and triangular-site thresholds are `1/2` | C5 | Exact controls |
| A rounded decimal is not a definition of `p_c` | governance | Enforced in the literature provenance layer |

Published threshold estimates remain method-specific in `data/literature_threshold_sources.json`; incompatible quoted uncertainties are not collapsed into one synthetic interval.

## Strong current finite-size evidence

| Claim | Level | Evidence | Current interpretation |
|---|---:|---|---|
| Primitive same-`N` Gaussian tori have a nonzero orientation-dependent matching signal | C3 | P31 | Independent-seed confirmation at five frozen sizes |
| Tested signs agree with `Delta cos(4 theta)` | C3 | P31/P37/P50 | Strong evidence for an odd square-harmonic sector; not unique H4 identification |
| Current finite-size data are compatible with `DeltaM ~ DeltaCos4 N^-13/8` | C3 | P31/P32/P37/P50 | Held-out and three prospective Gaussian lineages support the law over the tested range |
| Gaussian `1+i` fixed-coordinate doubling follows `-2^-13/8` on three lineages | C3 | P37/P50 | Parameter-free sign/radial transformation passed prospectively |
| Local residual-to-root conversion satisfies `-DeltaRoot*mean(M')/DeltaM ~= 1` | C3 | P35/P45 | Finite root movement is locally explained by the measured residual and slope |
| Angular-normalized root amplitude passes the frozen N=65/85 test | C3 | P45 | Clean high-stat support for the root-moving sector |
| Clean 100M full-curve doubling on 65->130 and 85->170 is compatible with the frozen thermal-even residual/root mechanism | C3 | P49 / PR #73 | Residual/root joint scores are about `4.45/2`; the 85->170 lineage supplies the visible ~2.1 SE tension |
| The raw center-slope ratio is already asymptotically `2^(3/8)` at these sizes | C2 negative refinement | P49 | False at 100M precision: the ~0.2% drift is decisively resolved and requires finite-size correction |
| The pure `P4[S'] ~ N^-5/4` law is sufficient | C2 negative refinement | P48/P49 | Fresh replication rejects the pure law; prespecified q=2 and Jordan-log corrections both survive current same-geometry tests |
| Wrapping-only GLS helps variance | C1 negative | server archive | No: tested matching-difference wrapping channels are configuration-identical |
| Single-geometry motif controls clear a robust multi-size production gate | C2 mixed/negative | P34 | No; paired same-`N` controls remain the relevant variance-reduction route |

### Key numbers

P31 independent confirmation:

```text
N=65   z=16.03
N=85   z=11.23
N=130  z=5.22
N=145  z=5.27
N=170  z=2.58
A4 pooled = 0.7885 +/- 0.0352
```

Three prospective fixed-coordinate Gaussian-doubling lineages:

```text
frozen ratio = -2^(-13/8) = -0.3242098887...
65 -> 130 = -0.31382 +/- 0.0908
85 -> 170 = -0.34095 +/- 0.1118
145 -> 290 child residual z = -0.483
```

P45 angular-normalized root amplitude:

```text
A_p(65) = 0.42034 +/- 0.02157
A_p(85) = 0.39495 +/- 0.03078
frozen  = 0.45101 +/- 0.02013
```

Clean P49 full-curve replication:

```text
thermal-even doubling: chi2 = 4.448 / 2 at u=0
raw root residual:      chi2 = 4.481 / 2
finite-slope root:      chi2 = 4.448 / 2
slope ratios: 1.2939835, 1.2943776
asymptotic target: 1.2968396
raw slope chi2: 6412.89 / 2
```

P49 fresh `P4[S']` correction discrimination on N=130/170:

```text
pure N^-5/4:  chi2 = 37.887 / 2
q=2 correction: chi2 = 1.790 / 2
Jordan log:      chi2 = 0.677 / 2
```

These correction scores are fresh-seed replications on previously used geometries, not new-geometry model selection.

## Open interpretations

| Interpretation | Level | What would move it forward |
|---|---:|---|
| Leading matching residual is asymptotically `N^-13/8` | C2/C3 | Prospective N=185/265 score, norm-5 multiplier, additional full-curve lineages |
| Leading angular harmonic is uniquely H4 | C2 | Norm-5 H4-vs-H12 test (#57); N=185/265 supplies an earlier weaker H12 lever |
| Root gap follows a clean angular-normalized `N^-2` law | C2 | Prospective new-size/full-curve tests and the third doubling lineage; raw `-1/4` remains compatible but not decisive |
| Thermal metric/slope corrections are explained by the minimal scalar+H4 relative `1/N` structure | C0 | Frozen prospective 145->290 prediction on PR #83; richer models only after that score |
| `P4[S']` drift is ordinary q=2 correction versus LCFT/Jordan logarithm | C0/C2 | Fresh N=185/265 full curves with chronological frozen scoring (#48/#72) |
| `x=21/4`, spin-4 thermal-family LCFT field is the mechanism | C0/C2 | Radial competitors, harmonic discrimination, parity controls and derivative spectrum (#37/#43/#44/#48/#57) |
| Matching/complement defines an RG parity involution | C0/C2 | Self-matching/self-dual controls and derivative-parity tests; a full OPE/interchiral automorphism is a strictly stronger hypothesis |
| `V_<1,4>` supplies the historical post-`L^-7` scalar correction | C0 | Conditional theory in PR #63; exact multi-angle H0 test is gated behind #57/#43 (#74) |
| `kappa3=-5/3` is universal/exact | C0 | Same-modulus controls and continuum bridge (#25/#54) |
| A simple algebraic/PSLQ formula gives `p_c` | C0 | Bounded post-provenance search (#1); low priority |

## Current execution priorities

If compute or attention is scarce:

1. **#43 prospective N=185/265 full curves** — genuinely unused geometries; report the original frozen endpoints before all secondary competitors. Production is frozen at 500M paired permutations per target.
2. **#57 norm-5 Gaussian multiplier** — highest-value discriminator of H4 versus H12 and also a parameter-free radial competitor test.
3. **#50 / PR #83 third-lineage full curve (145->290)** — test the resolved slope finite-size correction and induced root prediction without changing `y_t=3/4`.
4. **#44 exact C4 self-matching parity control** — direct control of the matching-parity mechanism.

The N=1105 four-angle campaign remains gated. PR #77/#74 provides an exact H0/H4/H8/H12 projector design, but it should run only after the cheaper norm-5 test and the N=185/265 score.

The synthetic design red-team merged in #71 shows that the old five-size geometry is much better at angular discrimination than at selecting among radial correction mechanisms. More replicas on the same five points are therefore lower information value than new Gaussian multipliers, new sizes, or full-curve levers.

## Reproducibility and covariance scope

#39 remains a precision/reproducibility task rather than a global scientific blocker. P37/P45/P49/P50 clean runs provide stronger provenance than the historical pilots. Cross-size covariance diagnostics from #69 are on `main`; the old hidden coupling does not explain the observed radial tension.

## Engineering status

- Huawei result archive and production tools are integrated into `main` via PR #21.
- N=185/265 threshold-rank production designs are on `main` via PR #66.
- Clean P49 full-curve results are on `main` via PR #73.
- Synthetic model-discrimination red-team is on `main` via #71.
- Literature provenance layer is on `main` via PR #62.
- CI covers Python 3.9/3.11/3.13 plus C++17 build/self-tests.

## Explicit non-claims

The project does **not** currently claim:

- a closed form or exact new value for square-site `p_c`;
- proof that `13/8` is the unique asymptotic radial exponent;
- proof that H4 is the unique angular harmonic;
- an exact bare `2^(3/8)` center-slope ratio at the tested finite sizes;
- proof of the `x=21/4` LCFT operator identification;
- proof that matching parity extends to a full local OPE/interchiral automorphism;
- universality of `kappa3=-5/3`;
- a rigorous new percolation bound.
