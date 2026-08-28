# Research Synthesis — 2026-08-28

This is the execution-facing synthesis of Matching One. It is intentionally more opinionated than the claim ledger: the purpose is to state the strongest evidence, the weakest links, and the next experiments with the highest information value.

## Working picture

The project is no longer primarily a search for a closed-form decimal expression for square-site `p_c`.

The most coherent current finite-size picture is:

```text
same-N orientation signal
    -> odd square-harmonic sector
    -> central matching residual approximately DeltaCos4 * N^-13/8
    -> local root shift through DeltaRoot ~= -DeltaM/M'
```

Gaussian-integer multiplication supplies exact transformations of area and microscopic orientation, making it a practical finite-size spectroscopy tool rather than just a convenient geometry generator.

This is a numerical/RG mechanism. It is not yet a unique H4 or LCFT-operator proof.

## Strongest evidence

### 1. Independent five-size same-N orientation confirmation

P31 used a fresh independent seed with 100 million paired replicas at each of

```text
N = 65, 85, 130, 145, 170.
```

Every orientation difference had the sign predicted by `Delta cos(4 theta)`. Pooling the available seeds gives

```text
A4 = N^(13/8) DeltaM/DeltaCos4 = 0.7885 +/- 0.0352,
chi-square = 1.53 / 4.
```

The important observation is the reproducible signed orientation sector across frozen geometries, not the exact fitted amplitude.

### 2. Three prospective parameter-free `1+i` Gaussian lineages

For an H4-like correction with radial exponent `13/8`, multiplication by `1+i` predicts

```text
DeltaM(2N)/DeltaM(N) = -2^(-13/8)
                     = -0.3242098886627524...
```

Fresh results:

```text
65 -> 130: -0.31382 +/- 0.0908
85 -> 170: -0.34095 +/- 0.1118
145 -> 290: prospective child residual z = -0.483
```

All three exact Gaussian genealogies preserve the predicted sign/radial transformation without fitting a target amplitude or exponent.

This strongly supports an odd square-harmonic semigroup sector, while still leaving H4/H12/H20 aliases open.

### 3. Root movement is locally tied to the measured residual

Threshold-rank curves satisfy

```text
C = -DeltaRoot * mean(M') / DeltaM ~= 1
```

on the tested finite systems. P45 independently confirms the angular-normalized root amplitude at N=65/85:

```text
A_p(65) = 0.42034 +/- 0.02157
A_p(85) = 0.39495 +/- 0.03078
frozen  = 0.45101 +/- 0.02013
```

So the orientation residual is not merely correlated with the pseudo-critical root: it moves the root through the expected local linear mechanism.

### 4. Clean 100M full-curve doubling partially closes the mechanism and exposes the next correction

PR #73 adds clean threshold-rank full curves at N=130/170 and scores 65->130 and 85->170.

For the frozen thermal-even residual and root tests:

```text
X_even(u=0)        chi2 = 4.448 / 2
raw root residual  chi2 = 4.481 / 2
finite-slope root  chi2 = 4.448 / 2
```

The visible strain is concentrated in the 85->170 lineage at about 2.1 sampling SE. This is compatible but no longer an effortless closure.

The raw slope multiplier, however, is decisively too simple:

```text
slope ratios = 1.2939835, 1.2943776
target       = 1.2968396 = 2^(3/8)
chi2         = 6412.89 / 2
```

The point-level discrepancy is only about 0.2%, but the 100M data resolve it sharply. This is evidence for a finite-size correction to the thermal metric, not evidence that `y_t=3/4` should be replaced.

### 5. The derivative-parity pattern survives, but `P4[S']` needs a correction

Fresh P49 replication confirms that normalized `P4[S]` follows its pure law very well, while `P4[D]` and `P4[D']` show modest strain. `P4[S']` decisively fails the pure

```text
P4[S'] ~ N^-5/4
```

law.

Prespecified secondary models on the fresh N=130/170 data give

```text
pure N^-5/4:     chi2 = 37.887 / 2
q=2 correction:  chi2 =  1.790 / 2
Jordan log:       chi2 =  0.677 / 2
```

These geometries were not new, so this is mechanism development/replication rather than prospective selection between q=2 and Jordan log.

## What is still weak

### H4 is not uniquely identified

The `1+i` sign reversal distinguishes odd from even square harmonics, but

```text
H4, H12, H20, ...
```

all reverse under a `pi/4` rotation. This remains the clearest empirical ambiguity.

The norm-5 multiplier in #57 is the cheapest exact experiment that separates H4 from H12. It also supplies different no-fit ratios for the main lower-dimensional radial competitors, so one campaign can do angular and radial spectroscopy simultaneously.

### `13/8` is strongly supported over the tested range but not uniquely asymptotic

P31/P32 plus three exact `1+i` lineages are hard to explain as a random fit accident, but the synthetic red-team in #71 shows that the old five-size design has weak power to distinguish competing radial mechanisms. More samples on the same five points do not solve that design problem efficiently.

This is why new N, new Gaussian multipliers and full curves now outrank another five-size amplitude fit.

### The `x=21/4` operator identification remains conditional

The thermal-family level-4 spin-4 Virasoro candidate exists and has the right dimension, but lattice coupling, matching parity, logarithmic mixing and lower competitors remain open.

A useful logical separation is:

1. a matching/complement involution on RG theory space can produce empirical even/odd tangent sectors;
2. promoting that involution to a full local OPE/interchiral automorphism is a stronger assumption.

Numerical support for the first must not be described as proof of the second.

### Full-curve semigroup closure now requires finite-size slope structure

The original triptych

```text
DeltaM_(2N)/DeltaM_N       = -2^-13/8
mean(M')_(2N)/mean(M')_N   =  2^3/8
DeltaRoot_(2N)/DeltaRoot_N = -1/4
```

is too naive at current slope precision. The residual/root mechanism still works approximately, but the slope requires a controlled finite-size correction.

PR #83 freezes a minimal scalar+H4 relative-`1/N` correction from the two clean source lineages and makes a no-target-refit prediction for the still-unseen 145->290 full-curve slope/root result. This is the right next test before adding richer slope models.

## If we can run only three next experiments

### 1. Prospective N=185/265 full curves — #43

Highest immediate priority because these are genuinely unused geometries.

Production is frozen at 500 million paired permutations per target. Score the original frozen `DeltaM/DeltaS` endpoints first, then the already-frozen radial/harmonic competitors in chronological order. Do not allow secondary P48 correction work to contaminate the #43 primary reveal.

These two N also have materially different H12 leverage, so they provide the first prospective H12 constraint before the more expensive norm-5 campaign.

### 2. Norm-5 Gaussian spectroscopy — #57

This remains the highest-value dedicated operator discriminator.

For the frozen lineages 65->325 and 85->425, pure H4 and H12 predict opposite child signs. The same run also distinguishes the main fixed H4 radial competitors through different parent/child ratios.

The synthetic red-team says this kind of new semigroup lever is much more informative than simply shrinking the SE on the old five sizes.

### 3. Third full-curve doubling lineage — #50 / PR #83

Run/score the 145->290 threshold-rank curve. The primary residual law remains first; then test the already-failed raw slope baseline and the frozen two-sector finite-size slope correction. The induced root prediction is part of the same test.

A successful third-lineage slope correction would turn the P49 slope 'failure' into a reproducible finite-size structure rather than a loose caveat.

## Next control after those three

Run the exact C4 self-matching parity control (#44). It has higher theory value than adding another ordinary square-site size because it separates generic C4 anisotropy from a matching-odd central residual at an exact self-matching center.

## Gated later work

### N=1105 scalar/harmonic decomposition

Do not run N=1105 before #57 and #43.

PR #77 / #74 provides an exact four-angle decomposition of H0/H4/H8/H12 at N=1105. The H0 weights have modest norm and exactly annihilate H4/H8/H12 within the declared truncation. This is a much cleaner eventual test of the conditional `V_<1,4>` scalar mechanism than a two-angle H4-nulling combination.

It does not cancel H16 and higher harmonics, so it is a finite harmonic projector, not an all-orders scalar theorem.

### P48 `S'` correction robustness

The q=2 and Jordan-log alternatives were already frozen on #48 before the N=185/265 targets. PR #76/#72 adds only a stricter training-only robustness split using N=65/85/130; it must be scored after the earlier freezes and must not reorder them.

## What to deprioritize

For now, do not spend major effort on:

- broad PSLQ searches;
- another large-Pell scan;
- simply multiplying replicas on the same old five N to choose a radial mechanism;
- N=1105 before the norm-5 discriminator;
- GPU work whose end-to-end information-per-wall-time advantage has not been measured;
- new operator notes that do not produce a new frozen prediction or exact control.

## Current project thesis

A concise thesis consistent with the evidence is:

> Square-site/matching finite-size corrections contain a reproducible orientation-dependent odd square-harmonic sector. Its central residual over the tested range is consistent with `DeltaCos4*N^-13/8`, and Gaussian-integer multiplication produces the prescribed sign/radial transformation on three prospective fixed-coordinate lineages. Full-curve data confirm the local residual-to-root mechanism while resolving finite-size corrections in the thermal slope and in the `S'` derivative channel. The next decisive work is new-geometry and new-multiplier spectroscopy, not another free exponent fit on the old sizes.

That is the line to optimize around until a prospective experiment breaks it.
