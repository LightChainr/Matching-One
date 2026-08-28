# Research Synthesis — 2026-08-28

This is the execution-facing synthesis of Matching One. It is intentionally shorter and more opinionated than the claim ledger: the purpose is to identify the strongest evidence, the weakest link, and the next experiments with the highest information value.

## Working picture

The project is no longer primarily a search for a closed-form decimal expression for square-site `p_c`.

The most coherent current picture is that the finite-size square-site/matching observable contains a lattice-orientation correction whose leading measured component behaves approximately as

```text
Delta M_N(theta_1, theta_2)
    ~ A * Delta cos(4 theta) * N^(-13/8),
```

and that Gaussian-integer multiplication acts on this correction in a nearly parameter-free way by changing both area and microscopic orientation.

This is a finite-size numerical mechanism, not yet an operator proof.

## Strongest evidence

### 1. Independent five-size same-N orientation confirmation

P31 used a fresh independent seed with 100 million paired replicas at each of

```text
N = 65, 85, 130, 145, 170.
```

Every orientation difference had the sign predicted by `Delta cos(4 theta)`:

```text
N=65   z=16.03   A4=0.8093
N=85   z=11.23   A4=0.8666
N=130  z= 5.22   A4=0.9330
N=145  z= 5.27   A4=0.7501
N=170  z= 2.58   A4=0.6277
```

Pooling the available seeds gives

```text
A4 = 0.7885 +/- 0.0352,
chi-square = 1.53 / 4.
```

The important fact is not the exact chi-square. It is that five frozen geometries and an independent seed reproduce one signed orientation sector with compatible scaled amplitudes.

### 2. Parameter-free Gaussian `1+i` doubling on two lineages

Multiplication by `1+i` maps `N -> 2N` and rotates the microscopic square lattice by `pi/4`. For an H4 correction with radial exponent `13/8`, the frozen prediction is

```text
Delta M(2N) / Delta M(N) = -2^(-13/8)
                         = -0.3242098886627524.
```

Fresh P37 data gave

```text
65 -> 130: -0.31382 +/- 0.0908
85 -> 170: -0.34095 +/- 0.1118
```

with both sign reversals correct. No amplitude or exponent was fitted to those four target points.

This is currently more informative than another free-exponent fit because it couples radial scaling and angular transformation in one test.

### 3. A third prospective Gaussian-doubling lineage passes

P50 froze the lineage-specific `145 -> 290` prediction before the target run. The fresh child result was

```text
Delta M_290 = -0.000160648 +/- 0.000040542
frozen target = -0.0001376564 +/- 0.000024997
residual z = -0.483.
```

The zero benchmark is about four sampling standard errors away. This gives a third exact Gaussian lineage with the predicted sign and compatible magnitude.

Three lineages make accidental geometry-specific agreement less plausible than the original two-lineage result.

### 4. The residual-to-root conversion works locally

Threshold-rank full curves show

```text
C = -DeltaRoot * mean(M') / DeltaM
```

within roughly `2e-4` of one over the tested pilot sizes. A clean high-stat P45 run further gives

```text
A_p = -N^2 DeltaRoot / DeltaCos4
N=65: 0.42034 +/- 0.02157
N=85: 0.39495 +/- 0.03078
frozen target: 0.45101 +/- 0.02013.
```

Direct roots and the independently reconstructed `A_M/B` agree closely. The orientation residual therefore really does move the finite matching root through the expected local linear mechanism.

## What is still weak

### H4 is not uniquely identified

The `1+i` sign reversal separates odd square harmonics from even ones, but

```text
H4, H12, H20, ...
```

all reverse under a `pi/4` rotation. The current evidence therefore supports an odd square-harmonic sector much more strongly than it supports the unique statement “the signal is pure H4.”

This is the single clearest empirical weakness.

### The `x=21/4` LCFT identification remains a theory hypothesis

The Virasoro calculation supplies an allowed spin-4 thermal-family candidate with the right scaling dimension, but it does not prove that the lattice observable couples uniquely to that field. Matching parity, logarithmic mixing, and competing fields remain unresolved.

### Full-curve Gaussian semigroup closure is incomplete

Fixed-`p` doubling works on three lineages. The stronger test is the full triptych

```text
DeltaM_(2N) / DeltaM_N       = -2^(-13/8)
mean(M')_(2N) / mean(M')_N   =  2^(3/8)
DeltaRoot_(2N) / DeltaRoot_N = -1/4.
```

The root-level ratio is especially useful because it links the observed angular sector directly to the classical `L^-4` finite-size root drift.

## If we can run only three next experiments

### 1. Full-curve Gaussian doubling triptych — #49/#50

Highest priority.

Use threshold-rank data on the three exact lineages

```text
65 -> 130
85 -> 170
145 -> 290
```

and score the residual, slope, and root ratios together. Do not spend effort making the covariance machinery philosophically perfect first; retain aligned batches and report both full-covariance and simple diagnostics.

A clean pass would turn the current fixed-coordinate observation into a much more coherent finite-size mechanism.

### 2. Norm-5 H4 versus H12 discriminator — #57

This is the most decisive harmonic test per unit compute.

The frozen norm-5 Gaussian multiplier predicts opposite signs for pure H4 and pure H12 on the chosen lineages. A sign-resolved result is more valuable than another moderate-N exponent fit.

If H4 wins, the project can substantially narrow the leading angular sector. If H12 wins, that is equally useful because it falsifies the current preferred harmonic assignment without destroying the broader odd-harmonic semigroup picture.

### 3. One exact parity control — #44 first, #42/#48 next

Run the C4 self-matching site triangulation control if implementation cost remains low. It can separate ordinary square-lattice anisotropy from the matching-odd central residual with an exact self-matching center.

This has higher theoretical discrimination value than adding many more square-site target sizes.

## Useful parallel work

These are worthwhile but should not displace the three tests above:

- prospective unused `N=185,265` score (#43);
- paired same-N motif variance reduction (#40);
- Yang-Zhou/Jacobsen-Reply primary-source table transcription (#4);
- finite-width annihilator correction-spectrum study (#47/#58);
- kappa3 continuum bridge (#54).

## What to deprioritize

For now, do not spend major effort on:

- broad PSLQ searches;
- N=1105 four-angle production before the cheaper norm-5 harmonic test;
- GPU work whose target statistic is already cheap enough on CPU;
- adding another theory note without a new falsifiable prediction;
- squeezing extra digits from the threshold before the finite-size mechanism is clearer.

## Current project thesis

A concise thesis that matches the evidence is:

> Square-site/matching finite-size corrections contain a reproducible orientation-dependent odd square-harmonic sector. Its measured radial behavior is consistent with `N^-13/8`, and Gaussian-integer multiplication produces the predicted sign/radial transformation on three prospective lineages. The next task is to determine whether the leading harmonic is genuinely H4 and whether the same semigroup law closes at the full-curve root level.

That is the line to optimize around until a decisive experiment breaks it.
