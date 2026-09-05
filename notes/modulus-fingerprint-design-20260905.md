# The rectangular-torus fingerprint: what it can decide, and how much of it to run

**Date:** 2026-09-05  **Claim level:** C0 (exact design input, no lattice claim)
**Artifact:** `results/modulus-shape-discrimination/latest.json`
**Script:** `scripts/modulus_shape_discrimination.py`

`docs/ROADMAP.md` item 2 asks for "the smallest scalar-cancelled two/modulus score
that separates ordinary analytic correction, generic log mixing and Q4-Jordan
shape." You cannot pick that score without knowing what the alternatives predict
for it, so this note computes that first. Everything here is arithmetic about
candidate shapes. Nothing is claimed about the lattice.

## The design

A lattice `L = ω₁Z + ω₂Z ⊂ C` has modulus `τ = ω₂/ω₁` and orientation `θ = arg ω₁`.
Area-normalizing a weight-`k` lattice amplitude leaves

```text
Â_k(L) = (|ω₁|/ω₁)^k · (Im τ)^{k/2} f(τ) = e^{-ikθ} · ĝ_k(τ).
```

The orientation enters only through `e^{-ikθ}`. That is where the spin lives, and it
is the whole design:

> Run each torus twice — once axis-aligned, once with the lattice turned 45° — and
> take the **difference**.

Period vectors `(1,0),(0,r)` and `(1,1),(-r,r)` give the same 1:r torus at `θ = 0`
and `θ = π/4`. A 45° turn multiplies a weight-`k` amplitude by `e^{-ikπ/4}`, so

| weight mod 8 | factor | fate in the difference |
|---|---|---|
| 0 | `+1` | **cancels exactly** |
| 4 | `-1` | doubles |

This is the axis/diamond pairing the orientation programme already uses, applied to
a rectangle instead of a square. It is free, and it is exact rather than statistical.

The score is

```text
D(r) = Â_axis(r) - Â_diagonal(r),     reported as D(r)/D(1).
```

## What it separates

| candidate | weight | `D(2)/D(1)` |
|---|---:|---:|
| **E4 — the Q4-Jordan prediction** | 4 | **11/4 = 2.75** |
| E8 (= E4²) | 8 | — cancels identically |
| E6 | 6 | — vanishes at the square point |
| Δ | 12 | 0.125 |
| E4³ | 12 | 20.796875 |
| E12 | 12 | 32.515625 |
| spin-4 amplitude ∝ area, no modular structure | — | 4 |

Two things follow.

**Modular weight is not the hard part.** Weight 8 and everything congruent to it
disappears from the channel for free, and the nearest surviving modular competitor
sits about a factor of 12 away. No precision argument is needed for either.

**The only real competitor is plain area scaling**, and against it the entire
discriminating content is the single constant `1/E4(i) = 0.68692…`. The measurement
is 11/4 against 4 — a gap of 45% of the predicted value, needing the ratio to about
**15% relative** for a 3σ separation.

## One rectangle is the whole experiment

`E4(ri) → 1` very fast, so for `r ≥ 2` the prediction is `0.687·r²` to high accuracy:

| r | `D(r)/D(1)` | `r²` | ratio |
|---:|---:|---:|---:|
| 2 | 2.75 | 4 | 0.6875 |
| 3 | 6.18233512793 | 9 | 0.686926125 |
| 4 | 10.9908008589 | 16 | 0.686925054 |

The ratio column is the only thing distinguishing the hypothesis from the area law,
and it is already constant to 0.08% at `r = 3`. **Aspect ratios 3 and 4 test nothing
that `r = 2` has not already tested.** Buying them would be exactly the kind of
same-purpose repetition the roadmap tells us to stop.

So the production ask is one aspect ratio, two orientations, and enough statistics to
reach ~15% on the ratio — not a ladder of geometries.

## What this does not settle

The 11/4 is **conditional, not a prediction**. The additive shape `A~(τ)` is not
fixed by the Jordan relation — it shifts under `q̃ → q̃ + αq` — so the ratio is clean
only if the normalization removes the same block. That is the open question in
`docs/astra/Q2-additive-shape-ambiguity.md` and it is unchanged by anything here.

What the difference channel does contribute to it: any **scalar** part of the unfixed
additive shape has weight ≡ 0 mod 8 and drops out of `D` exactly. A spin-4 part does
not. So the channel narrows Q2 without closing it, and the narrowing is exact.

Also not established: that any lattice observable equals one of these shapes; that
the measured amplitude is the root-normalized log slope rather than a leading
amplitude with the same symmetry; that the channel is free of lattice-level spin-4
contamination of non-modular origin.
