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
| 0 | `+1` | **cancels** |
| 4 | `-1` | doubles |

This is the axis/diamond pairing the orientation programme already uses, applied to
a rectangle instead of a square.

**One caveat, stated up front, because it limits the realized design and not the
idealization.** Two Gaussian integers of the *same* norm can never differ by exactly
45°: multiplying by `1+i` doubles the norm. So no realizable same-`N` pair is the
idealized pair of the table above. What a realizable pair does cancel exactly is
**spin 0** — any orientation difference does that, since a scalar is
orientation-independent, and that is where all scalar finite-size corrections live.
Weight 8 is *not* cancelled by the realized pair; it leaks at a coefficient computed
exactly below.

The score is

```text
D(r) = Â_axis(r) - Â_diagonal(r),     reported as D(r)/D(1).
```

## What it separates

| candidate | weight | `D(2)/D(1)` |
|---|---:|---:|
| **E4 — the Q4-Jordan prediction** | 4 | **11/4 = 2.75** |
| E8 (= E4²) | 8 | — cancels for the idealized pair; see the leakage below |
| E6 | 6 | — vanishes at the square point |
| Δ | 12 | 0.125 |
| E4³ | 12 | 20.796875 |
| E12 | 12 | 32.515625 |
| spin-4 amplitude ∝ area, no modular structure | — | 4 |

Two things follow.

**Modular weight is mostly not the hard part.** The nearest surviving modular
competitor sits about a factor of 12 away, so no precision argument is needed there.
Weight 8 needs the caveat above and is quantified below.

**The only real competitor is plain area scaling**, and against it the entire
discriminating content is the single constant `1/E4(i) = 0.68692…`. The measurement
is 11/4 against 4 — a gap of 45% of the predicted value, needing the ratio to about
**15% relative** for a 3σ separation.

## One rectangle is the whole experiment

> **Superseded 2026-09-05, same day.** The reasoning below is correct *against the
> area competitor* and the arithmetic stands. The overreach is the heading: it treats
> the competitor list as complete. The N=290 run landed at `1.880 ± 0.177`, which is
> none of the listed hypotheses, and against a linear-in-`r` law longer rectangles
> separate a great deal (at `r=4`: 4.0 against 10.99 against 16). See
> `notes/modulus-fingerprint-n290-result-20260905.md`.

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

So the production ask *against the area competitor* is one aspect ratio, two
orientations, and enough statistics to reach ~15% on the ratio. That was the right
first run. It was not the whole experiment.

## The concrete run: N = 290, four period matrices

A square torus of `N` sites needs a Gaussian integer of norm `N`; a 1:2 rectangular
torus of the *same* site count needs one of norm `N/2`. Both families need two
representations to separate the spin-4 amplitude `A` from the scalar part `C` in
`O(θ) = C + A·cos 4θ`. Sweeping `N` up to 1200, 46 sizes qualify. The best of them
is one this repository has already run:

| family | modulus | lattice | `w` | period vectors | `cos 4θ` |
|---|---|---|---|---|---|
| square | `i` | `<w, i w>`, `|w|² = 290` | `17+i` | `(17,1),(-1,17)` | `+20447/21025` |
| square | `i` | | `13+11i` | `(13,11),(-11,13)` | `-19873/21025` |
| rectangular | `2i` | `<w, 2i w>`, `|w|² = 145` | `12+i` | `(12,1),(-2,24)` | `+19873/21025` |
| rectangular | `2i` | | `9+8i` | `(9,8),(-16,18)` | `-20447/21025` |

All four have exactly 290 sites. The two families sample **the same pair of `cos 4θ`
values with the roles exchanged**, so their angular leverage is the identical rational
`8064/4205 = 1.9177…` — within 4% of the maximum possible 2, and equal between
numerator and denominator of the score, so the ratio estimator pays no relative
variance penalty.

That is not a coincidence. Multiplication by `1+i` maps the norm-145 representations
to the norm-290 ones — `(12+i)(1+i) = 11+13i`, `(9+8i)(1+i) = 1+17i` — and that map
*is* the 45° turn, which flips the sign of `cos 4θ`. The rectangular family is the
square family rotated.

Estimator: within each family, two members determine `C` and `A` exactly, so
`A = (O₁ - O₂)/(cos 4θ₁ - cos 4θ₂)`, and the score is `A_rect / A_sq`, predicted
`11/4` against `4`.

### The spin-8 leakage, exactly

Because the realized pair is not a 45° turn, a spin-8 component of amplitude `A₈`
enters the projector at `(cos 8θ₁ - cos 8θ₂)/(cos 4θ₁ - cos 4θ₂)`. At N=290 that is

```text
square:      +1148/21025 = +0.054602
rectangular: -1148/21025 = -0.054602
```

— exactly equal and opposite, for the same `1+i` reason the leverages are equal. So a
spin-8 component biases the score by roughly `-0.055·(A₈/A₄ + A₈'/A₄')`: up to ~11%
if the two amplitudes were comparable, which is the same order as the precision we
are trying to reach. The committed H4-beats-H8 results (`H4 0.4163/2` against
`H8 16.0120/2`) bound `A₈/A₄` well below 1, so the real bias is expected to be a few
percent — but it is a **systematic on the score, not a statistical error**, and it
does not shrink with samples.

The fix is three orientations per family, which determine `C`, `A₄` and `A₈`
together. The smallest size where both families have three representations is
**N=650** (square `|w|²=650`: `25+5i`, `23+11i`, `19+17i`; rectangular `|w|²=325`:
`18+i`, `17+6i`, `15+10i`). That is 2.24× the sites of N=290 and needs three coupled
runs rather than two, so it is the upgrade to make only if N=290 comes back
interesting.

Two things to be honest about. Two members per family determine `C` and `A` with
nothing left over, so this size cannot check the `cos 4θ` form itself — that rests on
the orientation programme at other sizes. And the engine couples two period matrices
per run, so this is two runs; treating them as independent for the ratio is
conservative but gives up the shared-field variance reduction.

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
