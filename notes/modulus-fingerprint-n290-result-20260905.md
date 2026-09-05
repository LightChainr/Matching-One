# N=290 modulus fingerprint: the weight-4 shape is excluded, and so is everything else we named

**Date:** 2026-09-05  **Claim level:** C3 negative (prospective, frozen before the run)
**Frozen design:** `predictions/modulus_fingerprint_n290_v2_20260905.yaml`
**Artifact:** `results/modulus-fingerprint-n290/latest.json`
**Runner:** `scripts/score_modulus_fingerprint_n290.py`

## What was measured

The spin-4 amplitude `P4_S_prime` of the matching-odd slope, on a square torus
(`τ = i`) and on a 2:1 rectangular torus (`τ = 2i`) with the **same 290 sites**, each
family realized at two lattice orientations whose `cos 4θ` values are the same pair
with the roles exchanged, so both families divide by the identical rational leverage
`8064/4205`.

| family | τ | samples | `P4_S_prime` | |z| |
|---|---|---:|---|---:|
| square | `i` | 200,000,000 | `+0.002768 ± 0.000197` | 14.1 |
| rectangular | `2i` | 100,000,000 | `+0.005203 ± 0.000307` | 16.9 |

```text
ratio  =  1.880 ± 0.177
```

## The result

Every prediction on the frozen competitor list is excluded:

| hypothesis | predicted | z |
|---|---:|---:|
| **Q4-Jordan, area-normalized weight 4** | 2.75 | **−4.93** |
| no modulus dependence | 1 | +4.98 |
| plain area scaling | 4 | −12.01 |
| weight-12 Δ | 0.125 | +9.93 |
| weight-12 E4³ | 20.797 | −107 |
| weight-12 E12 | 32.516 | −173 |

The run was frozen, including its sample counts, seed, replica offset and stopping
rule, before any data at these period matrices existed. It is a prospective test and
it came back negative for the shape it was built to look for.

The other two spin-4 channels agree in direction and are consistent with the primary
one: `P4_S` gives `1.236 ± 0.304`, `P4_D_prime` gives `1.714 ± 0.414`. None is near
2.75 or 4.

## What is actually excluded

The hypothesis under test is a **conjunction**: that the amplitude carries the
area-normalized weight-4 shape `g₂(τ)`, *and* that the normalization removes the same
block so that the ratio is `Ê4(2i)/Ê4(i)`. Rejecting a conjunction does not say which
conjunct failed, and the second one was already flagged as an assumption — the
additive shape `A~(τ)` is not fixed by the Jordan relation
(`docs/astra/Q2-additive-shape-ambiguity.md`).

So this is a C3 negative for **the fingerprint as constructed**, not for the Q4-Jordan
module. What it removes is the cheap version of the test.

## The known systematic cannot explain it

The realized orientation pair leaks spin-8 into the projector at exactly
`±1148/21025 = ±0.0546`, equal and opposite between the families, so it biases the
score. Solving for the contamination that would move 2.75 to 1.880 requires

```text
A8/A4 = 3.44
```

in both families. The committed H4-beats-H8 results (`H4 0.4163/2` against
`H8 16.0120/2`) put `A8/A4` well below 1, and even `A8/A4 = 1` moves 2.75 only to
2.465, still 3.3σ from the measurement. The systematic is real but far too small.

## The thing to be careful about

`1.880 ± 0.177` is compatible with **the aspect ratio itself**, `r = 2`
(`z = −0.68`). That is a **post-hoc observation, not a result.** `r¹` appeared in my
first competitor sweep this morning and I dropped it when I reframed the design
around the difference channel, so it was not on the frozen list. It has no
prospective standing and must not be reported as one.

It does, however, correct something I wrote earlier today. The design note concluded
"one rectangle is the whole experiment", on the grounds that longer rectangles do not
separate weight 4 from area scaling — which remains true. The overreach was treating
that competitor list as complete. Against a linear-in-`r` law, longer rectangles
separate a great deal: at `r = 4` the three hypotheses predict **4, 10.99 and 16**.

## The follow-up, and an arithmetic obstruction

**Aspect ratio 3 is impossible.** A 3:1 rectangular torus has `N = 3|ω|²`, and 3 is
inert in `Z[i]`, so `3 | a²+b²` forces `9 | a²+b²` and `N` keeps an odd power of 3 —
never a sum of two squares. No 3:1 rectangle shares a site count with a square torus
on this lattice. The reachable ladder is `r ∈ {2, 4, 5, 8, 9, 10, …}`: those `r` that
are themselves sums of two squares.

**Aspect ratio 4 at N=3380 is the designed follow-up**, and it is better than N=290 in
two independent ways:

- square `|w|² = 3380`: `58+4i`, `52+26i`, `44+38i`;
  rectangular `|w|² = 845`: `29+2i`, `26+13i`, `22+19i` — **three** orientations in
  each family, which fit `C`, `A₄` and `A₈` together and remove the spin-8 systematic
  instead of bounding it;
- it separates linear-in-`r` (4.0) from weight 4 (10.99) from area (16) by margins no
  precision argument is needed for.

It is 11.7× the sites of N=290, so it is a real production block rather than an
afternoon. It should be frozen separately, with the linear-in-`r` law named as a
competitor *before* the run, since this run is where it came from.

## Not established

- identification or exclusion of the Q4 Jordan module;
- that the measured amplitude is the root-normalized log slope rather than a leading
  amplitude with the same symmetry;
- anything about shapes carrying the same E4 factor;
- the linear-in-`r` law, which is a post-hoc reading of one point.
