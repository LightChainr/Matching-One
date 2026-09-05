# The N=580 ladder, rescored without a denominator

**Date:** 2026-09-06
**Tickets:** #579 (method), #567 (the run), #575 (the covariance)
**Claim level:** C2 — a reanalysis of a committed run, no new samples
**Artifact:** `results/aspect-ladder-n580-projective/latest.json`

The frozen design scored one number, `A4(4i)/A4(i)`, and deliberately kept the
r=2 rung out of the decision because it carries a spin-8 leakage of the opposite
sign. Two corrections have already been applied to that scoring: the ratio
z-test became a Fieller contrast (#577), and the covariance it needs was
measured by deterministic replay (#575). This note applies the third and last
one — stop forming a ratio at all — and the result is not a refinement of the
previous verdict. It replaces it.

## The three-entry test

A model that predicts proportions predicts a **ray**. Testing it is the
covariance-weighted distance from the measured response to that ray,
`min_a (y − a v)ᵀ S⁺ (y − a v)`, and no coordinate is nominated as a
denominator. For two entries and one ray this is Fieller's z, squared — verified
to 1e-14 across six orders of magnitude in the predicted ratio, which is the
anchor test in `tests/test_projective_inference.py`. For three entries it does
what Fieller cannot: use every rung at once.

| competitor | σ | σ across the missing covariance | required \|A8/A4\| to reach r=2 |
|---|---:|---|---:|
| bare_aspect_ratio | 2.74 | [2.14, 4.77] | 7.7 |
| q4_jordan_weight4 | 7.00 | [5.93, 8.86] | 32.1 |
| plain_area_scaling | 7.13 | [6.10, 8.93] | 32.0 |
| no_modulus_dependence | 9.22 | [9.21, ∞] | **0.4** |
| weight8_E8 | 10.47 | [10.03, 10.97] | 182 |
| weight12_E4_cubed | 11.25 | [11.13, 11.38] | 783 |
| weight12_E12 | 11.25 | [11.13, 11.38] | 785 |
| weight12_delta | ∞ | — | 222 |

Seven of the eight verdicts hold across **every** positive-definite value of the
covariance entry the first scoring run threw away. Only `bare_aspect_ratio`
moves with it, so that entry is worth one deterministic replay and nothing more.

## The shape, stated without any denominator at all

The second divided difference on r = 1, 2, 4 is a **linear functional** of the
response,

```text
f[1,2,4] = (m(4) − 3 m(2) + 2 m(1)) / 6
```

so it needs no ratio, no matrix inverse, and nothing about it degrades when one
rung sits close to zero. On `r²` it returns exactly 1; on any line, exactly 0.

```text
measured   f[1,2,4] = −4.663e−04 ± 1.530e−04     z = −3.05
                                    across the admissible range: z ∈ [−3.79, −2.62]
```

Every frozen competitor predicts this to be **zero** (the two families linear in
r) or **strictly positive** (every modular family: +0.79, +1.00, +16.7, +211,
+331). The measurement is negative, and the sign does not flip anywhere in the
admissible range.

**The response is concave in the aspect ratio. Nothing in the frozen list is.**

## The dichotomy, and the assumption it lands on

The r=2 rung is contaminated by spin-8 leakage, so the natural objection is that
the concavity is that leakage. Fix each model's amplitude from the r=1 and r=4
rungs, whose leakages share a sign, then predict r=2 and read the gap as spin-8.
That is the last column above, and it produces a clean trap:

- **`no_modulus_dependence`** is the only competitor the r=2 rung can accept
  within the assumed bound — it needs \|A8/A4\| = 0.4. It is also the one the
  clean r=1 / r=4 pair excludes outright, at 9.2σ.
- **Every other competitor** needs \|A8/A4\| between 7.7 and 785.

So under the frozen design's own assumptions, no competitor fits all three
rungs. The escape hatch is a spin-8 amplitude **comparable to or larger than**
the spin-4 amplitude it is supposed to perturb.

That is exactly the assumption the design ruled out, and here is where it came
from. `predictions/modulus_fingerprint_n290_v2_20260905.yaml` says the committed
H4-beats-H8 results "bound `A8/A4` well below 1", quoting `H4 0.4163/2` against
`H8 16.0120/2`. Those numbers are a model-selection result on a **homology
character** — χ²/df of 0.21 against 8.0. They are not a measurement of the
**angular** spin-8 to spin-4 amplitude ratio. The step from one to the other is
a plausibility argument, and it is not quantified anywhere in this repository.

## What this does not separate, and what would

A spin-8 amplitude comparable to spin-4, and a `C + A4 cos4 + A8 cos8` form that
is simply wrong, give the same arithmetic here. This run cannot tell them apart,
and the frozen design said so before it ran: `not_established_by_this_run`
already listed "the `C + A4 cos4 + A8 cos8` form itself, which two orientations
determine exactly and therefore cannot check."

Two orientations per rung determine `C` and `A4` exactly with nothing left over.
Three determine `C`, `A4` and `A8` together. The design already costed that:
**N=650**, square `|w|² = 650` (`25+5i`, `23+11i`, `19+17i`), rectangular
`|w|² = 325` (`18+i`, `17+6i`, `15+10i`). The note that costed it said this was
"the upgrade to make only if N=290 comes back interesting."

It came back interesting, for a reason nobody anticipated.

## What changed methodologically, which may matter more

The ratio test could not have found this. It formed `A4(4i)/A4(i)` and excluded
the r=2 rung from the decision — and the entire signal is in the r=2 rung. The
contamination that motivated dropping that rung is now the thing under
suspicion, and dropping the rung is precisely what made it unfalsifiable.

Carrying a known systematic as an extra basis direction costs one degree of
freedom. Dropping the entry costs the entry, and with it any chance of noticing
that the systematic is not what you assumed.

## Not established

- which of "large spin-8" and "wrong angular form" is the case;
- that the response is concave at more than about 3σ — one end of the admissible
  covariance range gives 2.62;
- any replacement law. Concave in r is a description of three points, not a
  family;
- anything about `|A8/A4|` itself. This note establishes that the frozen bound
  is unsourced, not that it is false.
