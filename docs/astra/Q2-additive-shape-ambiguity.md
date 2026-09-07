# Q2 — What is the largest amplitude-free modular observable, given that the additive shape is not fixed?

**Blocks:** `docs/ROADMAP.md` item 2, the modulus/shape fingerprint intended to
separate a Jordan module from ordinary logarithmic mixing.

---

## Setup

A rank‑2 pair `(q, q~)` with `D q = x q`, `D q~ = x q~ + q` transforms under a finite
dilation as `q~ ↦ s^{-x}(q~ - log(s) q)`. A top-field response on a torus of modulus
`τ` and linear size `L`, with `N` sites and `log L = ½ log N`, therefore takes the
form

```text
L^{x-2} δO_top(L, τ) = A~(τ) - λ_top A_q(τ) log L,
```

giving the logarithmic coefficient

```text
B_logN(τ) = -(λ_top / 2) A_q(τ).
```

The ordinary torus Ward identity fixes the ratio of the two amplitude functions,

```text
A_q(τ) / A_ε(τ) = (493/96) g₂(τ),
```

so before an unknown non-universal overlap constant the log coefficient carries the
rational factor `-493/192`. If the residual-to-slope/root normalization removes the
same thermal-primary block as elsewhere in our construction, then

```text
B_root(τ) = C_J · Re g₂(τ)      (or C_J · g₂ for a chiral projection),
```

which freezes, in the area-normalized convention,

```text
B_root(e^{iπ/3}) = 0,          B_root(2i) / B_root(i) = 11/4.
```

Those two are the fingerprint we intended to measure.

## The obstruction

The **additive shape `A~(τ)` is not fixed by the Jordan relation.** It shifts under
the redefinition `q~ → q~ + α q` — which is the residual gauge freedom of any rank‑2
pair — and it also depends on generic‑`c` derivatives of the colliding energy and
hull blocks. Nothing in the Jordan structure alone selects an `E₂`, a `log η`, a
weight derivative or a modular derivative.

So the predicted quantity we can actually observe is contaminated by a function we
cannot compute, and the clean-looking `0` and `11/4` above are only available if the
normalization really does remove the same block — which is an assumption, not a
theorem.

## We have now measured it, and it is not 11/4

**Added 2026-09-05, after this question was written.** The two-modulus block described
below was frozen and run. The measurement, on a square torus (`τ = i`) and a 2:1
rectangular torus (`τ = 2i`) with the same 290 sites, of the spin-4 amplitude of the
matching-odd slope, is

```text
ratio  =  A(2i) / A(i)  =  1.880 +/- 0.177          (prospective, frozen before the run)
```

against `Ê4(2i)/Ê4(i) = 11/4 = 2.75`. That is **4.9 sigma low**. Plain area scaling
(4) is excluded at 12 sigma and no modulus dependence (1) at 5.0 sigma, so the number
is not consistent with the obvious alternatives either. A known spin-8 leakage in the
orientation projector (`±1148/21025`, exactly equal and opposite between the two
families) would need `A₈/A₄ = 3.44` to account for the gap, which this project's own
H4-beats-H8 results rule out.

So the conjunction *{weight-4 shape} ∧ {the normalization removes the same block}* is
false. Since the second conjunct is exactly the assumption this question is about,
the question is now sharper and more valuable, not less:

> **What additive-shape contribution `A~(τ)`, admissible under the `q̃ → q̃ + αq`
> gauge and the generic-`c` derivatives of the colliding blocks, would carry the
> ratio from `11/4` to `1.88 ± 0.18` at `τ = 2i` relative to `τ = i`?**
>
> Equivalently: is `1.88` reachable inside the rank-2 module at all, or does it
> falsify the module rather than the normalization?

That is a number to explain, not a design to approve, and it is worth more than the
original framing. Answer it together with the question as originally posed below.

One further datum, for whatever it is worth: `1.880 ± 0.177` is compatible with the
bare aspect ratio `r = 2` (`z = -0.68`). We record that as a post-hoc observation with
no prospective standing — it is a reading of a single point — and we are **not**
asking to have it rationalized. If a linear-in-`r` law falls out of a correct
treatment, that is informative; if it does not, please say so plainly rather than
constructing one.

We have since frozen the design that gives that reading its one chance to lose:
an aspect ladder `r = 1, 2, 4` at a single site count, with the bare aspect ratio
named as a competitor before any block runs
(`predictions/aspect_ladder_n1300_20260905.yaml`). At `r = 4` the weight-4 shape
predicts `10.99` and the bare ratio `4.00`. So the measurement is being taken
either way, and an answer here is worth more before it than after.

## The question

> Fix the two-parameter ambiguity `(α, C)`: `α` the Jordan gauge `q~ → q~ + α q`, and
> `C` an unknown non-universal overlap constant multiplying the whole response.
> Over a set of moduli `{τ₁, …, τ_k}`, **what is the largest set of functionally
> independent combinations of the observable data `{A~(τ_i), B_logN(τ_i), B_root(τ_i)}`
> that is invariant under both?**
>
> And: is there any additional constraint — modular covariance of the full torus
> one-point function, the `c → 0` limit procedure itself, or a canonical choice
> forced by the collision — that removes the `α` freedom, rather than merely
> quotienting by it?

The first half is an algebra question with a definite answer and tells us what to
measure. The second half is the physics question and would give a sharper test.

## Why the answer decides something

We are choosing between spending a production block on a two-modulus shape score and
not spending it. The score is only worth running if some invariant combination
actually separates a Jordan module from ordinary logarithmic mixing.

| Answer | What we do |
|---|---|
| **Invariants exist and separate** | Freeze the smallest such combination as the score and run the two-modulus block. This is `ROADMAP` item 2, and it becomes executable. |
| **Invariants exist but do not separate** | The shape route cannot identify the module. Drop item 2, and say so in the status ledger rather than leaving it "ready". |
| **`α` is removable by a stated principle** | The `0` at the CM point and the `11/4` rectangular ratio become genuine predictions instead of conditional ones — and since `11/4` is now excluded at 4.9 sigma, that outcome would falsify the module in this channel rather than confirm it. Either way it is decisive. |

The block has since been spent: it cost about 70 minutes and returned the number
above, so the table's first two rows are now about what to do with a measurement
rather than whether to take one.

## What we already have, and do not need re-derived

- The `-493/192` coefficient and the `A_q/A_ε = (493/96) g₂` Ward ratio, both
  already checked against exact Ward and Hecke fractions in this project.
- The CM zero and the `11/4` rectangular ratio, verified independently by a
  90-decimal direct `E₄` `q`-series. The competitor shapes and what each predicts for
  the same ratio are in `results/modulus-shape-discrimination/latest.json`; the
  measurement is in `results/modulus-fingerprint-n290/latest.json`.
- The general theory of `g₂`, `E₂`, `E₄`, quasi-modularity and the `E₂` anomaly.
- The dilation cocycle and the `log L = ½ log N` bookkeeping shown above.

## Do not spend output on

- Reviewing modular forms, or listing properties of `g₂`.
- Re-deriving the boxed relations, which we hold as exact.
- Which lattice geometry realizes which `τ` — we handle that.
- Statistical design of the score, or how the measurement was taken. The number
  above is frozen, prospective and ours to defend.

## Provenance of the framing above

`notes/q4-jordan-log-slope-shape.md`, including its own "Boundary" section, which is
the source of the statement that `A~(τ)` is not fixed. The verified rational
constants are frozen in `predictions/q4_jordan_log_slope_shape_20260829.json`.
