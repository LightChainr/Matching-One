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
| **`α` is removable by a stated principle** | Best case: the `0` at the CM point and the `11/4` rectangular ratio become genuine predictions instead of conditional ones, and the test is much sharper than planned. |

A negative answer here saves a production block, which is the cheaper of the two
outcomes and the more likely one.

## What we already have, and do not need re-derived

- The `-493/192` coefficient and the `A_q/A_ε = (493/96) g₂` Ward ratio, both
  already checked against exact Ward and Hecke fractions in this project.
- The CM zero and the `11/4` rectangular ratio, verified independently by a
  90-decimal direct `E₄` `q`-series.
- The general theory of `g₂`, `E₂`, `E₄`, quasi-modularity and the `E₂` anomaly.
- The dilation cocycle and the `log L = ½ log N` bookkeeping shown above.

## Do not spend output on

- Reviewing modular forms, or listing properties of `g₂`.
- Re-deriving the boxed relations, which we hold as exact.
- Which lattice geometry realizes which `τ` — we handle that.
- Statistical design of the score.

## Provenance of the framing above

`notes/q4-jordan-log-slope-shape.md`, including its own "Boundary" section, which is
the source of the statement that `A~(τ)` is not fixed. The verified rational
constants are frozen in `predictions/q4_jordan_log_slope_shape_20260829.json`.
