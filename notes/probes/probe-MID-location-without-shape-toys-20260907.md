# Probe MID — Location without shape: explicit toy families

**Assign to:** ordinary (adequate) reasoning, ordinary CPU. Analysis + numpy/mpmath. No percolation engine, no arXiv, no N=725, no #612 rescoring.

This is a **construction** probe. The deliverable is families you invent, not a match to a repository JSON.

---

## Standing

- Does **not** enter `docs/STATUS.md`. Does not close #613, #618, #622, #276.
- Does not search arXiv (#620). Does not prove square-site rates (#618). Does not define the percolation `S` (#622); you supply the *analytic* counterweight that #622’s W1 asked for.
- Submit a **new** PR against `claude/matching-one-workspace-pwr5pv`. Comment the URL on the issue. Leave the issue open. Not against #616.

Frontier: `claude/matching-one-workspace-pwr5pv` @ `8b5f9d1a`.
Theorem L as a cited fact: `notes/p613-quantile-convergence-20260907.md` — location, not rate. You do not re-prove it.

---

## North-star

Exhibit explicit sequences of CDFs `F_N` (or quantile functions `Q_N`) on `[0,1]` such that

```text
Q_N(u) → p_*     uniformly on every compact [ε, 1−ε] ⊂ (0,1)
```

for some `p_* ∈ (0,1)`, **while** the Aff(1)-invariant shape

```text
Z_N(u) = (Q_N(u) − Q_N(a)) / (Q_N(b) − Q_N(a))
```

(with declared anchors, default `a=0.2`, `b=0.8`) does one of the following, **by construction**:

| target | what you must exhibit |
|---|---|
| T1 | `Z_N` fails to converge (oscillation, two accumulation points, or `N`-dependent anchors that still oscillate after a fixed-anchor definition) |
| T2 | `Z_N` converges, and the limit `Z_∞` can be prescribed: for any continuous strictly increasing `ζ` with `ζ(a)=0`, `ζ(b)=1`, there is a family with that `Z_∞` |
| T3 | `Z_N` converges, but the limit **depends on the anchors** `(a,b)` in a way that is not a reparametrisation (then “the” shape is not well-posed even for toys) |

If T2 is achieved, W1 of #622 becomes a theorem about location statements in general: **Theorem L cannot pin percolation’s shape**. That is the point of this probe.

If you believe T2 is false for monotone CDFs with a single crossing of every `u∈(ε,1−ε)`, prove that restriction and still do T1.

---

## Constraints on the toys (so they are fair analogues)

Each family must be:

- a CDF in `p` for each `N`: `F_N(0)=0`, `F_N(1)=1`, nondecreasing (cadlag allowed; say so);
- strictly increasing on a neighbourhood of `p_*` for large `N`, so `Q_N(u)` is a singleton for `u∈[ε,1−ε]`;
- **not** a mixture of two distant jumps that “cheats” location by parking mass at `0` and `1` only — at least one of the families must have a window of width `w_N → 0` about `p_*` that carries all interior quantiles (this is what Theorem L actually looks like);
- parametrised in closed form, or by a 20-line sampler of a named distribution. No black-box neural CDF.

Write the two Aff(1) actions and say which one your `Z` quotients:

```text
on p:   p ↦ αp+β
on Q:   Q(u) ↦ α Q(u)+β
```

`Z` as written quotients the **Q-action**. One family should also be inspected under a warp of the `p`-axis, to show the two quotients disagree (a numerical example is enough; the theorem is #622 W5).

---

## Directions (all three targets attempted; T2 is the prize)

### D1 — Logistic / probit window (baseline that *does* pin a shape)

```text
F_N(p) = σ( (p − p_*) / w_N ),    w_N → 0,   σ a fixed logistic or Gaussian CDF.
```

Compute `Z_N`. It will converge to a **fixed** logistic/probit shape independent of `w_N`. This is the family people silently imagine when they say “sharp threshold ⇒ universal shape”. Record it as the **non-example**: location plus a *fixed* window profile pins `Z`. Theorem L does not give you a fixed profile.

### D2 — Prescribed shape (T2)

Replace `σ` by a sequence `σ_N`, or by a mixture of two window profiles whose weights oscillate, or by an `N`-dependent skew (e.g. GEV / skew-logistic with a parameter `γ_N`).

Goal: pick any target `ζ` (give three: symmetric logistic; a strongly skew one; a piecewise-linear “kink”), build `F_N` with window `w_N = N^{-1}` or `N^{-3/4}` (the exponent is free in a toy), prove `Q_N(u)→p_*` uniformly on `[ε,1−ε]`, and prove `Z_N → ζ`.

If a monotone one-parameter exponential family cannot do T2, say so and use a two-parameter window (location + skew), still with width → 0.

### D3 — Oscillation (T1)

Two profiles `σ` and `τ` with different `Z_∞`, switched on even/odd `N`, or a skew parameter `γ_N = sin(log N)`. Prove location still holds. Then `Z_N` has no limit. This kills “sharp ⇒ shape converges”.

### D4 — Anchor pathology (T3)

Using one D2 family, recompute `Z` with anchors `(0.1,0.9)` vs `(0.3,0.7)` vs `(0.2,0.8)`. If the three `Z_∞` differ by more than reparametrisation of `u` (make that precise: they do not lie on one orbit under increasing maps that fix `{a,b}` — or they do). Either “anchors are gauge” or “the invariant is not well-posed”. One page.

### D5 — What extra input would pin `ζ` in percolation

After T1/T2, one page, no literature search. Named extra inputs already in the #613/#606 notes (RSW, four-arm / F1, self-duality). For each: would it pin `σ` (the window *profile*), or only the *width*? Write Y/N/unknown. Unknown is allowed. Do not quote new papers.

---

## Compute

Grid `u = 0.05(0.05)0.95`, `N` in `{2^k}` over at least two decades of `w_N`. JSON of `Q_N`, `Z_N` for each family. Plots optional as committed png under `results/`, not required.

No Monte Carlo of percolation. No L=3 enumerator.

---

## Stop rules

- If T2 succeeds, **stop inventing more families**. Write the theorem-shaped statement and D5. That is a complete probe.
- If T2 fails for all monotone window families you tried, prove a restricted negative (“one-parameter location-scale windows have unique Z_∞ determined by σ”) and still deliver T1.
- Do not fit these toys to #582 amplitudes. Do not mention 1.55 except as “retired, irrelevant here”.

---

## Deliverables

```text
notes/probe-location-without-shape-YYYYMMDD.md
  formulas, proofs of location, T1/T2/T3 verdict, D5 table
scripts/probe/toy_cdf_families.py
results/probe-location-without-shape/latest.json
```

## Interface

#622 uses this as the W1 laboratory. #618 may cite T2 as “location does not imply a common ω for the 9-vector” — only if you actually proved T2. #620 is not involved.
