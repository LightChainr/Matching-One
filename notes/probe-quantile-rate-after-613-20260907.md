# Probe MID — After `Q_N(u) → p_c`: which rates are theorems from already-named inputs?

**Date:** 2026-09-07
**Assign:** Agent, proof-skeleton probe. No production, no N=725 rescoring, no optimizer, no arXiv search.
**Base:** `claude/matching-one-workspace-pwr5pv` @ `8b5f9d1a` (frontier).

**Tag legend:** **T** = theorem from named inputs; **C** = conditional on a named open hypothesis; **F** = would be a fit (not performed); **X** = not well-posed.

**Working suspicion (proved / falsified / replaced below):**

> No polynomial rate for `Q_N(u) − p_c` follows from H1–H3 plus subcritical exponential decay. The empirical `L^{-4}` root shift is a product of two independent open inputs (F1, F2), not a third mechanism. #582's 9-vector has no right to a single exponent from location alone. An Aff(1)-invariant shape `Z_N` is invisible to Theorem L.

**Verdict:** the suspicion is **confirmed** at the level of tags: R1 is a T (no power), R5/R6 are T, and every positive rate claim is C on F1/F2 (or blocked). No step required a source outside the closed list.

---

## R0 — Quote sheet (no paraphrasing)

All quotes are copied from `notes/p613-quantile-convergence-20260907.md` (PR #614), which is the primary. Section/pointer labels are the note's own.

**Theorem L (verbatim statement, abridged to the operative hypotheses):**

> **Theorem L (repository lemma).** Let `T_N = R^2 / P_N Z^2` be honest periodic square-cell tori with `N` vertices and shortest nonzero Euclidean period `ell_N = min{|v| : v in P_N Z^2, v != 0}`, with
>
> (H0) `ell_N > sqrt(2)` and `ell_N / log N -> infinity`.
>
> ... Assume
>
> (H1) `r_G + r_hat = 2` configuration-wise (digital Alexander duality)
>
> (H2) for each of the two transitive finite-range graphs G (NN) and Ghat (NN+NNN), and every p below its own site critical point, `P_p(0 <-> distance n)` decays exponentially in n
>
> (H3) `p_c^site(G) + p_c^site(Ghat) = 1`
>
> ... for every fixed `u in (0,1)`, `Q_N(u) := F_N^{-1}(u) -> p_c^site(G)`, uniformly on every compact subinterval of `(0,1)`.

**H1 (verbatim):** "`r_G + r_hat = 2` is the repository's own digital Alexander duality (`notes/digital-alexander-duality-proof.md`), proved there for honest periodic square-cell tori. It is **assumed** here, not sourced from the literature."

**H2 (verbatim quotes, with the note's pointers):**

> **Theorem 1.1, item 3** (verbatim): "If `(J_{x,y})_{x,y∈V}` is finite-range, then for any `β < β_c`, there exists `c = c(β) > 0` such that `P_β[0 <-> Λ_n^c] <= e^{-cn}` for all `n >= 0`."
>
> **§1.2, paragraph "Site percolation"** (verbatim): "As in [AB87], the proof may be adapted to site percolation on transitive graphs."

Source line the note gives: H. Duminil-Copin & V. Tassion, *A new proof of the sharpness of the phase transition for Bernoulli percolation and the Ising model*, Commun. Math. Phys. **343** (2016) 725–745; arXiv:1502.03050; DOI 10.1007/s00220-015-2480-z. `[AB87]` = Aizenman–Barsky, Commun. Math. Phys. **108** (1987) 489–526; classical companion Menshikov, Dokl. Akad. Nauk SSSR **288** (1986) 1308–1311.

**H3 (verbatim, with the note's warning):**

> The relation (1.1) of G. R. Grimmett & Z. Li, *Percolation critical probabilities of matching lattice-pairs*, Random Structures & Algorithms **65** (2024) 832–…; DOI 10.1002/rsa.21226; arXiv:2205.02734, reads:
>
> (1.1) `p_c^site(G) + p_c^site(G*) = 1`
>
> (1.2) `p_c^site(G*) <= p_c^site(G)` (trivial)
>
> (1.3) `p_u^site(G) + p_c^site(G*) = 1`
>
> "So (1.1) is **not** a theorem of that paper; the theorem is (1.3), proved in the companion Grimmett & Li, *Hyperbolic site percolation*, arXiv:2203.00981 ... The chain we actually use is (1.3) + `p_c^site(G) = p_u^site(G)` for amenable `G` [Burton–Keane (1989), Comm. Math. Phys. 121, 501–505] => (1.1)."

**van den Berg counterexample (verbatim, with the note's warning):**

> "**⚠ Warning that must stay in the record.** van den Berg's abstract states that Sykes and Essam 'suggested that the above relation holds for all mosaics … **we have constructed a counterexample**', and that the derivation holds only for a restricted class. So (H3) is a genuine hypothesis with **known counterexamples in greater generality**; it is satisfied here specifically because `Z^2` is amenable ... Do not carry (H3) to an arbitrary planar matching pair." (J. van den Berg, J. Math. Phys. **22**(1) (1981) 152–157, DOI 10.1063/1.524747.)

**DKS parent (verbatim):**

> "Paul Duncan, Matthew Kahle & Benjamin Schweinhart, *Homological percolation on a torus: plaquettes and permutohedra*, arXiv:2011.11903, Ann. Inst. H. Poincaré Probab. Statist. (2025) already prove that giant cycles ... have a sharp threshold ... In `d = 2`, `i = 1`, 1-dimensional plaquette percolation **is** bond percolation on the square torus ... Their site-percolation model is the **triangular** lattice, not the square lattice, so the square-**site** case ... is not literally in DKS."

**Union bound + radius (verbatim):**

> "`R_N = floor( ell_N / 4 )` ... `P_p(r_G > 0) <= N * P_p^{Z^2}( 0 <-> Λ_{R_N}^c ) <= N * A(p) * exp(-c(p) R_N) <= N * A(p) * e^{c(p)} * exp(-c(p) ell_N / 4)` ... `N exp(-c'(p) ell_N) = exp( log N - c'(p) ell_N ) -> 0 <=> ell_N / log N -> infinity` ... **Note the rate `c'(p)` depends on `p` and degenerates as `p ↑ p_c`, so the conclusion is per-fixed-`p`, never uniform up to `p_c`.**"

**No-rate / no-value non-claims (verbatim):**

> "* No **rate**. At `p = p_c` the argument says nothing; `F_N(p_c)` may converge to any value in `[0,1]`, and `c'(p) -> 0` as `p ↑ p_c`." and "* No **value**. Theorem L locates the limit at `p_c^site(Z^2)`. It does not compute it ... and says nothing about correction-to-scaling shape."

**Quote-sheet defects found:** none. Every quote claimed in the #613 note is present in the note's own text with a theorem/equation pointer. (No `UNVERIFIED_IN_NOTE` marks needed.)

---

## R1 — No polynomial rate from the union bound (main theorem-shaped negative)

**Theorem R1 (tag T).** Under H1–H3 and DCT as used in #613, for each fixed `p < p_c`,

```
P_p(r_G > 0) <= N · A(p) e^{c(p)} exp(−c(p) ℓ_N / 4) → 0
```

whenever `ℓ_N / log N → ∞`, with `c(p) > 0` allowed to tend to `0` as `p ↑ p_c`. Consequently this bound does **not** imply `|Q_N(u) − p_c| ≤ C_u N^{−ω}` for any `ω > 0`, any `u ∈ (0,1)`, or any `C_u` independent of `N`.

**Proof.** The union bound is copied verbatim from R0; its validity is the note's §1.4 and is not re-derived. The claim to prove is only the negative: that the bound, with `c(p)` an arbitrary positive function vanishing at `p_c`, gives no uniform power.

Fix `ω > 0` and a prospective window `p_N := p_c − N^{−ω}` (any fixed multiplicative constant is irrelevant; take it to be 1). To conclude `Q_N(u) ≥ p_N` for all `u ≥ u_0` one would need `F_N(p_N) → 0`; to conclude `|Q_N(u) − p_c| ≤ N^{−ω}` on compact `u` one needs `F_N(p_N) → 0` and `1 − F_N(p_c + N^{−ω}) → 0` uniformly. It suffices to kill the first.

From R0,

```
P_{p_N}(r_G > 0) ≤ exp( log N + log A(p_N) + c(p_N) − c(p_N) ℓ_N / 4 ).
```

Since `r_G ∈ {0,1,2}`, `F_N(p_N) = E_{p_N}[r_G]/2 ≤ (3/2) P_{p_N}(r_G > 0)`, so `F_N(p_N) → 0` would follow if `c(p_N) ℓ_N / log N → ∞`.

Now `H2` (DCT Thm 1.1(3)) supplies, for each fixed `p < p_c`, a *number* `c(p) > 0`; it supplies **no lower bound** on how `c(p)` behaves as `p ↑ p_c`. The hypothesis is therefore consistent with, e.g.,

```
c(p) = (p_c − p)^{2/ω}.
```

For this admissible `c`, at `p_N = p_c − N^{−ω}` one has `c(p_N) = N^{−2}`, and with the primitive-Gaussian `ℓ_N = √N` (the note's own example, `ℓ_N = √N`),

```
c(p_N) ℓ_N / log N = N^{−2} · √N / log N = N^{−3/2} / log N → 0,
```

so the bound gives `P_{p_N}(r_G > 0) ≤ exp( log N · (1 + o(1)) ) → ∞`, which is vacuous. No conclusion about `F_N(p_N)` follows.

This is the receding-window failure made explicit: the `N` at which the bound becomes small recedes faster than any polynomial, because the only control on `c(p)` is `c(p) → 0`, and that vanishing may be slower than any `(p_c − p)^α`. Since the calculation holds for **every** `ω > 0` with a different admissible `c(p)` each time, no single `ω` is forced.

Equivalently: Theorem L's argument is a family of *pointwise* (fixed-`p`) exponential bounds; pointwise exponential decay of connection probabilities with a rate that is arbitrary in `p` is compatible with any polynomial (or slower) window. A uniform power `|Q_N(u) − p_c| ≤ C_u N^{−ω}` is a statement about a `u`-uniform, `N`-uniform window, and it is exactly what the union bound does not provide. ∎

**Tag:** T (uses only H1–H3 + DCT + the note's `ℓ_N = √N` example). No `ν`, no RSW. **This is the main negative product.**

---

## R2 — What Theorem L *does* give (sharp, T-tagged)

Each item ≤ 10 lines, no extras.

1. **Interior quantile location.** (T) `Q_N(u) → p_c` uniformly on `[ε, 1−ε]` for each `ε > 0`. This is Theorem L verbatim.

2. **Mixed-sign bounded weights inherit location.** (T) Let `w_k` bounded with `Σ w_k = 1`. Then `Σ_k w_k Q_N^{(k)}(u) → p_c`. Three-line proof from the note's "combine quantile functions, not CDFs": quantile functions add pointwise and are all squeezed to `p_c`; a bounded signed combination of `u↦Q_N^{(k)}(u)`, each converging uniformly to the *same* constant `p_c`, converges to `p_c`. (Boundedness keeps the linear combination from amplifying the finite error.)

3. **Unbounded weights fail.** (T) If `sup_N max_k |w_k| = ∞`, the combination can amplify the `O(1)` finite-size error without bound; location is not controlled. This is the note's own amplification argument, stated fully: `Σ w_k (Q_N^{(k)} − p_c)` has no bound if the coefficients are unbounded.

4. **`u_N → 0` or `1` not claimed.** (T, as a scope statement) For `u_N ↓ 0` the proof needs a lower bound on `u_N` against the super-exponential tail `F_N(p_c − ε) ≤ N C e^{−c' ℓ_N}`; not available. `Q_N(u_N)` is only constrained to the correct side of `p_c`.

5. **Thin tori fail.** (T, as a scope statement) `ℓ_N = O(log N)` ⟹ the bound `N e^{−c'ℓ_N}` does not tend to zero; the model is quasi-one-dimensional and has no phase transition at `p_c(Z^2)` at all.

6. **`ℓ_N ≤ √2` not honest.** (T, as a scope statement) For `ℓ_N ≤ √2` the quotient is not honest: a unit cell can have two corners identified and an NN+NNN edge closes into a loop, so (H1) is not stated.

7. **No statement about `F_N(p_c)`, no rate, no operator.** (T, as a scope statement) Copied from R0's two non-claims. `F_N(p_c)` may converge to anything in `[0,1]`; `c'(p) → 0`; no correction shape.

Any corollary not on this list needs a proof from the closed list or is blocked.

---

## R3 — `p_L^H` *is* `Q_L(1/2)`, and `M(1/2)` is the matching-odd content at the self-dual parameter

**Identity (tag T).** Under the declared contract `F_N(p) = [1 + M_N(p)]/2` with `M_N` continuous, strictly increasing, `M_N(0) = −1`, `M_N(1) = +1` (census), the unique root `p_L^H` (i.e. `M_N(p_L^H) = 0`) satisfies `F_N(p_L^H) = 1/2`, hence `Q_L(1/2) = p_L^H`.

**Proof.** `F_N(p_L^H) = (1 + 0)/2 = 1/2`, and `F_N` strictly increasing ⟹ `F_N^{-1}(1/2) = p_L^H`, which is `Q_L(1/2)`. ∎

**Consequence.** The slogan "`p_L^H` is not the median" is **false** under the contract, if "median" means `Q(1/2)`. What *is* true, and what F1/F2 use:

- `p_L^H ≠ 1/2` ⟺ `M_N(1/2) ≠ 0`.
- `Q_L(1/2) − 1/2` is **not** a self-dual vanishing.
- Self-duality / self-matching would force `M_N(1/2) = 0` and `p_L^H = 1/2` (ledger's unconditional self-matching theorem).
- `M_N(1/2)` is the **matching-odd content at the self-dual parameter** — a coordinate different from `Q_L(u) − p_c`.

**First-order expansion (tag T, the ledger's identity).** Taylor at `p_c`:

```
M_N(p_c) + M_N'(p_c)(p_L^H − p_c) + R = 0
⇒  p_L^H − p_c = − M_N(p_c)/M_N'(p_c)  +  higher order.
```

**Conditional consequence (tag C), copied from the ledger with its remainder hypothesis exactly as stated:**

> If `M_L'(p_c) ≍ L^{3/4}` (F1) and `M_L(p_c) ≍ L^{-13/4}` (F2), then `p_L^H − p_c ≍ L^{-4}`, provided the ledger's remainder/error hypotheses hold (numerator `M_L(p_c) = b_L A(τ) + o(b_L)`, denominator `a_L M_L'(p_c + λ a_L) → M'(λ,τ)` uniformly with `M'(0,τ) ≠ 0`, and depth `δ_L/a_L → 0`).

**Tag Outcome J:** this C may **not** be treated as a T. F1 and F2 are not simultaneously available in any current rigorous model (Outcome J).

**Census numbers (cite PR #606's `results/homological-balance-exact-torus/latest.json`; not re-run here):**

| L | `p_L^H` | `M_L(1/2)` |
|---|---|---|
| 3 | 0.586511455113 | `−21/64` |
| 4 | 0.590672112331 | `−13757/32768` |

(Hole: if that JSON is not yet on the frontier, the numbers above are quoted from PR #606's description and the JSON path is left for the exact-controls interface.)

---

## R4 — Conditional statements, locked to named opens

One row per extra hypothesis **already named** in #606/#613. No rows added by searching.

| extra (already named) | implied rate | applies to | tag | still open because |
|---|---|---|---|---|
| F1 only (`M' ≍ L^{3/4}` ⇔ `α_4 ≍ L^{-5/4}`, `ν = 4/3`) | denominator of root shift | `M_L'(p_c)` | C | square-site conformal invariance |
| F2 only (`M(p_c) ≍ L^{-13/4}`, `Q_4 ε`, `x = 21/4`) | numerator of root shift | `M_L(p_c)` | C | operator identification |
| F1 + F2 + ledger remainder | `p_L^H − p_c ~ L^{-4}` | the root, **not** the 9-vector | C | Outcome J |
| RSW / box-crossing as *named* in #606/#613 (triangular / Ising, not square) | whatever those notes already say | **not** square-site `Q_N` | C or X | square site not in the named theorems |
| Russo formula (in-repo, exact) | `M'` is a pivotal mass, no rate | the derivative, not `Q` | T | — |

Row that must **not** be written: "RSW on square site ⇒ quantile window `N^{-3/4}`". That extra is **not named as a square-site theorem** in the closed list.

**BLOCKED_ON_LITERATURE_PROBE:** the implication "already-proved planar RSW ⇒ a quantile-window power, even off square site" — see R10, item 1. Not invented here.

---

## R5 — One exponent for the 9-vector is not implied by location

**Tag T (from R1 + R2 only).**

**Statement.** Theorem L + "no polynomial rate" (R1) ⇒ the 9-vector `(Q_N(0.1), …, Q_N(0.9))` has no theorem-given common `ω`.

**Proof (two lines, the point being the logic).** Uniformity on compact `u` (Theorem L) is a statement about **limits**, not about a common power. R1 already shows that not even a single fixed `u` carries a guaranteed power from the union bound. A fortiori, even if each `u` individually had a rate, the rates need not agree across `u`, and no common `ω` follows. ∎

**Consequence for #582.** #582's "one transferable direction" is therefore **not** a corollary of location. (It may still be fitted as F — that is the HIGH probe's stop rule, not this probe's.)

---

## R6 — `Z_N` is invisible to Theorem L

Let `Z_N(u) = (Q_N(u) − Q_N(0.5)) / (Q_N(0.8) − Q_N(0.2))` (or any Aff(1)-invariant of `Q_N`).

**Lemma (tag T, elementary analysis + Theorem L).** Theorem L constrains `Q_N(u) → p_c` for each fixed `u`, hence **both** numerator and denominator of `Z_N` tend to `0`. The limit of `Z_N` is a `0/0` form and is **not determined** by Theorem L.

**Proof.** Theorem L gives `Q_N(u) → p_c` uniformly on compact `u`; in particular `Q_N(0.5), Q_N(0.2), Q_N(0.8) → p_c`, so numerator `→ 0` and denominator `→ 0`. An Aff(1)-invariant ratio of two quantities that both vanish carries no limit information from their shared limit alone; different subleading behaviors give different `Z_N` limits, and Theorem L is silent on subleading behavior (R1). ∎

**Interface note:** this is the lemma the HIGH probe cites when refusing to interpret `Z_N` shape as a location statement. No `Z_N` is computed on production data here.

---

## R7 — Scaling function as a *definition*, if F1 were true

**Tag C.** Grant F1 hypothetically. Define `τ = L^{1/ν}(p − p_c) = L^{3/4}(p − p_c)` and

```
Φ_L(u, τ) := F_L( p_c + L^{-3/4} τ )   (or the inverse for Q).
```

"`Φ_L → Φ` for a single scaling function `Φ(u, τ)`" would mean:

1. **Adjacent-size first differences at fixed `u`** are then given by `Φ(u, τ)`-differences: `Q_L(u) − Q_{L'}(u)` is determined by `Φ` on the two `τ` values, so first differences are `Φ`-geometry, not free parameters.
2. **Three-size second differences** are the discrete second derivative of `Φ` in `τ`; a genuine `Φ` predicts them, and a one-parameter `Φ` gives a *relation* among them (this is where #609's 1.55 ratio would have to live).
3. **Why a chart (HIGH probe) can still fake a curvature even if `Φ` exists:** a nonlinear chart on `Q` (quantile inversion, CDF↔quantile, affine renormalization of the histogram) contributes its own second-order term `D^2O(θ',θ')` (observer curvature, R8-adjacent), so a nonzero second divided difference is **not** evidence against `Φ` unless the chart is fixed first.
4. **Why F2 (matching-odd at `p_c`) is a *different* coordinate:** F2 is the `τ = 0` value of `M`, i.e. `M_L(p_c)`, a **value** at criticality, not a `τ`-derivative of `F`. It is the numerator of the root shift, decoupled from the thermal denominator F1 that `Φ` encodes.

**Discipline:** no `ν` estimate, no N=725, and #582's amplitudes enter only as "this is the object that would have to be a vector field on `Φ`" — nothing is fitted.

---

## R8 — Russo / pivotal: what it cannot do

From the exact in-repo formula `M'(p) = pivotal_primal(p) + pivotal_matching(1−p)`:

- **(T)** `M'` is a sum of two positive masses; its sign is `+`; this is the strict increase used in R3's identity.
- **(T)** This formula does **not** give `M'(p_c) ≍ L^{3/4}` — that is F1, open.
- **(T)** Integrating `M'` from `1/2` to `p_c` recovers `M(p_c) − M(1/2)`, which is **not** F2 (F2 is the value `M(p_c)`, not the integral against the derivative; the integral identity is exact but circular for F2).
- **(X)** Using this formula to "identify Q4" is a category error: `M'` is a finite pivotal-mass sum, while Q4 is a continuum operator identification; the former has no operator content.

No Monte Carlo of pivotals.

---

## R9 — #615 escape D is this route

**Tag T (scope statement).**

> The homological-balance root is #615's **escape D** (signed observable `X = r−1` + a sign-change theorem: `M_L(p) = P_2 − P_0` crosses from −1 to +1), **not** escape A (critical-sector completeness) and **not** a bounded-task-rank statement. Theorem L is the **location half** of D — quantiles squeezed to `p_c`. F1/F2 are the **rate half** of D and are open. P398 balanced order remains irrelevant to `p_c`.

Concretely, #615's dichotomy: bounded task rank has no threshold content without CSC (its Theorem 1); the escape that actually locates `p_c` is a signed observable whose zero has an independent monotonicity/duality argument. Route D is exactly that. The no-go of #615 is **not** reopened.

---

## R10 — Punch list for the literature probe

Collected `BLOCKED_ON_LITERATURE_PROBE` questions, phrased as questions (not paper requests):

1. **Does any already-proved planar RSW imply a quantile-window power, even off square site?** Specifically: given box-crossing lower bounds (triangular site, or FK-Ising) at criticality, is there a *proved* consequence of the form `|Q_N(u) − p_c| ≤ C_u N^{−ω}` for some `ω > 0`, or only a characteristic-length `L(p)` with `L(p_c + N^{-ω}) ≫ N`? This is the missing premise for the (blocked) "square-site quantile window `N^{-3/4}`" row of R4.

2. **Is there a theorem that inverse-CDFs inherit sharp-threshold windows?** If a CDF `F_N(p)` has a threshold window `|p − p_c| ≤ a_N` (i.e. `F_N(p_c − K a_N) → 0`, `F_N(p_c + K a_N) → 1`), does `Q_N(u) − p_c = O(a_N)` follow uniformly on compact `u`, or does inverse-image distortion add a factor? (The note's quantile passage uses only `F_N(p_c ± ε) → 0/1` for **fixed** `ε`, not a window.)

3. **What finite-size *correction* (not location) is proved on triangular site / FK-Ising for crossing probabilities?** I.e. is there a proved leading correction `N^{-Δ}` to `P_cross` on a model with proved conformal invariance, which could serve as a positive control for the (unproved, square-site) F1/F2 product? This matters because R7 needs a `Φ` where both the denominator (F1) *and* a numerator correction are rigorously available in the same model, and the ledger's Outcome J says no such model is currently named.

(These are handover questions only; none are answered here.)

---

## Stop-rule and claim-boundary record

- R1 did **not** fail: a polynomial rate does **not** follow, so the probe proceeds (no paper-killing required).
- No source outside the closed list was used; the single blocked implication (R4/R10 item 1) is recorded as `BLOCKED_ON_LITERATURE_PROBE`, not answered.
- No `docs/STATUS.md` edit, no issue close, no Monte Carlo, no `ν` fit.
- Tags: R1 T, R2 T (×7), R3 T (+ one C), R4 table (T/C/X + one blocked), R5 T, R6 T, R7 C, R8 T (×3) + X, R9 T, R10 handover.

---

*Sources (closed list, already named): #613 note `notes/p613-quantile-convergence-20260907.md` (PR #614), #606 ledger `notes/homological-balance-root-ledger-20260906.md`, #615 note `notes/bounded-task-rank-threshold-no-go-20260907.md`, in-repo Russo formula `M'(p) = pivotal_primal(p) + pivotal_matching(1−p)`.*
