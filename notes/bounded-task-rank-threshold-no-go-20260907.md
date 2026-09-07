# Bounded finite-horizon task order carries no threshold content without completeness

**A no-go theorem, its strengthening, and its escape conditions.**
Probe analysis of `LightChainr/Matching-One`, issue #594 Q3 (and #598/#599/#607 context).

---

## 0. Setup and definitions

**System.** For each finite size `L` we have a state space `X_L` (`|X_L| = n_L`),
a Markov generator `G_L(p)` (row-stochastic: nonnegative off-diagonal, zero row
sum), a source matrix `B_L ∈ R^{n_L × m}` and a readout matrix `C_L ∈ R^{q × n_L}`.
The finite-horizon response is

```
R_L(t,p) = C_L exp(t G_L(p)) B_L,    0 ≤ t ≤ T.
```

We say the family is **analytic in `p`** if every entry of `G_L(p)` is an
analytic function of `p` on `[0,1]` (finite matrices of analytic entries are
analytic; spectral gaps are then analytic away from eigenvalue crossings, and
here we will work with a family whose gap has an explicit closed form).

**Task order.** Two notions, which we separate explicitly.

* **Exact task order** `r_ex(G,B,C)` = dimension of the controllable-and-observable
  (Kalman) subspace, i.e. the minimal realization order of the map
  `u ↦ {C exp(tG) B u : t ∈ [0,T]}`.
* **Finite-horizon balanced/effective order** `r_bal(G,B,C,T; ε)` = the number of
  Hankel singular values of the finite-horizon Hankel operator above tolerance `ε`.

**Lemma 0 (analyticity collapses horizon).** *If `G` is a finite matrix then for
any `T > 0` the restriction of the response to `[0,T]` determines the full
response, because `t ↦ C exp(tG) B` is real-analytic; equality on `[0,T]` implies
equality of all Taylor (Markov) coefficients `C G^k B`. Hence*
`r_ex(G,B,C)` *is independent of `T > 0`, and exact equality of two responses on
`[0,T]` is equivalent to equality of all Markov parameters.*

Consequently the no-go below holds verbatim for the exact task order; for the
balanced order the same statement holds as an inequality `r_bal ≤ dim(X_vis)`
at every tolerance and every horizon.

**Definition (closing-gap threshold).** Let `γ_L(p)` be the spectral gap of
`G_L(p)` (the slowest nonzero relaxation rate). We say `p_c ∈ (0,1)` is a
**closing-gap threshold** for the family if

1. for each finite `L`, `p ↦ γ_L(p)` is analytic on `[0,1]` and `γ_L(p) > 0`;
2. the pointwise limit `γ_∞(p) := lim_{L→∞} γ_L(p)` exists, and
   `γ_∞(p) > 0` for `p < p_c` with `lim_{p ↑ p_c} γ_∞(p) = 0`.

This is exactly the abstract form of a subcritical exponential-decay rate
vanishing at `p_c` (the percolation `ξ → ∞` signature). The singularity lives
only in the `L → ∞` limit; every finite `L` is smooth.

---

## 1. Main verdict

**BOUNDED_TASK_RANK_HAS_NO_THRESHOLD_CONTENT_WITHOUT_COMPLETENESS.**

More precisely, the three options reduce to the first one, with a sharpening:

> *For a declared task `(G_L, B_L, C_L)`, the statement `sup_L r_task(L) < ∞`
> implies **nothing** about any threshold `p_c` — unless the critical modes are
> forced to be controllable-and-observable by an extra assumption (critical-sector
> completeness, §4). Without that assumption one can attach, to any bounded-rank
> task, a hidden sector carrying a closing-gap threshold at an **arbitrary**
> prescribed location `p*`, with the response unchanged to machine precision.*

The question is therefore **not malformed**; it is one-sided: bounded task rank
is a statement about *compressibility of a declared I/O map*, while `p_c` is a
statement about *the full generator's spectral edge*. These are orthogonal
coordinates unless a completeness condition couples them.

---

## 2. Strongest theorem / counterexample

### 2.1 The Kalman fact that drives everything

**Proposition 1 (hidden sector ⟺ off-CO).** *Let `W ⊆ X_L` be a `G_L`-invariant
subspace. The response `R_L(t,p)` is independent of the restriction of `G_L` to
`W` (i.e. of `G_L|_W`) if and only if*
`range B_L ∩ W^⊥-orthogonal` *and* `W ⊆ ker C_L` *— concretely, `W` lies in the
uncontrollable-unobservable (Kalman) subspace of `(G_L, B_L, C_L)`. Under this
condition every mode in `W` contributes a partial-fraction coefficient zero to
`R_L`, for all `p, t`.*

This is the exact content of "hidden sector": a hidden mode is one outside the
controllable-and-observable subspace.

### 2.2 Theorem 1 (no-go, direct sum)

**Theorem 1.** *There exist two analytic Markov families `{G_L^{(1)}(p)}` and
`{G_L^{(2)}(p)}` with common source/readout `(B_L, C_L)`, such that*

```
R_L^{(1)}(t,p) = R_L^{(2)}(t,p)   for all L, p, t ∈ [0,T],
```

*yet the families have different closing-gap thresholds*
`p_c^{(1)} = 1/2 ≠ 1/3 = p_c^{(2)}`. *Consequently*
`sup_L r_ex(L) = 2 < ∞` *for both, while `p_c^{(1)} ≠ p_c^{(2)}`.*

**Proof (constructive; verified in `scripts/no_go_theorem.py`).**

*Visible sector.* A fixed 2-state generator
`G_vis = [[-1,1],[1,-1]]` with `B_vis = (1,0)^⊤`, `C_vis = (1,0)`.
Its response `(1+e^{-2t})/2` has task order 2 and is `(L,p)`-independent.

*Hidden sector.* On `{1,…,L}` take the biased nearest-neighbour walk with
reflecting boundaries: right-jump rate `a`, left-jump rate `b`. Its generator is
tridiagonal and its spectral gap has the exact closed form

```
γ_L(a,b) = a + b − 2√(ab)·cos(π/L).
```

This is a connected, positive, nearest-neighbour (hence *local*) Markov chain,
irreducible for every finite `L`, and entrywise analytic in `(a,b)`.

Parameterise `a = e^{φ}`, `b = e^{−φ}` with `φ = p − p_c`. Then

```
γ_L(p) = 2 cosh(p − p_c) − 2 cos(π/L)  →  γ_∞(p) = 4 sinh²((p−p_c)/2),
```

which is `> 0` for `p ≠ p_c` and vanishes exactly at `p = p_c`. Each finite-`L`
gap `γ_L(p) > 0` and is analytic in `p`; the closing occurs only in `L → ∞`.

*Assembly.* Take
`G_L^{(j)} = G_vis ⊕ H_L(p − p_c^{(j)})` with
`B_L = (B_vis ; 0)`, `C_L = (C_vis, 0)`. The response is exactly the visible one,
independent of `j, L, p`, while the hidden gap closes at `p_c^{(j)}`. Two families
(`p_c = 1/2` and `p_c = 1/3`) then have identical responses and different
thresholds. ∎

**Corollary 1 (the no-go).** *`sup_L r_task(L) < ∞` places **no** constraint on
where a closing-gap threshold sits. Equivalently: task compressibility
(`r_bal ≈ 3–4`) and threshold identifiability (`p_c`) are orthogonal observables
in the absence of a completeness assumption.*

> **Remark (why this is not "obvious").** The content of Theorem 1 is not that a
> sector *can* be hidden — Proposition 1 already says that — but that the hidden
> sector can host a **thermodynamic singularity that exists only in the `L → ∞`
> limit while every finite `L` stays analytic, irreducible and local**, and that
> this can be done at an **arbitrary** `p*`. That is the exact statement needed to
> rule out "bounded rank ⇒ some `p_c` information" without an extra hypothesis.

---

## 3. Is the direct sum merely a trivial counterexample?

Direct sum is the *most economical* witness, but the no-go is not an artifact of
block-diagonality. Three levels, in increasing nontriviality:

**(a) Direct sum — the clean theorem.** This is Theorem 1. It is "trivial" only in
the sense that the hidden sector is an explicit invariant summand. It is not
logically trivial: it is the precise certificate that the required missing
hypothesis is *completeness* (off-CO modes carry zero response information).

**(b) Symmetry-protected dark sector — a non-block-diagonal but reducible
precedent.** If a finite group `K` acts on `X_L` and `G_L` commutes with `K`, then
`G_L` is block-diagonal in the isotypic decomposition of `K`. A source/readout
confined to one isotypic type sees none of the other types. **P398 is exactly
this**: a large reflection-odd sector is *strictly uncontrollable/unobservable*
for reflection-even source/readout (the parity selection rule; verified in the
earlier probe rounds to ~1e-17). In the original join/detach basis `G` is **not**
block-diagonal — the hiding is representation-theoretic, not a coordinate
coincidence. So the no-go survives "no explicit block decomposition" in the
relevant sense.

> **Important structural fact.** An irreducible Markov generator cannot carry a
> *nontrivial* commuting finite group action with more than one isotypic type:
> commuting with a nontrivial `K` forces a nontrivial invariant decomposition,
> hence reducibility. Therefore **symmetry-protected hiding is necessarily a
> reducible mechanism**. "Irreducible" and "nontrivial symmetry" are mutually
> exclusive here.

**(c) Irreducible — exact masking is fine-tuned but real.** On a connected,
positive walk one can still make the *slowest mode* exactly invisible by choosing
the readout `C` orthogonal to the slow right eigenvector `v` (`C v = 0`) or the
source `B` orthogonal to the slow left eigenvector `w` (`w^⊤ B = 0`): the
partial-fraction coefficient of `e^{λ₁ t}` in the response is `(C v)(w^⊤ B)`,
which vanishes exactly. Verified on a 4-state biased walk:

| configuration | `(C v)(w^⊤ B)` |
|---|---|
| generic `C=e₁, B=e₁` | `4.84e-1` |
| `C v = 0` | `4.13e-17` |
| `w^⊤ B = 0` | `−2.31e-17` |
| both | `−1.97e-33` |

This is *fine-tuning* (generic readout/source see the slow mode), not a robust
sector, and — crucially — it masks a single slow mode rather than letting its
threshold be relocated while the response stays constant. So:

> **Precise verdict on "irreducible / local / positive / connected".**
> *Exact* hiding survives only as (i) a symmetry-protected reducible sector
> (P398), or (ii) a fine-tuned single-mode orthogonality (irreducible). What
> survives **generically** in the irreducible/local/positive class is **weak
> coupling**: a hidden sector coupled to the visible sector with strength `ε`
> perturbs the response by `O(ε)` and its balanced order by `O(ε)`, so at any
> finite tolerance the no-go still holds in the *balanced/effective* sense, with
> the threshold again freely placeable.

Hence the no-go is not an artifact of direct sum; it is the exact statement that
**threshold information can only enter a declared task through controllable-and-
observable modes** — and neither symmetry, locality, positivity nor irreducibility
provides such a mode by itself.

---

## 4. Escape conditions — what is actually necessary

The candidates, assessed against "necessary vs sufficient".

**A. Critical-sector completeness (CSC) — necessary.**

> `H_crit(L) ⊆ H_CO(G_L, B_L, C_L)` for each `L`, where `H_crit(L)` is the span of
> the modes whose relaxation rate `→ 0` as `L → ∞`.

*Status:* **necessary.** Proposition 1 is exactly the contrapositive: a critical
mode outside `H_CO` contributes coefficient zero to the response and its threshold
is invisible. This is the *minimal* hypothesis that closes the no-go — without it
Theorem 1 applies verbatim. *Not sufficient*: CSC guarantees the critical modes
are *reachable*, not that a *low-rank* realization carries them; a task can be
CSC-complete yet need high order to resolve the slow edge.

**B. Horizon growing with correlation length (`T_L → ∞`, `T_L ≍ ξ_L`) — sufficient,
not necessary.**

For the **exact** order `T` is irrelevant (Lemma 0), so this cannot rescue the
exact-rank question. For the **balanced/effective** order, growing `T` raises the
Hankel singular values of the slow (but CO-reachable) modes, so a fixed tolerance
eventually *reveals* them — *provided* they are already in `H_CO`. `T_L → ∞`
therefore **amplifies** CSC-complete critical modes but cannot create coverage
that is not there. It is a sufficient amplifier, not a necessary condition.

**C. Uniform spectral approximation — sufficient, strictly stronger than A.**

Requiring `|R_L − R_L^{(r)}| → 0` *uniformly on a window containing the slowest
mode* is (A) + "the low-order model retains the slowest mode's dynamics". It is
sufficient, but it is a *tautological strengthening*: it assumes the slow mode is
both in `H_CO` **and** kept by the approximation. This is not the minimal
assumption.

**D. The observable itself defines the transition — a different logic entirely.**

`M_L(p) = P_2 − P_0` with `M_L(p_L^H) = 0` and an independent
topology/matching argument for `p_L^H → p_c` does **not** go through rank at all.
It works via *monotonicity / subcritical decay / matching complement / digital
Alexander identity / topological phase separation* — i.e. a **zero of a signed
observable with a sign-change argument**, not compressibility. This is the route
that actually yields `p_c`, and it is orthogonal to (A)–(C). The two statements

* "the observable has a threshold theorem", and
* "the observable is easy to realize at low order",

are **logically independent**. One can be true while the other is false.

**E. Experiment-language completeness — the right framing, reducible to B.**

`r_d(k) = min(d+1, k+1)` (from #607) says a *fixed experiment depth* stays low-rank
while a *richer future language* opens hidden complexity. Threshold identifiability
is the same phenomenon at the level of the whole family: the critical information
lives at depth/horizon `≍ ξ_L`, invisible to fixed-depth languages. The formal
version is exactly B with "depth" replaced by "horizon": critical modes become
task-visible only when the declared experiment reaches the critical scale.

> **Bottom line.** The single necessary condition is **A (critical-sector
> completeness)**. B and E are amplifiers of A; C is a stronger sufficient form;
> D is an independent mechanism that sidesteps rank. None of B–E removes the need
> for A.

---

## 5. P398 consequence

**What P398 establishes.**

> `r_balanced ≈ 3–4`, nearly constant across width, is the statement that the
> *declared task* — the specific `(sources, readouts, horizon)` chosen — has low
> effective input–output complexity. It is a property of a **triple `(G,B,C)`**,
> not of `G` alone.

**What P398 does not establish.**

> It does **not** establish that the underlying process has a low-dimensional
> critical state, nor that any thermodynamic threshold is small, large, or even
> present. By Theorem 1, a bounded task order is compatible with an arbitrary
> hidden critical point.

**Is P398 therefore "useless for percolation"?** No — the correct statement is
scope-separated:

| use | does P398 help? |
|---|---|
| threshold *estimation* of `p_c` | **No direct contribution** (Theorem 1) |
| observability / controllability diagnostics | **Yes** — it is a clean, exact *dark-sector* testbed |
| state semantics ("what is a state") | **Yes** — it separates microscopic / quotient / task / balanced notions (the 1430 > 750 > 209 > 32 > 3–4 ladder) |
| experiment design (which readouts open which sector) | **Yes** — the parity selection rule is a worked example of the exact tensor criterion |
| hidden-sector selection rules | **Yes** — the finite-group response theorem is precisely this |

**One-sentence scope statement (ready for a paper):**

> *"A bounded finite-horizon balanced order of a declared task certifies the
> compressibility of that task's input–output map and the existence of an exact or
> symmetry-protected uncontrollable/unobservable sector; it carries no information
> about any thermodynamic threshold unless the critical modes are separately shown
> to be controllable and observable."*

---

## 6. Square-site consequence

An observable/task is entitled to carry threshold information **only if**

1. **critical-sector completeness** holds: the modes whose rate closes at `p_c`
   lie inside `H_CO(G_L, B_L, C_L)` — otherwise they are invisible by Proposition 1;
2. preferably, the observable **itself has a sign-change / monotonicity theorem**
   (the `M_L(p) = P_2 − P_0` route), because that is what converts "reachable" into
   "locates `p_c`" — rank alone never does.

The real threshold route in Matching-One (`X = r−1`, `M_L(p) = E[X] = P_2 − P_0`,
homological-balance root `M_L(p_L^H) = 0`, `p_L^H → p_c`) does **not** draw its
strength from low Hankel rank. It draws it from: monotonicity of `M_L`, subcritical
exponential decay, the matching complement, a digital Alexander identity, and
topological phase separation. Those are *spectral-edge / signed-observable* facts,
not *compressibility* facts.

This is the strict meaning of

```
task compressibility  ≠  threshold identifiability.
```

The first is a statement about `(G,B,C)`; the second is a statement about `spec(G_L)`
and about a signed observable with a sign-change argument.

---

## 7. #609 companion recommendation

#609's single-exponent amplitude law `A(N) ~ N^{−0.970}` under-predicts by ~55%
in second differences with ratios `1.5368, 1.5574` in two lineages. The open
question is whether the residual is (1) a genuine second scaling mode or (2) a
step-size-dependent quantile-reconstruction bias. A changed-readout control
separates them, and the observability framing above says exactly how.

**Principle.** (1) is a property of the *generator's spectrum* (process-level);
(2) is a property of the *readout map* (observer-level). A readout change that
keeps the critical mode's projection on `H_CO` fixed but changes the
reconstruction procedure will move (2) and leave (1) fixed.

**Minimal design.**

* Hold the process and `N` fixed.
* Vary only the **readout coarsening**: e.g. change the quantile grid / bin
  definition / `N`-normalisation used to reconstruct the law — keeping the same
  critical observable but a *different* reconstruction operator `C → C′`.
* Compute the second-difference discrepancy ratio under `C` and `C′`.

**Decision rule.**

* If the ratio `≈ 1.55` is **readout-stable** → genuine second mode (process-level).
* If it **moves** with `C′` → observer-induced curvature (reconstruction bias).

The "changed readout" that is most diagnostic is one that *does not commute* with
the reconstruction's `N`-dependence: keep the critical observable identical, alter
the finite-`N` quantile/binning recipe. This is the minimal experiment, and it
costs no new Monte Carlo — only a re-read of existing data through a second
reconstruction operator.

---

## Appendix. A canonical coverage quantity (optional)

The quantity that directly measures whether the declared experiment "covers" the
critical sector is *not* the balanced order; it is the projection of the critical
modes onto the controllable-and-observable subspace. Canonically:

```
c_crit(L) := ‖ P_CO(L) · Π_crit(L) ‖
```

where `Π_crit(L)` is the projector onto the slow subspace and `P_CO(L)` the
Kalman CO projector. Equivalently, `c_crit` is the largest Hankel singular value
attributable to the slow modes. `c_crit(L) = 0` is exactly the no-go regime;
`c_crit(L) → 1` is exactly CSC-completeness. This is a derived, non-fabricated
indicator: it is the object the no-go theorem says *must* be checked before any
rank statement may be interpreted as threshold information.

---

### Files

```
notes/bounded-task-rank-threshold-no-go-20260907.md   (this document)
scripts/no_go_theorem.py         Theorem 1 + witness (exact, reproducible)
scripts/irreducible_masking.py   §3(c) fine-tuned exact masking
results/no-go/latest.json
results/irreducible-masking/latest.json
```
