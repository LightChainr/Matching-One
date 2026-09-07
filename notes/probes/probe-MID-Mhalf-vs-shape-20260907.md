# Probe MID — Does `M(1/2)` pin the shape on exact tori?

**Assign to:** ordinary (adequate) reasoning, ordinary CPU. Python + the in-repo L=3,4 enumerator. Bond L=3 (`2^{18}`) is in budget. Bond L=4 is not. No arXiv, no production, no N=725.

This is a **kill/promote** probe for #622’s W4: the odd moment `M(1/2)` either governs `Z_L` or it does not. The census numbers `M(1/2) = -21/64` etc. are **inputs**, not deliverables. Do not re-verify PR #606 as the job.

---

## Standing

- Does **not** enter `docs/STATUS.md`. Does not close #606, #619, #622, #608.
- Import `scripts/homological_balance/exact_torus_enum.py` (PR #606). Do not duplicate it. If dest paths are wrong, fix `parents[2]` only as needed to run.
- If #619 has already produced `M` polynomials and `Q(u)` tables, **import them** and do not recompute; start at D2.
- Submit a **new** PR against `claude/matching-one-workspace-pwr5pv`. Comment the URL on the issue. Leave open. Not against #616.

Frontier: `claude/matching-one-workspace-pwr5pv` @ `8b5f9d1a`.

---

## North-star

On an honest torus the physical measure at parameter `p` is product Bernoulli. That gives **one** law `F_L(·; p)` and one shape `Z_L(·; p)`. W4 claims the matching-odd number `M_L(1/2)` *governs* that shape.

Kill or promote W4 by answering:

```text
Among measures you are allowed to form from the same 2^N configs,
is Z a function of M(1/2)?
```

Allowed measures (closed list):

| id | measure | why it is fair |
|---|---|---|
| μ_p | product Bernoulli(p) | the physical family |
| μ_{1/2} | uniform on configs | the self-dual parameter |
| ν_β | exponential tilt of `r_b` at p=1/2: `∝ e^{β r_b}` | moves the odd sector without leaving the config space |
| ρ_t | exponential tilt of `n_black` at p=1/2 (i.e. vary p — this is μ_p again, check consistency) | sanity |
| π_even | uniform on `{r_b = r_w = 1}` | even sector only |
| π_odd± | uniform on `{r_b=2}` and on `{r_b=0}` | pure odd atoms |

`Z` is undefined if the inverse-CDF is not strictly increasing. In that case report the CDF’s jumps and skip `Z`; that is itself a verdict (a two-atom law has no interior shape).

W4 **dies** if two allowed measures have the same `M(1/2)` (or the same `E[r_b]`) and different `Z` on `{0.1,…,0.9}` beyond bisection error (`> 10^{-8}` after you say the tolerance).

W4 **survives L=3,4** if every pair of allowed measures with a common `M(1/2)` has the same `Z`. That is not a theorem for percolation; it is a finite-L fact #622 must then absorb.

---

## Directions

### D1 — Load the config table once

From the enumerator, a table of `2^N` rows is not required in git. Stream configs, accumulate:

```text
counts of (n_black, r_b, r_w)
```

at L=3 and L=4. Assert `r_b+r_w=2` on the stream (if this fails, stop the probe and report). Do **not** make “dual_fail=0” a results headline; it is a precondition.

### D2 — Physical `Z_L(u; p)` along p

Using the joint `(n_black, r_b)`:

```text
M(p) = E_p[r_b] - 1
F(p) = [1+M(p)]/2
Q(u) = F^{-1}(u)     (bisection, ≤1e-14)
Z(u; anchors 0.2/0.8, also 0.1/0.9)
```

Report, as **new** tables, not as a recap of `p_L^H`:

- `Z_L(u)` at the physical measure, for L=3 and L=4, at `p = p_L^H` and at `p=1/2` (these are different points because `M(1/2)≠0`);
- `‖Z_3 − Z_4‖_∞` on `{0.1,…,0.9}` at each of those p-choices.

Two sizes do not make a limit. The number `‖Z_3−Z_4‖_∞` is the first datum #622 has on W1 for *percolation*, as opposed to toys. If it is `O(10^{-3})`, constancy is not killed; if it is `O(10^{-1})`, constancy is dead at these sizes. Either is a result.

Also `Q(u)+Q(1-u)-1` at these p (single-model self-symmetry). If #619 already has this at generic p, only add the rows at `p=1/2` and `p=p_L^H`.

### D3 — The kill test for W4

At **fixed L=3** (512 configs; exact rational weights):

Compute `(M(1/2), Z)` — wait: `M(1/2)` is a number of μ_{1/2}. For tilted measures, the analogue is `E_ν[r_b]−1` at that measure, call it `m(ν)`.

Sweep `β` in the `r_b`-tilt `ν_β` (10–30 values). Plot/table `m(ν_β)` vs `Z(u)` at three u (0.2, 0.5, 0.8).

- If `Z` moves while `m` is held near `-21/64` (you will need a two-parameter tilt to hold `m` and move something else: tilt `n_black` and `r_b` together, or mix `π_even` with a small `π_odd` mass at a **fixed** odd expectation), that is the kill.
- Concrete kill construction: mix `π_even` with `π_odd+` and `π_odd−` at weights that keep `E[r_b]` fixed at the physical `E_{1/2}[r_b]`, and vary the even-sector internal law… but `π_even` is a single atom at L=3? **No**: `r_b=1` has many configs, with different `n_black`. Reweight *inside* `{r_b=1}` by `n_black` (exponential tilt of `n_black` conditional on `r_b=1`), **holding the odd masses P(r=0), P(r=2) fixed**. Then `m` is fixed (it depends only on the rank law) while the p-axis occupancy inside the even sector changes, which **can** move the Bernoulli-`p` inverse-CDF if you then evaluate Z of a **p-family built from those reweighted configs**.

Write this carefully. There are two different games; do both:

**Game A (rank law fixes M, occupancy still moves F).** Reweight configs at a **fixed p** (say p=1/2, so weights are just config weights, not `p^{n}(1-p)^{N-n}`). Then `F` is not the physical Bernoulli F; it is the CDF of a fake observable. That game is only about whether `Z` of a fake law is pinned by `m`. Useful but easy.

**Game B (physical Bernoulli family, configs reweighted by a tilt that is not n_black).** Replace product Bernoulli by `P(σ) ∝ p^{n}(1-p)^{N-n} e^{β r_b(σ)}`. Then `M_β(p)` and `Z_β` are both functions of `(p,β)`. Ask: at the `p` such that `F_β(p)=1/2` (the analogue of `p_L^H(β)`), does `Z_β` depend on `β`? If yes, W4 dies for the physical interpolation. If no, W4 survives this tilt family.

**Game B is the one that matters.** Game A is a page. Do Game B on L=3 fully (β-grid) and on L=4 at β=0 and two nonzero β.

### D4 — Self-dual comparison, square-bond L=3

Enumerate square-bond L=3 (`2^{18}`). Observable: ambient homology rank `r` as in #608’s exact L=3 path (or a 150-line self-contained rank). Physical p=1/2 is self-dual, `M(1/2)=0` if the observable is duality-odd in the usual way — **check**, do not assume.

Compute `Z` of this bond law at p=1/2. Compare `‖Z_bond − Z_site‖_∞` at L=3. If the two experiments are incomparable, say so (different edge sets, different `N`). If `Z_bond` is well-defined and far from `Z_site`, W2’s “same limit plus odd contamination” is not visible at L=3. Time box 20 CPU-min; no MC substitute.

Do **not** redo #619’s wrap/`X` identity. You may import it as a sanity check that you are on the same `r`.

### D5 — One-paragraph verdict for #622

```text
W4 at L=3,4:  KILLED / SURVIVES-THIS-FAMILY / ILL-POSED
evidence:     (Game B number: ΔZ at fixed m, or at p_L^H(β))
Z_3 vs Z_4:   ‖·‖_∞ = …
bond vs site: …
```

No ninth mechanism. No exponent.

---

## Stop rules

- D1 Alexander fails → stop.
- Game B on L=3 is enough to kill; if killed, D4 is optional.
- If Game B cannot be coded without a new percolation engine, stop and say `NOT_CODED`; do not invent MC.
- Do not fit `Z_3,Z_4` to `L^{-θ}`.

---

## Deliverables

```text
notes/probe-Mhalf-vs-shape-YYYYMMDD.md
scripts/probe/mhalf_vs_Z_tilt.py          (Game B)
scripts/probe/mhalf_Z_physical_tables.py  (D2; skip if #619 imported)
scripts/probe/bond_L3_Z.py                (D4, or skip with reason)
results/probe-Mhalf-vs-shape/latest.json
```

## Interface

#622 W4 lives or dies by Game B. #619 polynomials are a library. #608 `r` is the bond observable. #618/#620 not involved. Do not score productions.
