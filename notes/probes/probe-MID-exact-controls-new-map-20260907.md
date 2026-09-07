# Probe MID — Exact algebraic controls the new map stands on

**Assign to:** one Agent with ordinary (but adequate) reasoning and ordinary CPU. Python 3 stdlib + numpy. L=4 site (65536) is already in-repo and must be reproduced. L=3 bond (`2^{18} = 262144`) is in budget. **No production, no arXiv, no 9-vector residual fitting, no N=725.**

**Literature is out of scope.** Cite in-repo notes only.

---

## Standing

This is a closed-checklist *exact* probe. The new map currently leans on four short facts. They must become re-runnable identities with JSON, tests, and failure counts, so the HIGH and rate probes can cite numbers instead of folklore.

- Does **not** enter `docs/STATUS.md`.
- Does **not** close #276, #321, #606, #608, #581, #615.
- Negative answers (an identity fails) are the most useful outcome.
- Cite vs claim: numbers you recompute are claims; numbers you only read are cites.
- If C1 (census reproduction) fails, **stop the whole probe**.

Prefer a PR against **frontier** (`claude/matching-one-workspace-pwr5pv` @ `8b5f9d1a`). If a file is missing, stop and name the path. Do not invent Monte Carlo substitutes.

You will need:

```text
PR #606   scripts/homological_balance/exact_torus_enum.py
          scripts/homological_balance/verify_independent.py
          results/homological-balance-exact-torus/latest.json
          tests/test_homological_balance_exact_torus.py
PR #608   scripts/score_qtangent_empirical.py and any L=3 exact path
PR #615   notes/bounded-task-rank-threshold-no-go-20260907.md   (cite only)
PR #614   notes/p613-quantile-convergence-20260907.md           (cite only)
```

---

## Why this exists

Load-bearing supports of the new map:

```text
(1) M_L(1/2) ≠ 0 on honest square-site tori
    L=3: −21/64     L=4: −13757/32768
    ⇒ NN square-site is not self-matching.

(2) Q_G(u) + Q_G(1−u) − 1  is not ~ 0.
    The matching involution is a two-model relation, not a
    single-model even/odd that fitted g could be an irrep of.

(3) p_L^H is the root of M. Under the contract F=(1+M)/2 and
    strict increase, Q_L(1/2)=p_L^H. The scientifically live
    statement is p_L^H ≠ 1/2, i.e. M(1/2)≠0, not “not the median”.

(4) On self-dual square bond, X=r−1, O=wrap_either:
        O = 1 + (X − X^2)/2     (configuration-wise, or kill it)
        Cov(O,X) = P(r=2) = Var(X)/2
    #608’s stable coupling is partly algebra + Harris, not a
    new mechanism.
```

North-star: **make (1)–(4) code, then push every exact consequence that still fits in ordinary CPU.** Do not generalise to observer groupoids or rates.

---

## Programs (all required unless a stop rule fires)

### C1 — Reproduce the L=3,4 site census

```bash
python3 -m unittest tests.test_homological_balance_exact_torus -v
python3 scripts/homological_balance/verify_independent.py
python3 scripts/homological_balance/exact_torus_enum.py
```

Must match:

| L | configs | dual fail | (0,2)/(1,1)/(2,0) | `M(1/2)` | `p_L^H` |
|---|---:|---:|---|---|---|
| 3 | 512 | 0 | 259 / 162 / 91 | `−21/64` | 0.586511455113 |
| 4 | 65536 | 0 | 36559 / 19932 / 9045 | `−13757/32768` | 0.590672112331 |

Two winding algorithms, mismatch = 0. If dest paths write to `scripts/results/`, fix `parents[2]` and say so.

**Stop the probe if any row breaks.**

### C2 — `M_L` as an exact polynomial

`M_L(p) = ∑_k (2 r_b(k)/2 − 1) p^{n_k} (1−p)^{N−n_k}` over configs, equivalently from rank-pair counts *and* from black-count × rank joints if you have them.

Deliver, L=3 and L=4:

- `M` as a polynomial with **exact integer / rational coefficients** (degree N=9 and N=16). If the enumerator only stored rank pairs, extend it to accumulate `#{configs : (n_black, r_b)}` — 512 and 65536 are small.
- exact `M(0), M(1), M(1/2)` as fractions;
- exact derivative `M'(p)` as a polynomial, plus `M'(1/2)` as a fraction;
- a uniqueness check: one root in `(0,1)` by Sturm or by exact Descartes + bisection with rational bounds.

This is what the rate probe’s R3 expansion actually sits on. Without the polynomial, “`M'(p_c) ≍ L^{3/4}`” has no finite-L left-hand side you can even *write*.

JSON: coefficients as lists of `{k, c}` with `c` a string fraction.

### C3 — Quantile table, and the self-symmetry killer

Using exact `M` from C2 (bisection on a rational polynomial, tolerance `≤ 1e-14` stated):

For L=3 and L=4, `u ∈ {0.1, 0.2, …, 0.9}`:

```text
Q_L(u)
Q_L(u) + Q_L(1−u) − 1          # must not be ~ 0
Q_L(0.5)  versus  p_L^H        # should agree under C2’s uniqueness
F_L(1/2) = [1 + M(1/2)]/2      # exact fraction
```

If `Q(u)+Q(1-u)-1` is consistent with 0 at all listed u within `10^{-12}`, the self-symmetry kill is **false** and you stop the map-support claims (still deliver the table).

Also report `Q_L(u) − p_c` using the same truncated `p_c = 0.5927460` the #606 JSON already uses, **as a cite of a convention**, not as a theorem.

### C4 — Joint `(n_black, r_b, r_w)` and Alexander

For every config, `(r_b, r_w)` with black 4-connect, white 8-connect, as the enumerator already does. Add:

- full joint histogram `(n_black, r_b)` and `(n_black, r_w)`;
- confirmation `r_b + r_w = 2` on every config (C1 already; keep the failure list if any);
- the six hand configs in `verify_independent.py` plus three more you declare *before* running (e.g. 2×2 block, checkerboard, full row plus full column). If a new hand config disagrees with intuition, the intuition loses.

### C5 — Two algorithms, every rank, not just mismatch=0

`verify_independent.py` reports mismatch counts. Extend:

- per-rank confusion table L=3 (512) and L=4 (65536);
- max `|r_A − r_B|` (must be 0);
- runtime, so L=4 staying in budget is documented.

If you find a single mismatch, dump the configuration and stop.

### C6 — Square-bond L=3, wrap / `X` configuration-wise

Self-dual **square-bond** torus, `p` is not needed: enumerate all `2^{2L^2} = 2^{18}` bond configs at L=3. **Do not enumerate L=4 bonds (`2^{32}`).**

Let `r ∈ {0,1,2}` be ambient homology rank, `X = r−1`, `O = wrap_either` (the same observable #608 uses; if the name is ambiguous in code, pin it to the function in `score_qtangent_empirical.py` or the exact L=3 gate).

Check configuration-wise:

```text
O  ?=  1 + (X − X^2)/2
```

Report failure count. 0 = identity is load-bearing. `>0` = #608 must not cite it as exact; dump ≤ 20 failing configs.

Consequences, from the identity plus duality `X ↦ −X` if you use it (say so):

```text
Cov_p(O,X) at p=1/2  ?=  P(r=2)  ?=  Var(X)/2
```

At p=1/2 every config is equally weighted, so these are exact fractions from counts. Compute them from the same enumeration, *and* from the identity (should agree).

Also exact, at p=1/2, L=3 bond:

- `P(r=0), P(r=1), P(r=2)`;
- `E[O], Var(O), Cov(O,X), Cov(O, X^2)`;
- if `B_even` is defined in the L=3 exact gate, the triple joint `(O, X, B_even)` as a 3×3×(range) table. This is the typed-control skeleton #608 actually needs. Do not interpret it as a mechanism.

Time box: if L=3 bond will exceed ~20 CPU-minutes, stop, report partial, **do not switch to MC**.

### C7 — Harris / nested-reveal on L=3 bond, exact

#608’s “single-signed scale ladder” was reduced to Harris association for increasing observables. On L=3 bond, with the nested reveal filtration the empirical scorer uses (or the exact `_conditional_means` path, which #608 says *does* run at L=3):

- confirm telescoping `∑ Γ_j = Cov(Y)` to `1e-16` (cite #608 if you reproduce, claim if you recompute);
- for `(O, X)`, check all cross-increments `Γ_j^{O,X} ≥ −ε` with `ε=1e-15`. If a negative increment exists, Harris is not the whole story and you dump it.

Do not run L=4/8 MC. Cite #608’s empirical table if present; do not treat it as exact.

### C8 — Matching involution, three-line lemma + one numerical check

Lemma (hypotheses: planar matching pair, honest torus, H1):

```text
M_Ĝ(p) = − M_G(1-p)
F_Ĝ(p) = 1 − F_G(1-p)
Q_Ĝ(u) = 1 − Q_G(1-u)
```

Negative: this does **not** imply `Q_G(u)+Q_G(1-u)=1`. C3 is the numerical support.

If the enumerator can flip to matching-adjacency (NN+NNN / 8-connect as “black”) on the same L=3 site configs, compute `M_match(p)` and check `M_match(p) + M_site(1-p) = 0` as a polynomial identity (up to the convention of which side is Ĝ). If the code cannot flip adjacency without a rewrite larger than ~100 lines, **do not rewrite a percolation engine**; record `NOT_CODED` and keep the lemma. The site 4-connect / white 8-connect already in the enumerator *is* this pair — check whether `M` computed from `r_w` at `1-p` equals `−M` from `r_b` at `p`. That check is in budget:

```text
M_black(p) + M_white(1-p)  ?=  0     as polynomials, L=3 and L=4
```

because `r_b + r_w = 2` implies `E[r_b(p)] + E[r_w(p)] = 2`, and white at p is black-matching at 1-p under the usual coupling. Write the coupling you actually use. If the identity fails, Alexander and the involution have been conflated and you stop.

### C9 — Protocol card, ≤ 40 lines, no scoring

Paste-ready card for later scorers. You do **not** apply it to N=725.

```text
declared chart     = span{1, Q_base, g_frozen}  unless HIGH probe overwrites
transport          = #612 identity before any residual
weightings         = spin0 AND equal, always together
forbidden          = second exponent; chart change after residual
Q_N(u) → p_c       = location only (H1–H3); not a rate
M(1/2) ≠ 0         = square-site is not self-matching   (C1/C2/C3)
Q(u)+Q(1-u) ≠ 1    = single-model self-symmetry is false (C3)
#608 Cov(O,X)      = algebra + Harris, not a mechanism   (C6/C7)
p_L^H = Q(1/2) ≠ 1/2   under F=(1+M)/2 and strict increase (C2/C3)
```

### C10 — Tests in the tree

`tests/test_probe_exact_controls.py`:

- C1 numbers;
- C2 `M(0), M(1), M(1/2)` fractions;
- C3 `Q(0.5) = p_L^H` within `1e-12`, and `max_u |Q(u)+Q(1-u)-1| > 10^{-3}` (the kill is that it is *not* zero);
- C5 mismatch = 0;
- C6 failure count = 0, or skip if C6 was time-boxed;
- C8 polynomial involution if coded.

Keep L=4 tests as slow-marked if needed, but C1 L=3 must be fast.

---

## Out of scope

- Affine holonomy, 4%, Fieller, `Z_N` on productions → HIGH.
- RSW, arm exponents, F1/F2 *proofs*, arXiv → literature probe.
- Conditional expansions F1+F2 ⇒ `L^{-4}` → rate probe (they will cite your polynomials).
- P398 widths 9–10, Gate 3 labels, STATUS, closing tickets, any optimizer.
- L=5 site (`2^{25}`) — not requested. Do not “just run it”.

---

## Stop rules

- C1 fails → stop everything.
- C3 shows self-symmetry (Q-sum-1 ≈ 0) → stop map-support, still deliver tables.
- C6 failures → do not cite wrap/`X` as exact; still deliver the count.
- C6 > 20 CPU-min → stop C6, do not MC.
- C8 identity fails → stop and report; do not “fix” Alexander by changing connectivity conventions without saying so.

---

## Deliverables

```text
notes/probe-exact-controls-new-map-YYYYMMDD.md
scripts/probe/exact_controls_q_and_m.py          (C2, C3; imports #606 enumerator)
scripts/probe/exact_controls_joints.py           (C4, C5)
scripts/probe/square_bond_wrap_identity.py       (C6, C7; L=3 exact)
scripts/probe/matching_involution_poly.py        (C8, or a stub NOT_CODED)
tests/test_probe_exact_controls.py
results/probe-exact-controls/latest.json
results/probe-exact-controls/M_poly_L3.json
results/probe-exact-controls/M_poly_L4.json
results/probe-exact-controls/bond_L3_wrap.json
```

Do not duplicate `exact_torus_enum.py`. Import it.

## Interface

- HIGH cites C3 (`Q+Q−1`, `M(1/2)`) as the self-symmetry kill and C6 as #608 algebra.
- Rate probe cites C2 polynomials and C3 `Q(0.5)=p_L^H ≠ 1/2`.
- Literature probe should not need you; if they do, they are in the wrong file.
