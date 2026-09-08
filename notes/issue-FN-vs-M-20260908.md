# Issue #678 — Is `F_N = (1+M_N)/2` an identity of the rank triple, or only a coding convention?

**Date:** 2026-09-08
**Assign:** local (ASSIGNED_MACHINE: local). Parent #650. Follow-ups: PR #675 Fact 1, PR #676 A′, #613's `F_N = E[r_G]/2`.
**Constraint compliance:** read committed JSON only, no new census, no L=5, no Huawei, no `docs/STATUS.md`, no ticket closures, #636 not funded. Integers/rationals only — every check below is `fractions.Fraction` arithmetic; no float appears in any exact artifact or intermediate.

**Machine-checkable reader:** `scripts/probe/fn_identity_reader.py` (added in this PR) re-runs checks A–D from the committed artifacts exactly; checks E–F additionally run when PR #668's `results/issue665-joint-rank-wrap/axis-L{3,4}.json` is present in the tree (it is not yet on main — those two checks were run this session against the PR #668 branch artifacts via `git show pr/668:...` and pass bit-for-bit).

---

## 0. Verdict (one page)

**`F_L = (1 + M_L)/2` is a theorem (an exact polynomial identity of the rank triple), not a coding convention — on the committed axis L=3,4 tables.** Concretely, coefficientwise over the power basis with `Fraction` arithmetic, both at L=3 (N=9) and L=4 (N=16):

1. `P_0 + P_1 + P_2 = 1` (the rank layers partition the Bernstein measure);
2. `M_L = P_2 − P_0` — reproduces the committed `M_poly_L{3,4}.json` power coefficients bit-for-bit;
3. `F_L := E[r]/2 = (P_1 + 2 P_2)/2` equals `(1 + M_L)/2` as polynomials;
4. boundaries `F_L(0)=0`, `F_L(1)=1`, and `F_L(1/2) = (1+M_L(1/2))/2` exactly: `F_3(1/2) = 43/128`, `F_4(1/2) = 19011/65536`.

The reason it *had* to be an identity: `P_0 + P_1 + P_2 = 1` and `E[r] = P_1 + 2P_2` give

```text
F = E[r]/2 = (P_1 + 2P_2)/2 = (P_0 + P_1 + P_2 + P_2 − P_0)/2 = (1 + M)/2,
```

with `M = P_2 − P_0` exactly the #606/root-ledger definition. So Fact 1 (PR #675) and the #612 chart identity are the same `D`-side object at finite L: `D = 1{r=2} − 1{r=0}` has expectation `P_2 − P_0 = M` (both `both-two` and the rank-2×rank-2 cells are empty in the committed tables — PR #668 `five_name_by_rank_black.both-two = [0,0,0]`, `rank_pair_totals` has only `(0,2),(1,1),(2,0)`), and `E[(1+D)/2] = E[r]/2` because `r ∈ {0,1,2}`. No defect in Fact 1's `D` vs the Bernstein `a_k`: `a_k = #cross_k − #none_k` per layer reproduces PR #668's `collapsed_D` and the committed Bernstein integers at both L.

**Corollary (Q3):** since `F_L` is strictly increasing (committed Sturm: `M'_L > 0` on (0,1), one root), `F_L(p*_L) = (1 + 0)/2 = 1/2` at the unique root `p*_L` of `M_L`, hence `p*_L = Q_L(1/2)`.

**Reconciliation with W4 (Q3, kept distinct):** `M_L(1/2) ≠ 0` at both L (`−21/64`, `−13757/32768` — committed, exact), so `p*_L = Q_L(1/2) ≠ 1/2`. This is consistent with, and says nothing new about, PR #671's W4 verdict: `M(1/2)` is one pointwise value of the odd polynomial and does not pin the shape `Z`. `p*_L = Q_L(1/2)` is a *location* statement inside the committed chart `F = (1+M)/2`; W4 is a *shape* statement. They do not interact, and this note does not revive "M(1/2) governs the shape".

---

## 1. Sources (read-only, committed artifacts)

| Artifact | What is taken |
|---|---|
| `results/probe-exact-controls/joints.json` (PR #653 lineage) | `by_L.{3,4}.joint_black_count_rank` — per-`(k, r_b)` configuration counts, 512 / 65 556 configs, `dual_fail_count = 0`, `max_abs_rank_diff = 0` |
| `results/probe-exact-controls/M_poly_L{3,4}.json` (PR #653 lineage, generator `scripts/probe/exact_controls_q_and_m.py`) | committed power-basis `M_coeffs`, `M_half`, `Mp_half`, `sturm_roots_in_01 = 1` |
| `results/issue665-joint-rank-wrap/axis-L{3,4}.json` (PR #668) | per-row `(k, r_black, r_white, wrap, label5, count)` joint rows, `five_name_by_rank_black`, `collapsed_D`, `committed_bernstein`, `rank_pair_totals` |
| `results/probe-exact-controls/latest.json` | `C1_census` rank pairs `(0,2)/(1,1)/(2,0) = 259/162/91` and `36559/19932/9045`; `C3_quantiles` `Q_L(0.5)` |
| `results/homological-balance-exact-torus/latest.json` | `M_half`, `p_L^H` cross-check |

No artifact was re-enumerated; every integer below is quoted from the above files.

---

## 2. The check (exact, what the reader script does)

From `joint_black_count_rank` (and independently from PR #668's `joint` rows), for each layer `k` let `c_{k,r}` be the number of configurations with `k` occupied sites and ambient rank `r_b = r`. Each configuration has Lebesgue weight `p^k (1−p)^{N−k}` under the product measure, so the rank masses are the polynomials

```text
P_r(p) = Σ_k c_{k,r} · p^k (1−p)^{N−k},      r ∈ {0,1,2},
```

expanded to the power basis over `Fraction`. Checks (all coefficientwise, all exact):

- **(A)** `P_0 + P_1 + P_2 = 1` — i.e. `Σ_r c_{k,r} = C(N,k)` per layer, verified per layer before expansion;
- **(B)** `P_2 − P_0` equals the committed `M_coeffs` of `M_poly_L{L}.json` coefficient for coefficient;
- **(C)** `(P_1 + 2P_2)/2` equals `(1 + M_L)/2` coefficient for coefficient (this *is* Q1);
- **(D)** `F_L(0) = 0`, `F_L(1) = 1`, `F_L(1/2) = (1 + M_L(1/2))/2` evaluated in `Fraction` (no floats);
- **(E)** per-layer `a_k = c_{k,2} − c_{k,0}` equals PR #668's `collapsed_D` and `committed_bernstein`.

Results:

| L | N | (A) | (B) | (C) | (D) | F_L(1/2) | (E) |
|---|---|---|---|---|---|---|---|
| 3 | 9 | ✓ | ✓ | ✓ | ✓ | `43/128` | ✓ (`[-1,-9,-36,-78,-90,-36,36,36,9,1]`) |
| 4 | 16 | ✓ | ✓ | ✓ | ✓ | `19011/65536` | ✓ (`[-1,-16,-120,-560,-1812,-4272,-7448,-9424,-7874,-2896,1720,2832,1660,560,120,16,1]`) |

`F_3(1/2) = 43/128` matches `notes/probe-exact-controls-new-map-20260907.md` §2.5; `F_4(1/2) = 19011/65536` matches `(1 + M_4(1/2))/2 = (1 − 13757/32768)/2`. The earlier `probe-invariant-shape-limit-20260907.md` float `0.29508…` for `F_4(1/2)` is **not** this number; the exact value committed here is `19011/65536 = 0.290084…` — that note's `0.29508` appears to be a transcription of a different quantity and should not be cited for `F_4(1/2)` (flagged, not edited — that note belongs to PR #676 lineage).

**Why the identity is structural, not numerological (Q1, the algebra):** per configuration `r_b + r_w = 2` (I1, zero dual-fail violations in the committed tables), so `r_b ∈ {0,1,2}` and

```text
D := 1{r_b=2} − 1{r_b=0} = r_b − 1      (pointwise, per configuration),
```

because both sides are functions of `r_b` alone agreeing on {0,1,2}. Taking expectations: `E[D] = M_L`, and `(1 + D)/2 = r_b/2` gives `F = E[r]/2 = (1+M)/2`. The pointwise coincidence `D = r − 1` is exactly where the *appearance* of convention ends and the theorem begins; PR #675's Fact 1 (`a_k = #cross_k − #none_k`) is its per-layer version and check (E) confirms it bit-for-bit.

---

## 3. The rank triple in five-cell language (Q2)

From PR #668's `five_name_by_rank_black` (configuration counts; per-rank triples `[rank0, rank1, rank2]`):

| cell | L=3 | L=4 | rank |
|---|---|---|---|
| `none` (`n×b`, no wrap) | 259 | 36 559 | 0 |
| `x` (one-axis, dir0) | 78 | 9 406 | 1 |
| `y` (one-axis, dir1) | 78 | 9 406 | 1 |
| `both-same` rank-1 (spirals) | 6 | 1 120 | 1 |
| `both-same` rank-2 (exclusive crosses) | 91 | 9 045 | 2 |
| `both-two` | 0 | 0 | — |

So, exactly:

```text
P_0 = mass(none),                              (the (0,2) exclusive-cross cell)
P_1 = mass(x) + mass(y) + mass(spiral),        (rank-1: both one-axis cells plus spirals)
P_2 = mass(exclusive cross),                   (the (2,0) cell inside both-same)
```

and this matches the PR #668 `rank_pair_totals` at both L: `x+y+spiral = 78+78+6 = 162 = #(1,1)` and `9406+9406+1120 = 19932 = #(1,1)`; `none = #(0,2)`; exclusive cross = `#(2,0)`. **`P_1` is confirmed to be the rank-1 mass: spirals plus the two one-axis cells** (the issue's guess is right, verified at both L). `both-two` is empty at both L; the `n×b` / `b×n` coarse pairing of #668 carries the same information since `both-same × both-same` is pure `(1,1)` (PR #676 A′), so `d0 + d1` / `b×b` refinements do not alter the rank projection.

---

## 4. Corollary and reconciliations (Q3)

- **`p*_L = Q_L(1/2)`.** `M'_L > 0` on (0,1) (committed Sturm, `sturm_roots_in_01 = 1` and `Mp_half > 0`), so `F_L` is a strictly increasing bijection [0,1]→[0,1] and `F_L^{-1}(1/2)` is the unique root of `M_L`. Cross-check: committed `C3_quantiles.Q_L(0.5)` equals `p_L^H` bit-for-bit at both L (`0.5865114551126757`, `0.5906721123310283`).
- **`Q_L(1/2) ≠ 1/2`.** `M_L(1/2) = −21/64` (L=3) and `−13757/32768` (L=4), both exact and nonzero, so the median quantile level does not sit at the self-complementing occupation point. This is the site non-self-duality already recorded in the root ledger (`M_L(1/2) ≠ 0`).
- **Not W4.** PR #671's W4 killed "M(1/2) pins shape `Z`". The corollary above is purely a location identity inside the chart `F = (1+M)/2`; it neither uses nor constrains `Z`. The two claims are kept distinct as demanded by the ticket.

---

## 5. Boundaries honored

- No new census, no re-enumeration: all integers read from committed JSON (`joints.json`, `M_poly_L{3,4}.json`, PR #668 `axis-L{3,4}.json`).
- No L=5 anywhere. No Huawei; `NEED_HUAWEI` not needed.
- No `docs/STATUS.md` edit; #613/#625/#635/#658 not closed; #636 not funded.
- Integers/rationals only: the reader script uses `fractions.Fraction` exclusively; the only floats in this note are the two labeled `Q(0.5)` cross-check strings quoted verbatim from committed JSON.

Related: #650, #613, #625, #635, #658, #672, #678; PRs #653, #668, #671, #675, #676.
