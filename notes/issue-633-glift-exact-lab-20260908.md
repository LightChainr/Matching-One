# #661 — Compatible g-lift vs exact-lab ΔZ (after the #633/#653 repairs)

Leftover of #633, not a new production, not #622. One comparison, one notes
file. No N725 rerun, no new Monte Carlo, no STATUS, draft PR, no merge.

**Objects used (all existing):**

| object | source |
|---|---|
| frozen g (9 deciles, 725 not in the freeze) | `results/type582-residual/latest.json` `consensus_g`, cross-checked identical to `results/p612-n725-score/latest.json` spin0 consensus |
| repaired exact laboratory | `results/probe-invariant-shape/census-exact.json` (PR #653; bond `dual_fail = 0`; PR #628's broken bond Z **not** copied) |
| PR #655 lineage | `results/probe-invariant-shape/n725-zflow-corrected.json`, `scripts/threshold_quantile_lineage.py` (brought here unchanged from `analysis/633-n725-spin0`) |
| withdrawn historical reading | `results/probe-invariant-shape/blindness-and-glift.json` (93.9° raw / 159.6° affine-removed — withdrawn in #633) |

**Engine:** `scripts/probe_invariant_shape/issue661_glift_exact_lab.py` →
`results/probe-invariant-shape/issue661-glift-exact-lab.json`.

---

## Formula and anchors

```text
DZ_Q[g](u) = (g(u)−g(a))/W − (Q(u)−Q(a))(g(b)−g(a))/W²
W = Q(b)−Q(a)
```

Shared anchors **a = 0.2, b = 0.8** (both exist on g's 9-decile grid; no
interpolation of g, so the 1/4–3/4 sensitivity is not even computed).
Angles are **unoriented** (|cosine|), per the #633 retraction of the
159.6° oriented reading.

**Chart exactness.** The census stores Q only at the 1/4, 3/4 anchors and
Z on the 9 deciles; Q on the deciles is recovered as `Q(u) = Q(1/4) + (Q(3/4)−Q(1/4))·Z(u)`,
i.e. exactly up to an affine map of Q. The lift direction is invariant
under `Q ↦ cQ + d` (because `DZ_{cQ+d}[g] = DZ_Q[g]/c`), so the angle
below is chart-exact, not an interpolation. Nothing about the anchor pair
0.2/0.8 is missing on either side, so no angle is left undetermined.

**Comparison object.** ΔZ = Z_{L=4} − Z_{L=3} on the same deciles,
recomputed from the repaired census — the same object named in #633
(`blindness-and-glift.json` direction 6). The exact laboratory is
deterministic (a census), so it carries no sampling covariance; the N725
jackknife covariance lives in `n725-zflow-corrected.json` and is cited,
not recomputed.

---

## Result: g-lift vs exact-lab ΔZ, unoriented

| Q used for the lift | unoriented angle | \|cosine\| |
|---|---:|---:|
| site L=3 (primary) | **24.87°** | 0.9073 |
| site L=4 | 24.72° | 0.9083 |
| bond L=3 (dual_fail = 0) | 25.81° | 0.9003 |

Primary reading (spin0 g, Q = site L=3, ΔZ = Z_{L=4}−Z_{L=3}):

```text
DZ_Q[g]:   [-0.0459,  0.0000,  0.2328,  0.4543,  0.6099,  0.6550,  0.5133,  0.0000, -1.5582]
ΔZ_exact:  [-0.0068,  0.0002, -0.0010, -0.0035, -0.0055, -0.0057, -0.0031,  0.0049,  0.0258]
unoriented angle 24.87°  (|cos| = 0.9073)
```

**What this means.** The g-lift with anchors 0.2/0.8 is a *chart-compatible*
version of the #633 direction-6 question: g is lifted into the affine chart
of Q before the angle is taken, instead of comparing raw 9-vectors. On the
exact laboratory the lifted g and the lab tangent are close to parallel
(|cos| ≈ 0.91, unoriented ≈ 24.9°) — a materially different picture from
the withdrawn raw 93.9°/159.6° reading, and from the raw tiny-torus
unoriented 20.4°.

**Scope limits, stated not hidden:**

* The exact-lab/tiny-torus ΔZ is a toy-size tangent. **It cannot prove an
  asymptotic mismatch**, and the ~25° agreement here does not prove a lift
  in the N-regime either.
* On N725 there is still no second production size, so N725 g-vs-ΔZ
  remains **undetermined** — unchanged from #655. The 81.4° number in
  `n725-zflow-corrected.json` is g-lift vs Z_725 in the 0.2/0.8 chart (a
  same-size parallelism check), **not** the g-vs-ΔZ test and not used as
  an incompatibility claim.
* The `equal` consensus g (a different freeze vector, not just a
  reweighting) gives 25.10° / 24.97° / 25.95° across the same three Q's —
  sensitivity only, not a second vote (Δ < 0.3° on the primary).

---

## What was NOT done

* No N725 rerun; `n725_100m.hist.csv` untouched.
* No new production, no new S theory (#622 territory).
* No STATUS edit, no merge, #633/#622 stay open.
* The PR-#628 broken bond Z (dual_fail = 118133) is asserted against: the
  script refuses to run on a census whose bond `dual_fail ≠ 0`.

Reproduce:

```bash
python3 scripts/probe_invariant_shape/issue661_glift_exact_lab.py
python3 -m unittest discover -s tests -p "test_issue661*"
```
