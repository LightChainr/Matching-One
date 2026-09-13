# #635 identity rewrite after #668 — `both-same` is not rank-2

**Ticket:** LightChainr/Matching-One#672 (parent #650; follow-up to PR #662 from
#660 and PR #668 from #665). **Date:** 2026-09-08.
**Mode:** notes-first rewrite. No new census. Axis L=5 joint census is
**NEED_HUAWEI** and was not run. Integers only. No `docs/STATUS.md` edit;
no ticket closed or merged; draft PR against main.

This note consumes, and does not re-enumerate:

- **PR #662** — wrapping-form A, written before the #668 dictionary
  (`notes/issue-635-identity-20260908.md`). The numerical identity there is
  correct; its dictionary is not. This note replaces that dictionary.
- **PR #668** — joint `(k, r_b, r_w, wrap_black, wrap_white, label5_black,
  label5_white)` census at axis L=3 (`2^9`) and L=4 (`2^16`)
  (`notes/issue-665-joint-rank-wrap-20260908.md`,
  `results/issue665-joint-rank-wrap/axis-L3.json`, `axis-L4.json`).
- **PR #646** — probe635 sector map, 5-name labels, two-amplitude decomposition
  (`notes/probe635-sector-map-20260907.md`,
  `results/probe635-sector-map/latest.json`).
- **PR #653** — #632 bond-lab repair; repaired rank lab; wrapping-type 4×4
  tables (`notes/wrapping-type-census-l3l4-20260908.md`,
  `results/probe-invariant-shape/census-exact.json`).
- **PR #657** — axis L=5 (`2^25`) and diamond L=4 (`2^32`) 5-name/coarse
  censuses; #651 A-continues
  (`notes/wrapping-type-census-axis-L5-diamond-L4-20260908.md`).

All exact claims below are integers / `Fraction`. Conjectures are marked.

---

## 0. What was wrong with the #662 dictionary

PR #662 wrote wrapping-form A as a difference of `both-same` sector
amplitudes and mapped `both-same` onto the coarse `both` class, treating
`both-same` as the rank-2 (cross) class. PR #668's joint census kills that
identification with exact integers (axis L=3, `2^9`):

```text
both-same (black side) by rank:  rank-1 (spiral) 6,  rank-2 (cross) 91
both-same (black side) at L=4:   rank-1 (spiral) 1120, rank-2 (cross) 9045
```

So `both-same = spirals ⊎ crosses` — the label is a *cluster-level* statement
(one cluster winds both axes) and a rank-1 spiral cluster satisfies it. PR
#662's identity `a_k = Δ#(both-same)` is nonetheless bit-for-bit correct at
every tested size:

```text
reference bernstein (axis L=3): [-1,-9,-36,-78,-90,-36,36,36,9,1]
  A_white (m-both-same):        [-1,-9,-36,-78,-90,-45,  0,  0, 0, 0]
  A_black (p-both-same):        [ 0,  0,   0,   0,   0,   9, 36, 36, 9, 1]
```

The rewrite below states the same identity in the joint language so a later
reader cannot confuse the sector label with the rank class, and proves the
correction is harmless — the spirals cancel, not because every `both-same`
event is a cross, but because D never fires on a `(1,1)` rank pair.

## 1. The joint language (#668)

For each configuration C of the L×L axis torus with |C| = k occupied sites,
record the six-tuple of PR #668:

```text
(r_b, r_w)   ambient H_1 winding ranks of the black NN (primal) graph and
             the white NN+NNN (matching) graph
(w_b, w_w)   coarse wrap labels ∈ {neither, dir0, dir1, both}
(s_b, s_w)   #646 5-names ∈ {none, x, y, both-same, both-two}
```

and the observable

```text
D(C) = 1{black NN wraps either axis} − 1{white NN+NNN wraps either axis}
     = (w_b ≠ neither) − (w_w ≠ neither)
M(p) = Σ_k a_k p^k (1−p)^(N−k),   a_k = Σ_{|C|=k} D(C).
```

Structural law (PR #653, certified against an independent incidence-matrix
elimination): **`r_b + r_w = 2` on every configuration**, zero violations at
L=3 (512 configs) and L=4 (65 536 configs). Consequence: rank-2 on both
sides is impossible at any L on site; at most one side ever carries rank 2.

## 2. Wrapping-form A′ — the identity in the joint language

**A′ (rewrite of PR #662's wrapping-form A).** At every finite L on the
axis torus, in the (rank pair × wrap × 5-name) joint language:

```text
a_k = Σ_{|C|=k} D(C)
    = Σ_{|C|=k} [ 1{(w_b, w_w) = (both, neither)} − 1{(w_b, w_w) = (neither, both)} ]
```

i.e. **the D-carrying cells are exactly the two exclusive-cross cells**, and
within each of them the 5-name pair and rank pair are constant:

| cell `(w_b × w_w)` | 5-name pair | rank pair | D | axis L=3 total | axis L=4 total |
|---|---|---|---|---|---|
| `neither × both` | `none × both-same` | `(0, 2)` | −1 | 259 | 36 559 |
| `both × neither` | `both-same × none` | `(2, 0)` | +1 | 91 | 9 045 |
| `dir0 × dir0` | `x × x` | `(1, 1)` | 0 | 78 | 9 406 |
| `dir1 × dir1` | `y × y` | `(1, 1)` | 0 | 78 | 9 406 |
| `both × both` | `both-same × both-same` | `(1, 1)` | 0 | 6 | 1 120 |

These five cells are the entire support at L=3 and L=4 (PR #668 §"Cell
support"); no sixth coarse cell is populated, and within each populated
coarse cell the 5-name pair is constant. The name A′ carries the two
corrections to #662's A:

1. **D-carrying cells are named as exclusive crosses.** `(neither × both)`
   is the white side cross-wrapping alone (white rank 2, black rank 0) — a
   single cross cluster on the matching lattice, by definition of the rank.
   `both-same` there is doing cluster-level work (one cluster, both axes);
   the *rank* statement `r_w = 2` is the cross statement. A′ writes D on the
   rank/exclusive-cross cells, not on the `both-same` label.
2. **`both-same × both-same` is not rank-2 and carries D = 0.** It is
   exactly `(r_b, r_w) = (1,1)` — 6 configs at L=3 (k=6), 1 120 at L=4
   (k=8: 120, k=9: 416, k=10: 448, k=11: 128, k=12: 8). Both sides wrap,
   neither side cross-wraps: these are the spirals.

**Per-k form at axis L=3** (PR #668 `axis-L3.json`, integers): D fires only
at `neither×both` (k=0:1, 1:9, 2:36, 3:78, 4:90, 5:45) and `both×neither`
(k=5:9, 6:36, 7:36, 8:9, 9:1), plus the D=0 spiral cell at k=6 (6) and the
D=0 `x×x`/`y×y` cells (k=3:3+3, k=4:18+18, k=5:36+36, k=6:21+21). Collapsing
over the coarse 4×4 with the D above reproduces the committed Bernstein
integers `[-1,-9,-36,-78,-90,-36,36,36,9,1]` bit-for-bit — the same
integers PR #662 accepted, now with the correct cell bookkeeping. The
L=4 joint table (`axis-L4.json`) reproduces the committed L=4 rung the same
way (`collapsed_D` field, bit-for-bit with `[-1,-16,-120,-560,-1812,-4272,
-7448,-9424,-7874,-2896,1720,2832,1660,560,120,16,1]`).

## 3. The spiral-cancellation lemma

**Lemma (spiral cancellation; exact at L=3,4; conjecture for all L).** For
every configuration C with `(s_b, s_w) = (both-same, both-same)` —
equivalently (at every L where the census has been run, and structurally at
all L) `(r_b, r_w) = (1,1)` — the contribution of C to `a_k` is zero:
`D(C) = 0`. Hence the rank-1 spiral mass inside `both-same` on either side
contributes nothing to `Δ#(both-same)`, and

```text
a_k = Δ#(both-same)(k)  still holds,  a_k = Δ#(exclusive cross)(k)  is the
rank-honest form,  and the two agree because the spirals cancel.
```

**Proof at L=3 and L=4 (from the #668 integers).** By definition
`D(C) = (w_b ≠ neither) − (w_w ≠ neither)`. On the cell `both × both` both
indicators are 1, so `D = 0` — this is arithmetic, not a census fact. The
census facts are: (i) the coarse `both × both` cell contains exactly the 6
(L=3) / 1 120 (L=4) configurations tabulated, all `(1,1)`, all
`both-same × both-same`; (ii) no `both-same` configuration occurs outside
the cells `both × neither`, `neither × both`, `both × both`; (iii) the
`x`, `y` cells (`dir0 × dir0`, `dir1 × dir1`) are diagonal and carry D = 0
the same way. Summing D over the full per-k joint table therefore equals
summing D over the two exclusive-cross cells alone, and the collapse
reproduces the committed Bernstein integers bit-for-bit at both sizes
(tripwire 1 of PR #668). ∎ (at L=3,4)

**Conjecture 3.1 (all L).** On the axis torus at every L: the support of
the joint distribution is contained in the five cells above; every
`(both, both)` configuration has `(r_b, r_w) = (1,1)`; every `(neither,
both)` configuration has `(0,2)` and every `(both, neither)` has `(2,0)`;
and `both-two` has zero mass at every k. Evidence: exact at L=3 (2⁹), L=4
(2¹⁶); the 5-name/coarse projections of the same statement hold at axis L=5
(2²⁵) and diamond L=4 (2³²) (PR #657: five cells, `both-two` mass 0, no
sixth cell) — but the joint rank-splitting of PR #668 has **not** been run
at L=5 (axis L=5 joint census is NEED_HUAWEI).

**Conjecture 3.2 (structural half of 3.1, provable-looking).** The rank
statements of 3.1 do not need enumeration to fail: `r_b + r_w = 2` (proved
at L=3 over all 2¹⁸ bond configurations in PR #646's repaired lab and
certified at site L=3,4 in PR #653/#668) forces at most one rank-2 side,
so rank-2 × rank-2 never occurs at any L, and if the cell support of 3.1
holds then the rank exclusivity `(0,2)`/`(2,0)` on the D-cells is forced.
What enumeration is actually needed for is the *wrap-cell* support claim
(no sixth coarse cell, `both × both` always `(1,1)`), which is genuinely
per-L. Marked conjecture because `r_b + r_w = 2` itself is certified at
finite L only (it is a consequence of exact bond duality on the L=3 square
bond torus, PR #646; its all-L status is not claimed here).

## 4. A-continues survives the rewrite

**Question 3 of the ticket: did A-continues (#651: axis L=5 / diamond L=4,
`both-two = 0`) secretly depend on the wrong dictionary? No.**

A-continues is a statement about the 5-name projection alone: the `both-two`
cell has zero mass at every k at axis L=5 (2²⁵) and diamond L=4 (2³²), and
the coarse 4×4 support stays the five Mertens–Ziff cells (PR #657). None of
its inputs is the identification `both-same = rank-2`:

- The five-cell support claim lives at the level of `(w_b, w_w)` pairs and
  5-name pairs, which PR #668 showed are constant within each populated
  coarse cell at L=3,4 — the projection discards exactly the rank column
  that the wrong dictionary misidentified, so the projection is insensitive
  to the misidentification.
- The wrong dictionary entered only in PR #662's *interpretation* of the
  `both × both` cell (it called the spiral pairing a cross-pairing
  implicitly, by treating `both-same` as the rank-2 class). The verdict
  letter (A at finite L, wrapping-form) and the A-continues caveat-upgrade
  never used that reading: they used D, which is defined on the coarse wrap
  labels, and the bit-for-bit Bernstein collapse, which is rank-blind.
- Under A′, A-continues restates cleanly: at axis L=5 and diamond L=4 the
  five cells of §2's table keep their D signs (the two exclusive-cross
  cells carry ±1, the three diagonal cells carry 0), `both-two` never
  appears to add a sixth row, and the collapse still reproduces the
  committed Bernstein rungs (PR #657 tripwire vs PR #649, pass on both
  diamond L=4 kernels).

So A-continues survives verbatim; what changes is only that its "five
Mertens–Ziff cells" language should be read with §2's table — the
`both × both` row is the spiral row, `(1,1)`, D = 0.

**What A-continues does not settle** (unchanged from PR #662 §4, and still
not reachable without Huawei): the joint rank-splitting at axis L=5
(NEED_HUAWEI — comment and stop if a future ticket needs that rung), and
the Jacobsen-map gaps (order of limits; one pTL operator vs two lattices;
novelty boundary vs Mertens–Ziff 2016 / Jacobsen 2015). None of those is
affected by the dictionary rewrite.

## 5. Boundary: this rewrite still does not fund #636

PR #662's §5 recommendation stands unchanged, and the rewrite strengthens
it slightly: the sector amplitudes of PR #646 are **not rank-pure** (PR
#668), so any future sector-level transfer-matrix reading would first have
to split `both-same` by rank to be measuring one homology class — an
additional cost the old dictionary hid. The wrapping-form of verdict A
remains published (Mertens–Ziff 2016); the unpublished piece remains the
conceptual sector identification; degeneracy, if exact (Conjecture 3.1),
still makes sector machinery carry less structure than its cost. Nothing
here promotes #636.

**Tripwire standing rule (unchanged):** any claimed decomposition must
reproduce `[-1,-9,-36,-78,-90,-36,36,36,9,1]` at axis L=3 bit-for-bit,
integers / Fraction only, or it is unused.

---

**Claim boundary.** Integers only in exact artifacts; the per-k and
per-cell integers above are quoted from `results/issue665-joint-rank-wrap/`
(PR #668, axis L=3,4), `results/probe-invariant-shape/census-exact.json`
(PR #653), and the PR #657 committed tables — none re-enumerated here; the
only new computation on the Mac for this note was none (notes-first; the 2⁹
occupancy lists were not needed, the #668 JSON carries the per-k cells
already). Conjectures 3.1/3.2 are marked. Do not cite #628's 118133 bond
`dual_fail` as physics (repaired compound artifact, PR #653/#646).
`M(p) + M(1−p) ≠ 0` on the site side stands as in PR #662.

Related: #635, #636 (not funded), #650, #651, #660, #665; PRs #645, #646,
#649, #653, #654, #657, #662, #668.
