# #659 Integer structure of the five wrapping-type cells (2026-09-08)

Reads only the committed wrapping-type tables:

- `results/wrapping-type-census/axis-L3.json`, `axis-L4.json`, `diamond-L3.json` (PR #653, branch `repair/632-bond-lab-and-wrapping-l3l4`)
- `results/wrapping-type-census/axis-L5.json`, `diamond-L4-k1.json`, `diamond-L4-k2.json` (PR #657, branch `analysis/640-651-axis-L5-diamond-L4`)

No new enumeration. No L=6. Reader script: `scripts/wrapping_five_cell_reader.py` (reads the JSONs, re-verifies every identity in this note as exact integers).

Notation: `n×b` = #(neither×both), `b×n` = #(both×neither), `d0` = #(dir0×dir0), `d1` = #(dir1×dir1), `b×b` = #(both×both), all at black-mass `k`. `N = L²` (axis) or `2L²` (diamond).

---

## 1. Re-verified base facts (exact, all five tables)

- **Five-cell support.** All 16 cells vanish except `n×b, b×n, d0, d1, b×b`. Verified by scanning every `k` of every table.
- **Row sums:** `n×b + b×n + d0 + d1 + b×b = C(N,k)` at every `k` (all five tables).
- **Spiral pairing:** `d0 = d1` at every `k` (all five tables). Equivalent row-sum form: `n×b + b×n + 2·d0 + b×b = C(N,k)`.
- **MZ pairing:** `a_k = b×n(k) − n×b(k)` equals the committed Bernstein coefficient at every `k` on all five tables (axis L=3,4,5; diamond L=3,4).
- **No complement symmetry:** `n×b(k) ≠ b×n(N−k)` on every table (already noted in the #653 note; confirmed here cell-wise).

## 2. New exact identities (verified as integers on every committed table)

The five cells have sharp onsets. Everything below is checked bit-for-bit on all tables.

### (I1) Low-`k` binomial regime — `n×b(k) = C(N,k)` below the first wrap

- axis: `n×b(k) = C(N,k)` for `k < L` (L=3: k≤2; L=4: k≤3; L=5: k≤4)
- diamond: `n×b(k) = C(N,k)` for `k < 2L` (L=3: k≤5; L=4: k≤7)

Reading: below the smallest single-direction wrap (a full straight line: `L` cells on axis, `2L` on diamond), every configuration has black wrapping nothing and white wrapping both.

### (I2) High-`k` binomial regime — `b×n(k) = C(N,k)` above the blocking threshold

- **Both geometries:** `b×n(k) = C(N,k)` for `k ≥ N − L + 1`
  (axis: k≥7/13/21 for L=3/4/5; diamond: k≥16/29 for L=3/4), and strictly below at `k = N−L`.

Reading: a set fails to wrap dir0 only if it is contained in the complement of a blocking line of `L` cells; past size `N−L` every black set wraps both directions, hence white (its complement) wraps neither.

### (I3) Onset values

- axis `d0(L) = L` (3, 4, 5): the only dir0-wrapping `k=L` sets are the `L` full rows (columns are dir1).
- axis `b×n(2L−1) = L²` (9, 16, 25): the minimal both-wrap is a full row plus a full column (`2L−1` cells); `L²` crossing choices.
- diamond `d0(2L) = 4·C(2L,4)` (L=3: 60 = 4·15; L=4: 280 = 4·70) — clean but interpretation open.
- diamond `b×b(2L) = 2L` (6, 8).
- diamond `b×n` onset at `k = 3L−1` (L=3: k=8; L=4: k=11) with first value `4L²` (36, 64). The onset formula `3L−1` fits both diamond sizes; note axis onset is `2L−1`. Flagged: only two diamond points, treat as conjecture.

### (I4) Exact deficit decomposition below the both-wrap onset

- axis, for `L ≤ k < 2L−1`: `C(N,k) − n×b(k) = d0(k) + d1(k)` (no `b×b`, no `b×n` yet).
- diamond, for `2L ≤ k < 3L−1`: `C(N,k) − n×b(k) = d0(k) + d1(k) + b×b(k)`.

So in the pre-both-wrap window the collapsed deficit is carried by dir-paired cells alone — the L=3 sign-change analysis of PR #653 generalizes: every sign change of `a_k` in this window is two-term (`d0` vs `n×b`), not three-term.

### (I5) `b×b` onsets (raw, no closed form found)

- axis `b×b` starts at `k = 2L` with values 6, 120, 910 (L=3,4,5). No binomial-factor closed form fits all three (2·C(2L,·) patterns fail at L=4 or 5); recorded as gap.
- diamond `b×b` starts at `k = 2L` with value exactly `2L` (see I3).

## 3. Sequences: recurrences, OEIS, Catalan/Motzkin

Windowed sequences (support windows), per table:

- axis `d0`: L3 `[3,18,36,21]`; L4 `[4,48,280,992,2230,2976,2128,672,76]`; L5 `[5,100,1000,6500,30300,106185,285725,591150,923875,1054375,840585,448500,155450,33950,4325,255]`
- axis `b×b`: L3 `[6]`; L4 `[120,416,448,128,8]`; L5 `[910,9900,46800,123000,192700,179890,97350,29650,4800,350,10]`
- diamond `d0`: L3 `[60,612,2790,7038,10350,8496,3741,864,108,6]`; L4 `[280,6080,63104,416256,1955800,6939136,19184832,41955584,72799556,99548800,105751456,85723936,52177384,23625792,7955456,2000384,375004,51360,4880,288,8]`
- diamond `b×b`: L3 `[6,72,468,1464,2304,1584,378]`; L4 `[8,192,2528,22528,145392,697472,2516160,6861248,14122828,21708992,24420032,19560128,10804544,3964736,918848,122176,7192]`

OEIS queries (run 2026-09-08, exact integer search on oeis.org):

- `5,100,1000,6500,30300,106185,285725` — **no match**
- `910,9900,46800,123000,192700,179890,97350` — **no match**
- `4,48,280,992,2230,2976,2128,672,76` — **no match**
- `60,612,2790,7038,10350,8496,3741,864,108,6` — **no match**
- `120,416,448,128,8` — **no match**
- `280,6080,63104,416256,1955800,6939136` — **no match**
- `8,192,2528,22528,145392,697472` — **no match**
- `3,18,36,21` — **no match**

One near-hit worth recording: axis L=3 `n×b` for `k ≤ 5` coincides with row `n=3` of **A279445** (k-point placements on an n×n grid with ≤2 points per row/column): `1, 9, 36, 78, 90, 45`. This is a small-L coincidence — at L=3 "no full row/column" is the same as "no 3-in-line" — and the axis L=4 `n×b` row does **not** match A279445's row `n=4` (e.g. 1812 vs 528 at k=4). No Catalan, Motzkin, or nested-path sequence matches any of the five-cell windows.

No fixed-L, fixed-cell sequence tested is P-recursive of low order in the window either (windows are too short to certify a recurrence; recorded as gap, not a negative claim).

## 4. Generating functions / factorizations (Q3)

Converting the committed Bernstein coefficients to the power basis and factoring exactly over ℚ (SymPy `factor_list`, content ±1):

- axis L=3 (deg 9), axis L=4 (deg 16), axis L=5 (deg 25), diamond L=3 (deg 18), diamond L=4 (deg 32): **each M_L is irreducible over ℚ** — a single irreducible factor, no cyclotomic or binomial factors, no rational roots.

Consequence for proposed closed forms: any generating-function ansatz whose L-th diagonal specializations factor (e.g. products of low-degree terms) is dead at these sizes; M_L is already irreducible at L=3. The physical root drifts 0.58651 → 0.59067 → (axis L=5, root not recomputed here) 0.59…; diamond L=3 root 0.59425 — consistent with #638 and not re-litigated here.

## 5. What axis L=6 and diamond L=5 would have to show (Q4 — predictions only, no enumeration)

The identities above make sharp, integer-valued predictions. A single violated entry kills the corresponding claim; they are cheap to check against any future rung:

axis L=6 (`N=36`):
1. `n×b(k) = C(36,k)` for `k ≤ 5`, first deficit at `k=6` with `d0(6) = d1(6) = 6`;
2. `b×n(11) = 36` and `b×n(k)=0` for `k<11`;
3. `b×b` first mass at `k=12` (value unconstrained by this note);
4. `b×n(k) = C(36,k)` for `k ≥ 31`, strictly less at `k=30`;
5. `C(36,k) − n×b(k) = d0(k)+d1(k)` for `6 ≤ k ≤ 10`.

diamond L=5 (`N=50`):
1. `n×b(k) = C(50,k)` for `k ≤ 9`, first deficit at `k=10` with `d0(10) = 4·C(10,4) = 840` if I3 extends;
2. `b×b(10) = 10` if I3 extends;
3. `b×n` onset at `k = 14` with first value `100 = 4L²` if I3 extends;
4. `b×n(k) = C(50,k)` for `k ≥ 46`;
5. `C(50,k) − n×b(k) = d0(k)+d1(k)+b×b(k)` for `10 ≤ k ≤ 13`.

The `3L−1` diamond both-wrap onset and the `4·C(2L,4)` diamond onset value rest on two points (L=3,4) and are the most falsifiable items.

## 6. Method note / tripwire

Every identity in §1–§2 and the predictions in §5 are re-checked by `scripts/wrapping_five_cell_reader.py` against the committed JSONs; the script uses exact integer arithmetic only (no floats anywhere in the verification path). The factorizations in §4 were computed from the committed Bernstein integer lists via exact rational Bernstein→power-basis conversion.
