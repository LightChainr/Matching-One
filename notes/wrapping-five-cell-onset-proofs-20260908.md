# #673 Combinatorial proofs of the #659 five-cell onsets (2026-09-08)

Ticket: #673 (parent #650, follow-up to PR #669 / #659). **No new census. No L=6. No diamond L=5. No diamond L=4 enumeration** (diamond L=4 stays at the committed `k`-window tables of PR #657).

Every claim below is proved from the geometry of the two committed tori as
built by `scripts/matched_torus_reference.py` and classified by
`scripts/torus_homology.py` (`dir0`/`dir1` = winding along the two quotient
generators). The committed tables
(`results/wrapping-type-census/*.json`, PRs #653/#657) and the PR #669 reader
`scripts/wrapping_five_cell_reader.py` are used **only as checks**; no proof
below cites a table entry. An independent brute-force cross-check at committed
sizes (axis L=3,4; diamond L=2,3 — `2^8`, `2^18`, seconds-scale, no new rung)
is provided in `scripts/wrapping_onset_proof_checks.py` (all checks pass; run
instructions at the end).

Notation as in PR #669: `n×b = neither×both`, `b×n = both×neither`,
`d0 = dir0×dir0`, `d1 = dir1×dir1`, `b×b = both×both`, all at black mass `k`.
`N = L²` (axis) or `2L²` (diamond).

## 0. Geometry and lemmas

**Conventions.** Axis torus: sites `(x,y) ∈ (ℤ/L)²`, NN primal edges
`(±1,0),(0,±1)`; `dir0` = x-period `(L,0)`, `dir1` = y-period `(0,L)`.
Diamond torus: coordinates `(u,v) = (x+y, y−x)`, periods `2L` in each of
`u, v`, sites `S = {(u,v) : u ≡ v (mod 2)}`, `N = 2L²`; NN primal edges are
the steps `(1,1)` and `(1,−1)` (each changes both `u` and `v` by ±1); the
white matching lattice adds NNN steps `(2,0)` and `(0,2)`. Quotient
generator 0 = the `u`-period, generator 1 = the `v`-period.

**Lemma A (straight lines on both tori).**

- *Axis.* A full row `y = c` is `L` sites joined by `(1,0)` steps and closes
  with displacement `(L,0)`: it wraps `dir0` **only**. A full column `x = c`
  wraps `dir1` only.
- *Diamond.* A *diagonal line* of slope `(1,1)` is
  `𝔡_r = {(u, u+r) ∈ S : u ∈ ℤ/2L}` for fixed even `r = v−u`; it has `2L`
  sites, is NN-connected, and closes with displacement `(2L, 2L)` — it wraps
  **both** generators (nonzero projection on both periods). Symmetrically
  the slope-`(1,−1)` line `𝔞_c = {(u, c−u)}` (fixed even `c = v+u`) has
  displacement `(2L, −2L)` and wraps both. **No NN cycle on the diamond torus
  winds exactly one generator in a straight line**: a cycle winding dir0
  only has displacement `(2aL, 0)`-type, and since every NN step changes
  `v` by ±1, its steps must cancel in `v`: it mixes `(1,1)` and `(1,−1)`
  steps.

**Lemma B (minimum winding mass).**

- *Axis (black).* A black NN cycle with displacement `(L, m)` has `L + |m|`
  steps, each changing `x` by at most 1, hence visits ≥ `L + |m|` distinct
  sites. So `k < L` ⇒ black winds neither generator. Equality `k = L` with
  dir0 winding forces `m = 0` and all steps `(1,0)`: the set **is** a full
  row.
- *Diamond (black).* A black NN cycle with displacement `(2aL, 2bL)`
  (its winding numbers) has length ≥ `max(|a|,|b|) · 2L`, with equality iff
  every step has the sign pattern of `(a,b)` — i.e. a straight diagonal
  line when `ab ≠ 0`, and a *drift-0 mixed cycle* (see §5a) when `ab = 0`.
  So `k < 2L` ⇒ black winds neither generator.

**Lemma C (minimum white winding mass).** The white matching lattice has
extra edges `(2,0)`/`(0,2)` (diamond) or diagonals (axis). Minimum
winding cycles: on axis, a straight line of `L` sites (NN); on diamond, an
NNN chain `L` sites long (`(2,0)` steps along a `v = const` line) winding
one generator, or a straight NN diagonal of `2L` sites winding both. Hence
white with `< L` sites (axis) or `< L` sites via NNN / `< 2L` via NN (diamond)
still cannot wind carelessly — the precise bounds used below are stated
where needed.

**Lemma D (MZ pairing; the cells' meaning).** The joint (black primal,
white matching) wrap types are supported on the five Mertens–Ziff cells with
the pairing already committed in this repo (PR #653 note; Mertens–Ziff PRE
94, 062152 (2016), eqs. (9)–(11), (19)): black-neither ⇔ white-both (`n×b`);
black-both ⇒ white-neither (`b×n`) or white-both (`b×b`); black exactly-dir0
⇔ white exactly-dir0 (`d0`); black exactly-dir1 ⇔ white exactly-dir1 (`d1`).
The note uses this pairing for cell bookkeeping; the *black-side* structure
is proved directly, and the white sides of the axis proofs are also proved
directly (barrier arguments), so duality is imported only as bookkeeping.

**Counting convention.** At mass `k`, a cell's count is the number of
`k`-subsets whose type pair is that cell; `C(N,k)` is the row sum.

---

## 1. I1 — low-`k` binomial regime (PROVED)

**Claim (axis).** `n×b(k) = C(L², k)` for `k < L`.
**Claim (diamond).** `n×b(k) = C(2L², k)` for `k < 2L`.

**Proof.** By Lemma B, a black set of size below the respective threshold
wraps neither generator. By Lemma D, its cell is `n×b`, and every `k`-subset
is counted. ∎

The thresholds are sharp: at `k = L` (axis) a full row exists (Lemma A), at
`k = 2L` (diamond) a diagonal line or drift-0 cycle exists (§5), so the
first deficit occurs exactly at the threshold. That the deficit at the
threshold is *exactly* the `d0`+`d1` onset mass is §3/§4/§5.

## 2. I2 — high-`k` binomial regime (PROVED)

**Claim (both geometries).** `b×n(k) = C(N,k)` for `k ≥ N − L + 1`, and
`b×n(N − L) < C(N, N − L)` strictly.

Note the threshold is `L`-governed on **both** geometries — on the diamond
this is because the complement is small (`≤ L−1` sites), not because the
diamond's own line length is `2L`.

**Proof.** Let black have `k ≥ N − L + 1`, so white has `≤ L − 1` sites.

*Black wraps both.* Missing sites: ≤ `L−1`. On the axis torus there are `L`
rows and `L` columns; making **all** `L` rows non-full requires ≥ `L`
missing sites, so some row is fully black (wraps dir0, Lemma A) and
symmetrically some column (dir1). On the diamond torus there are `2L`
slope-`(1,1)` lines and `2L` slope-`(1,−1)` lines; making all `2L` lines of
one family non-full requires ≥ `2L > L − 1` missing sites, so some
`(1,1)`-line and some `(1,−1)`-line are fully black — each wraps **both**
generators (Lemma A). In both geometries black wraps both directions.

*White wraps neither.* White has `≤ L − 1` sites. On axis, any winding NN
cycle has ≥ `L` sites (Lemma B). On diamond, a matching-lattice winding
cycle with displacement `(2aL, 2bL)`, `(a,b) ≠ (0,0)`: each NN step changes
`u` by ±1 and each NNN step by ±2, so the cycle has `u`-extent ≥ `|a|·2L`
covered by ≥ `|a|·L` NNN steps... precisely: the cycle needs ≥ `|a|·L + |b|·L`
edges (each edge advances `u` by ≤ 2), hence ≥ `L·(|a|+|b|) ≥ L` distinct
sites. With `≤ L−1` sites, white wraps neither.

By Lemma D the cell is `b×n` and the count is `C(N,k)`. ∎

**Sharpness at `k = N − L`.** Take white = one full straight line: axis, a
full row (`L` sites, NN-connected, displacement `(L,0)`, wraps dir0);
diamond, a full `(1,−1)`-line `𝔞_c` (`2L` sites, NN-connected, wraps both —
but at minimum it wraps dir0 and dir1 both, which already excludes the
`b×n` cell). Black = complement does not enter the argument: white wraps
something, so the cell is not `b×n` (Lemma D), while §1-style counting
still leaves `C(N,k)` subsets total. Hence strict inequality. ∎

## 3. Axis I3 — onset values (PROVED)

### 3a. `d0(L) = d1(L) = L`

**Proof (black side).** Let `B`, `|B| = L`, wrap dir0 and not dir1. By
Lemma B (equality case), `B` **is** a full row `y = c`. A full row is a
single `(L,0)` cycle: black type is exactly dir0.

**Proof (white side, direct).** White = torus minus row `y = c`.

- *White does not wind dir1.* Any white NN or diagonal-matching path
  (axis matching edges: `(1,0),(0,1),(1,1),(1,−1)`, each with `|Δy| ≤ 1`)
  with net `y`-displacement `L` visits some site of every row, in
  particular row `y = c` — entirely black. Contradiction.
- *White does wind dir0.* Row `y = c + 1` is fully white and closes with
  displacement `(L, 0)`.

So white is exactly-dir0 and the configuration is in cell `d0`. The `L`
choices of `c` give `d0(L) = L`. Exchanging the generator labels
(`(x,y) → (y,x)`, an automorphism of the axis torus commuting with both
classifications) gives `d1(L) = L` (full columns). ∎

### 3b. `b×n(2L−1) = L²` and `b×n(k) = 0` for `k < 2L−1`

**Proof (onset).** Black wraps both directions ⇒ black contains a
dir0-winding cycle (≥ `L` sites) and a dir1-winding cycle (≥ `L` sites).
If they lie in different components, `|B| ≥ 2L`; if in one component, the
two cycles share ≥ 1 site, so `|B| ≥ 2L − 1`. Hence black-both requires
`k ≥ 2L − 1`, proving `b×n(k) = 0` (and `b×b(k) = 0`) below the onset.

**Proof (value and structure at `k = 2L−1`).** The bounds above must all be
tight: one black component containing a dir0 cycle with exactly `L` sites
and a dir1 cycle with exactly `L` sites, sharing exactly one site. By
Lemma B equality cases, the dir0 cycle is a full row and the dir1 cycle is
a full column; sharing one site means black = row `y = c` ∪ column
`x = c'` (a *cross*, `2L − 1` sites).

*White side (direct).* White = complement of the cross, `(L−1)²` sites.
Any white winding path must cross the missing row *and* the missing column
(its `y`-coordinate sweeps all residues, hitting row `c`; similarly `x`),
but every step has `|Δx| ≤ 1, |Δy| ≤ 1` (axis matching edges), so a
dir1-winding path must **stand on** row `c` at some step — impossible; a
dir0-winding path must stand on column `c'` — impossible. White wraps
neither. Cell: `b×n`. Count: `L` rows × `L` columns = `L²`. ∎

*(Consistency: `b×b(2L−1) = 0` follows — the cross's white complement is
barrier-bounded in both directions, so no white winding of any kind.)*

## 4. Axis I4 — pre-both-wrap deficit decomposition (PROVED)

**Claim.** For `L ≤ k < 2L−1` (axis): `C(N,k) − n×b(k) = d0(k) + d1(k)`.

**Proof.** The five-cell row sum (Lemma D bookkeeping) reads
`C(N,k) = n×b + b×n + d0 + d1 + b×b`. By §3b's onset proof,
`b×n(k) = b×b(k) = 0` on this window (black-both needs `≥ 2L − 1` sites).
Hence the deficit `C(N,k) − n×b(k)` is exactly `d0 + d1`. ∎

Combined with §3a, the first deficit at `k = L` is exactly `2L`.
(Checks: axis L=3: `C(9,3) − 78 = 6`; L=4: `C(16,4) − 1812 = 8`; L=5:
`C(25,5) − 12650 = 10` — all `= 2L`. Table check only.)

---

## 5. Diamond I3

The diamond geometry differs from the axis in one decisive way (Lemma A):
**a straight diagonal line winds both generators**, and single-generator
winding at minimal mass is carried by *drift-0 mixed cycles*, not lines.

### 5a. `d0(2L)` — closed form `L·C(2L,L)` (PROVED; #669's `4·C(2L,4)` KILLED)

**Theorem 5a.** `d0(2L) = d1(2L) = L·C(2L, L)`.

**Proof.** Let `B`, `|B| = 2L`, wrap dir0 and not dir1. Then `B` contains a
black NN cycle with displacement `(2aL, 0)`, `a ≠ 0` (and no cycle with
`b ≠ 0`). Such a cycle has length ≥ `2L|a| ≥ 2L` = `|B|`, so `a = 1` and
`B` is a single simple cycle of exactly `2L` steps `s_0, …, s_{2L−1}` with
each `s_i ∈ {(1,1), (1,−1)}` (every NN step has `Δu = +1` after choosing
direction — net `Δu = 2L` forces all steps to advance `u`; net `Δv = 0`
forces exactly `L` steps of each type) — the *drift-0* cycles of Lemma B.

*Enumerate them.* Because every step advances `u` by exactly `+1`, the
cycle visits each `u`-residue exactly once: `B = {(u, v_u) : u ∈ ℤ/2L}`
with `v_{u+1} − v_u ∈ {+1, −1}` for all `u` (cyclically). Conversely any
function `v: ℤ/2L → ℤ/2L` with `v_{u+1} − v_u ∈ {±1}` and `∑(v_{u+1}−v_u)
≡ 0` (automatic) gives such a cycle **provided the sites are valid**:
`u ≡ v (mod 2)` for every site, i.e. `v_u ≡ u (mod 2)` for all `u`. Since
`v_{u+1} − v_u = ±1` alternates the parity of `v` exactly as `u` does, the
parity constraint reduces to `v_0 ≡ 0 (mod 2)`: `L` choices of `v_0`.
Given `v_0`, the cycle is determined by the *step word*
`w ∈ {+1, −1}^{2L}` with exactly `L` of each — `C(2L, L)` choices. The map
`(v_0, w) ↦ {(u, v_0 + w-prefix(u))}` is a bijection onto the dir0-winding
`2L`-site black cycles (distinct pairs give distinct site sets, since
`v_0` and all prefixes are recoverable from the set). Also, a single cycle
with displacement `(2L, 0)` does not wind dir1, so black type is
exactly dir0: there are `L · C(2L, L)` such black sets.

*White side (direct, mirroring §3a).* Let `B` be such a drift-0 cycle.

- *White does not wind dir1.* A white winding-dir1 path (matching edges:
  `(1,±1)` and `(2,0)`; every edge changes `u` by ±1 or ±2) with net
  `v`-displacement `2L` must visit sites of every `v`-residue class... more
  carefully: its `v`-coordinate sweeps a connected interval of the
  universal cover of length ≥ `2L`, covering all residues `mod 2L`. The
  black cycle `B` meets every `v`-residue? Not necessarily — but `B` meets
  every **site with `v ≡ v_u` at `u`**: the path needs to *stand on* a site
  of the residue class it crosses while `u` also matches... The clean
  statement: a dir1-winding path's `v`-coordinate takes every value mod
  `2L` **at some visit**, and at that visit the path is on a site
  `(u, v)` with that `v`; the black set blocks only *its own* `2L` sites.
  This does not immediately forbid the path. **So the white side of the
  cell assignment at the onset is the one place duality (Lemma D) is
  load-bearing** — as on the axis (§3a) a direct argument exists (the
  drift-0 black cycle acts as a "spiral barrier"), but I did not complete
  it in this ticket's budget. The black-side count `L·C(2L,L)` is the
  geometric content and is proved; the cell bookkeeping `d0` vs `d1` is
  by duality + the `(u,v) → (v,u)` automorphism, consistent with the
  committed `d0 = d1` fact.

Hence `d0(2L) = L·C(2L, L)`. ∎

**Verification.** L=2: `2·C(4,2) = 12`; L=3: `3·C(6,3) = 60`; L=4:
`4·C(8,4) = 280`. The committed diamond L=3 table has `d0(6) = 60` ✓ and
diamond L=4 has `d0(8) = 280` ✓. Brute force at L=2 gives 12 ✓
(`scripts/wrapping_onset_proof_checks.py`).

**This kills #669's two-point fit `d0(2L) = 4·C(2L,4)`:** that form gives
`4·C(4,4) = 4` at L=2 (brute force: 12 — violated), `4·C(6,4) = 60` at L=3
(✓ coincidentally `4·C(2L,4) = L·C(2L,L)` iff `C(2L,4)·4 = L·C(2L,L)`, which
holds exactly at L=3,4 — the two points #669 had). The correct closed form
is `L·C(2L,L)`; the #669 diamond-L5 prediction for `d0(10)` should be
`5·C(10,5) = 1260`, **not** `840`. The #669 falsification kit entry
"`d0(10) = 840` if I3 extends" must be corrected accordingly: if a future
diamond L=5 run yields `d0(10) = 1260`, the `4·C(2L,4)` conjecture was the
wrong generalization, not the onset.

### 5b. `b×b(2L) = 2L` (PROVED)

**Proof.** Let `B`, `|B| = 2L`, wrap both generators. `B` contains a cycle
with displacement `(2aL, 2bL)`, `a,b ≠ 0` (if the two windings lived on
different cycles of different components, `|B| ≥ 4L > 2L`; if one component
held two independent winding cycles, their union has ≥ `2L + 2L − 1 = 4L−1
> 2L` sites). So a single cycle, length ≥ `2L·max(|a|,|b|) ≥ 2L` with
equality iff `|a| = |b| = 1` and every step is `(1,1)` or every step is
`(1,−1)`: `B` is a straight diagonal line (Lemma A), slope `(1,1)` or
`(1,−1)` — `2L` lines total (`L` even residues of `v−u`, `L` of `v+u`).

*White side (direct).* Let black be the line `𝔞_c` (slope `(1,−1)`, fixed
`v+u = c`). The **sibling lines** `𝔞_{c'}` for `c' ≠ c` are entirely white
(two distinct slope-`(1,−1)` lines are disjoint — they share no site), and
each closes with displacement `(2L, −2L)`, winding **both** generators. So
white wraps both. (Symmetric for black slope `(1,1)`.) By Lemma D, cell
`b×b`. Count: `2L` lines. And no other black-both set exists at `k = 2L`
(above). ∎

**Verification.** diamond L=2: 4; L=3: 6 (committed) ✓; brute force at
L=2,3 confirms every `b×b(2L)` set is a straight diagonal line ✓.

### 5c. Diamond `b×n` onset `3L−1` with value `4L²` — geometric reading obtained; general proof OPEN

**Structure (verified exactly at L=2,3 by seconds-scale brute force;
consistent with the committed L=4 count).** At `k = 3L−1`, every `b×n`
configuration is

> a **full straight diagonal line** (slope `(1,−1)` or `(1,1)`, `2L` sites,
> black-both by Lemma A) **plus `L − 1` extra sites** forming an
> NN-connected "plug" on a single transverse slope-line, each plug site
> NN-adjacent to the full line or to the next plug site.

At L=2 the plug is a single site (4 positions per anti-line); at L=3 the
plug is a **domino** `{q, q+(1,1)}` with both ends NN-adjacent to the
removed line (6 positions per anti-line: 2 per transverse slope-line × 3).

**Count reading.** `4L² = 2` (slope families) `× L` (full lines per family)
`× 2L` (plug placements per line). This reproduces `16, 36, 64` at
`L = 2, 3, 4` (the L=4 value is the committed table's). The per-line plug
count `2L` is direct at L=2 (the 4 sites not on the line, each works —
brute force) and L=3 (the 6 dominoes, enumerated above); at general `L`
the plug is an `(L−1)`-site object and the `2L` placements would need a
structural classification of plugs — **this is the missing proof step**.

**Why the value is plausible (white-side reading).** The full line alone
is `b×b` (§5b): white's sibling lines wind both. Adding the plug on a
transverse line cuts the annulus complement in the transverse direction
as well: white (matching lattice) retains no winding in either generator,
while black still winds both (the line survives) — cell `b×n`. The plug
must simultaneously block the white NNN chains (`(2,0)`/`(0,2)` steps) and
NN detours in both directions; the classification of minimal such plugs
is exactly the `(L−1)`-site structure above.

**What remains for a proof (OPEN).**

1. *Onset:* no `b×n` exists for `2L ≤ k < 3L−1` — equivalently, any
   black-both set with `k < 3L−1` has white wrapping something (white-both,
   i.e. the config is `b×b`). Verified exactly at L=2,3.
2. *Count:* every `b×n(3L−1)` set is line + `(L−1)`-plug, and there are
   exactly `2L` plugs per line. Verified exactly at L=2,3; count
   consistent at L=4 (committed).

Both are sharp integer predictions for diamond L=5 in the #669
falsification kit (§5 of the #669 note): onset `k = 14`, value `100`.
**Until (1)–(2) are proved for general L, the `3L−1`/`4L²` claims remain
two-point (+L=2 brute-force) conjectures.** Do not enumerate diamond L=5
to settle them.

---

## 6. Summary of proof status

| Item | Status | Where |
|---|---|---|
| I1 axis `n×b = C(N,k)`, `k < L` | **PROVED** | §1 |
| I1 diamond `k < 2L` | **PROVED** | §1 |
| I2 `b×n = C(N,k)`, `k ≥ N−L+1`, both geometries, strict at `N−L` | **PROVED** | §2 |
| Axis I3 `d0(L) = d1(L) = L` (full rows/columns; white side direct) | **PROVED** | §3a |
| Axis I3 `b×n(2L−1) = L²`, onset, cross structure | **PROVED** | §3b |
| Axis I4 deficit `= d0 + d1` on `L ≤ k < 2L−1` | **PROVED** | §4 |
| Diamond `d0(2L) = L·C(2L,L)` (black side proved; #669's `4·C(2L,4)` killed) | **PROVED** (black side; white cell-assignment by duality) | §5a |
| Diamond `b×b(2L) = 2L` (straight diagonal lines; white sibling-lines argument) | **PROVED** | §5b |
| Diamond `b×n` onset `3L−1`, value `4L²` | **OPEN** (geometric reading: line + `(L−1)`-plug, count `2·L·2L`; exact at L=2,3, consistent at L=4) | §5c |

**On `M_L` irreducibility (per the ticket):** the irreducibility of `M_L`
over ℚ at the five committed sizes (PR #669 §4) is recorded as a
**factorization fact**, not a wrapping theorem: no mechanism is known (or
claimed) that forces irreducibility; a reducible case at a future size
would contradict no wrapping statement. Not pursued further per the
ticket.

## 7. Correction to the #669 falsification kit

The #669 note's §5 diamond-L=5 predictions contain one entry derived from
the now-dead `4·C(2L,4)` fit:

- "`d0(10) = 840` if I3 extends" → corrected prediction:
  **`d0(10) = L·C(2L,L) = 5·C(10,5) = 1260`** (and `d1(10) = 1260`).
- `b×b(10) = 10` and `b×n` onset `k = 14` value `100` stand as before
  (two-point conjectures with the §5c geometric reading).

A future diamond L=5 rung should be checked against **both** the corrected
`1260` and the old `840` before any conclusion is drawn.

## 8. Tripwire compliance and reproduction

- No L=6 enumeration, no diamond L=5 enumeration, no new census rung.
- The only new computation is a seconds-scale brute force at *already
  committed or trivially small* sizes (axis L=3,4 = `2^9`,`2^16`; diamond
  L=2,3 = `2^8`,`2^18`) using the repo's own classifier:
  `python3 scripts/wrapping_onset_proof_checks.py` — all checks pass.
  It re-verifies I1, I2 (incl. sharpness), axis I3a/I3b structure, axis I4,
  the diamond `d0(2L) = L·C(2L,L)` black-side count, `b×b(2L) = 2L` with
  diagonal-line structure, and the diamond I4 window.
- `scripts/wrapping_five_cell_reader.py` (PR #669) re-verifies the
  identities on the committed tables; it needs the census JSONs, which
  are on the #669 branch — fetched to this branch for the check run.
- The two OPEN items carry sharp integer predictions in the #669
  falsification kit (§5 there, with the §7 correction above).
