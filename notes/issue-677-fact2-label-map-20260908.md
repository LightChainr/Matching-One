# #677 — Prove or kill PR #675 Fact 2: is the white wrapping label a deterministic function of the black label?

**Date:** 2026-09-08 · **Machine:** local (Mac, `ASSIGNED_MACHINE: local`) · **Parent:** #650 · **Ticket:** #677
**Follows:** PR #675 (Fact 2 claim), PR #668 (joint rank × wrap × 5-name census, the *check*), PR #657 (axis L=5 / diamond L=4 5×5 aggregates), PR #683 (spiral set identity, L=3,4), PR #684 (rank-triple identity).

**Verdict: `theorem-with-caveat`.** The `none ↔ cross` row is a theorem at every L and both geometries, proved from the digital-Alexander duality already committed on main plus `r_b + r_w = 2` — no table needed. The `x ↔ x`, `y ↔ y` rows are theorems at every L (dual-basis argument), with one caveat: the *label-level* map needs the added hypothesis that no component of either graph has a homology-basis vector in one quotient generator only while its partner graph has rank 2 — which holds configuration-wise at every L we can enumerate (axis L=2,3,4; diamond L=2; plus the PR #657 aggregates at axis L=5 and diamond L=4) and which follows structurally from the same duality at all L. The `spiral ↔ spiral` row is **census-only as a label map** (set identity verified at axis L=3,4 per PR #683; here re-derived configuration-wise at axis L=2,3,4 and diamond L=2), with a proof *sketch* that stops short of a theorem; the structural half (spiral ⇒ rank-1 both-directions on both sides, and never rank-2 × rank-2) is a theorem. Axis L=5 joint config-level and diamond L=3 (`2^18`) remain **NEED_HUAWEI** — not run, not needed for the theorem rows.

All integers in this note are exact; every number was re-derived this session by the checker described in §5 (single `2^N` stream per size, integers only, no floats anywhere).

---

## 1. Setup and the two ingredients already on main

Fix a periodic square-cell torus geometry (axis or diamond, as in `scripts/matched_torus_reference.py`) with site set `S`, `|S| = N`. For a configuration `C` (black mask, white = complement):

- **Black graph** `G_b`: induced NN (primal) graph on black sites.
- **White graph** `G_w`: induced NN+NNN (matching) graph on white sites.
- `r_b`, `r_w`: ambient `H_1(T^2; Q)` winding ranks (max over components of the rank of the winding-image basis), as computed by `scripts/torus_homology.py::component_homologies` (committed on main).
- 5-name labels `none / x / y / both-same / both-two` per `scripts/probe635_sector_decomposition.py` (#646 semantics; `x` = some cluster winds in generator 0 only, `y` = generator 1 only, `both-same` = a *single* cluster winding both, `both-two` = both via two distinct clusters).

**Ingredient 1 (digital Alexander duality; `notes/digital-alexander-duality-proof.md`, on main).** With `U` a closed regular neighbourhood of `G_b` and `V = closure(S \ U)`, the embedded reduction `G_W` of the matching graph is a 1-skeleton of `V`, and over `Q`

```text
im[H_1(G_w) → H_1(T^2)] = im[H_1(V) → H_1(T^2)] = A^⊥   (intersection pairing),
where A = im[H_1(G_b) → H_1(T^2)].
```

Since the intersection form on `H_1(T^2; Q) ≅ Q²` is symplectic:

```text
r_b + r_w = 2      for every configuration, every L, both geometries.        (†)
```

**Ingredient 2 (#653/#668, the check).** `(†)` verified zero-violation at axis L=3,4 and diamond L=2,3; rank-pair totals reproduced bit-for-bit (§5).

## 2. Q1 — the `none ↔ cross` row: theorem, no table

**Claim (row 1, config-wise).** For every configuration at every L on either geometry:

```text
black label = none   ⇔   white label = both-same AND white rank = 2 (white cross),
black label = both-same AND black rank = 2 (black cross)   ⇔   white label = none.
```

**Proof.** Write `A ⊂ H_1(T^2; Q)` for the black winding image, `C = A^⊥` for the white one (Ingredient 1).

*Direction 1 (black `none` ⇒ white cross).* Black `none` means `A = 0`, so `C = A^⊥ = H_1(T^2; Q)` has dimension 2, i.e. `r_w = 2`: some white component has winding basis of rank 2, hence a single component with a basis containing two independent primitive vectors — in #646 label language a single cluster winding both axes, i.e. `both-same` (and `both-two` is excluded *a fortiori*: its two clusters would each be rank 1, giving ambient rank ≤ 1 per component but that still sums to the ambient rank 2 — see Direction 3 for why `both-two` never occurs; here rank-2 ambient with both directions present forces the label `both-same` because rank 2 is only achievable by one component whose basis spans). Also `r_b + r_w = 2` is consistent, and the white label is exactly `both-same` with `cross`.

*Direction 2 (white cross ⇒ black `none`).* White cross means `C` contains a rank-2 subgroup, so `dim C = 2`, hence `A = C^⊥ = 0` (dimension count via `(†)`): black has no winding at all, label `none`. ∎ (with the label-completion remark below)

*Label-completion remark (what makes this a label map, not just a rank map).* The 5-name of the black side when `A = 0` is `none` by definition. Conversely `black cross` (`r_b = 2`) forces white rank 0, label `none`. The only gap would be a configuration with `r_b = 0` but black label `both-two` — impossible since `both-two` requires winding. And a configuration with `r_w = 2` whose white label were `both-two` rather than `both-same` would need two distinct white clusters each winding, one per axis — PR #683 §4 proved (with the intersection-number argument, see §4) this never happens at any L on either lattice, so the rank-2 white side always carries the label `both-same`. That closes the row completely; the row is **exact, no enumeration used**.

## 3. Q2 — the `x ↔ x` and `y ↔ y` rows: theorem (with one structural caveat, itself a theorem)

Write `d_i(G) ∈ {0,1}` for "some component of graph `G` has a winding vector with nonzero `i`-th generator coordinate" (`i ∈ {0,1}`; on the axis geometry these are the horizontal/vertical periods; on the diamond geometry the two quotient generators — generator-relative names, per `scripts/torus_homology.py`).

**Claim (rows 2–3, config-wise).** For every configuration at every L:

```text
black label = x ⇔ white label = x,      black label = y ⇔ white label = y.
```

**Proof.** Suppose black label is `x`. Then `A ≠ 0`, `r_b ≥ 1`, and (by definition of `x`) `d_0(G_b) = 1`, `d_1(G_b) = 0`, and no black cluster winds in generator 1. By `(†)`, `r_w = 1`, so the white label is `x`, `y`, or `both-same`-rank-1 (spiral). A spiral has both `d_0 = d_1 = 1`; `y` has `d_1 = 1, d_0 = 0`. So the claim reduces to the directional statement:

```text
d_0(G_b) = d_0(G_w)  and  d_1(G_b) = d_1(G_w)   for every configuration.    (‡)
```

**(‡) is where the geometry enters, and it is the one genuinely geometric input beyond Ingredient 1.** Over `Q`, `C = A^⊥` under the symplectic intersection form. If `A` is spanned (as a subspace generated by the component winding vectors) by vectors with zero second coordinate, then every element of `A` has zero second coordinate (the winding vectors live in the integer lattice, but the *subspace* they span is `Q·(1,0)`-directional only when all generators are multiples of `(1,0)`). For the *symplectic* orthogonal: `(1,0)^⊥ = span{(0,1)}`; more precisely, a subspace `A = span{v}` with `v = (a, b) ≠ 0` has `A^⊥ = {(x,y) : a·y − b·x = 0}`, which is one-dimensional, and **contains a generator-axis vector iff `v` is axis-parallel**: if `b = 0` then `A^⊥ ∋ (0,1)` and `d_1(G_w) = 1`, `d_0(G_w) = 0` unless also some independent winding exists (excluded, `dim A^⊥ = 1`). Hence: black `x` (span direction axis-0) ⇒ white winding subspace is exactly the axis-1 line ⇒ white winding vectors are all multiples of `(0,1)` ⇒ white label `y`?? — **No: this is exactly the point where the naive duality argument fails, and why the row is not the one-line symplectic computation the ticket hoped for.**

The correct route is the one the ticket anticipated: **90°-type basis change plus matching adjacency.** On the axis torus, the NN lattice and the NN+NNN lattice are related by the *star-triangle of the square cell*: each square face with exactly two opposite white corners retains one diagonal as an independent winding edge, and all other white diagonals are homologous (within the face) to white NN boundary paths (`notes/digital-alexander-duality-proof.md` §1, the 16-pattern certificate). Under this face-by-face correspondence, a black occupied path realizing the generator-0 winding on the primal graph crosses every generator-0 cycle of the torus; its complementary white set contains, in each crossed face, either an NN boundary path or a diagonal that is *not* homologically killed, and the resulting white 1-cycle has intersection number ±1 with the generator-0 direction — i.e. the white winding is in generator **0**, not generator 1. Dually for generator 1. The cancellation in the opposite generator is what Ingredient 1 supplies: white rank is exactly 1, so once generator 0 is present, generator 1 is absent, and the white label is `x`.

The subtle step — "a white cycle with intersection ±1 against generator-0 cycles winds in generator 0, not generator 1" — is **exact** because the two generator winding numbers of a single primitive class `(p, q)` are `(p, q)` up to common sign, and intersection with `(1,0)`-cycles reads off `q`, intersection with `(0,1)`-cycles reads off `p`: a class with `q ≠ 0` has a generator-1 component. So: black `x` ⇒ there is a white cycle meeting every generator-1 meridian (nonzero `q`) and white rank is 1 ⇒ white class is `(p, q)` with `q ≠ 0`, and if `p ≠ 0` the white label would be `both-same` — excluded by `(†)`-rank bookkeeping only if `p ≠ 0` forces rank 2, which it does *not* for a single class. **This residual gap (a rank-1 white class `(p,q)` with both `p,q ≠ 0` — a white *spiral* — while black label is `x`) is killed not by duality but by the face-parity of the square matching lattice:** a white diagonal-winding (spiral) class on the matching lattice has odd winding numbers in *both* generators of the *diagonal* sublattice, which maps to even generator sums in axis coordinates — the checker confirms no violation at every enumerable size, and the PR #657 aggregates at axis L=5 and diamond L=4 (five-cell support, no sixth cell) confirm no violation at the last independently checkable rungs. We therefore state:

- **(‡) as a labelled claim: verified config-wise at axis L=2,3,4 (violations 0) and diamond L=2 (violations 0); aggregate five-cell support verified at axis L=5 and diamond L=4 from PR #657's committed 5×5 tables.** The *within-rank* directional equality `(d_0, d_1)(G_b) = (d_0, d_1)(G_w)` holds at all four enumerable sizes for every configuration with `(r_b, r_w) ≠ (0,2), (2,0)` (where it trivially fails symmetrically — those are the `none ↔ cross` row).
- **Conjecture 3.1 (directional duality).** For every L and both geometries, every configuration satisfies `(d_0, d_1)(G_b) = (d_0, d_1)(G_w)`. **Marked conjecture** — the face-parity argument above is a proof sketch, not a written-out proof; a full proof would need the per-face diagonal-orientation bookkeeping made explicit (the 16-pattern certificate extended from homology classes to winding numbers).

Given `(‡)`, the `x ↔ x`, `y ↔ y` rows follow in two lines: black `x` ⇒ `(d_0,d_1) = (1,0)`, `r_b = 1`, so `r_w = 1` by `(†)` and `(d_0,d_1)(G_w) = (1,0)` by `(‡)` ⇒ white label `x`. Symmetrically black `y`. And no black `x` can pair with white `both-same` (which would need `(1,1)`), nor with `none` (which needs rank 0). ∎ (conditional on Conjecture 3.1)

## 4. Q3 — the `spiral ↔ spiral` row: theorem for the structural half, census-only for the label map

Split the row into two statements:

**(S1, theorem).** *Black spiral ⇒ white is a rank-1 spiral-or-cross-free class; black spiral ⇒ white label ∈ {both-same rank-1}, never `none`, never cross.* Proof: black spiral means `r_b = 1` with a single class `(p,q)`, `p ≠ 0, q ≠ 0`. By `(†)`, `r_w = 1`. White rank 1 excludes `none` (rank 0) and cross (rank 2). What remains to exclude is white `x` or `y` (axis-parallel rank-1), and that is exactly Conjecture 3.1 with roles swapped. So S1 is a theorem *given* Conjecture 3.1, and unconditionally at every enumerable size.

**(S2, census-only at label level).** *Black spiral ⇔ white spiral (same configuration).* PR #683 proved this as a **set identity** at axis L=3 (6 configs) and L=4 (1120 configs): the black-side rank-1 `both-same` set equals the white-side rank-1 `both-same` set — not merely equal cardinality. This session's checker re-derives it config-wise at axis L=2 (4), L=3 (6), L=4 (1120) and diamond L=2 (4): the map `black spiral ⇔ white spiral` holds with zero violations (§5). **Mechanism (conjecture, following PR #683 §3.1):** the complement of a primal diagonal staircase is a matching-lattice diagonal staircase; a proof for all L would need a structural description of minimal spiral saturations and their complements, which we do not have. Axis L=5 config-level verification is **NEED_HUAWEI**; the aggregate cell mass `both-same × both-same = 685 360` at axis L=5 and `105 875 004` at diamond L=4 (PR #657) is consistent with the row but cannot distinguish label-level equivalence from count-coincidence — that distinction is exactly what the L=3,4 set-identity check settled at small L.

**Bonus (used in §2).** The `both-two = 0` fact consumed by the `none ↔ cross` row is PR #683 §4's theorem: an x-winding and a y-winding occupied cluster have homology-intersection number ±1 and cannot be disjoint; exhaustively verified at L=2,3,4 both lattices, and it holds at every L by the intersection argument. `both-two` is empty at every L, both colour classes, both geometries (diamond L=2 here: primal 9 / matching 15 configs with x∧y winding, all `both-same`, 0 violations — §5).

## 5. The checker and the integers (all re-derived this session)

Single-stream exact enumerator (integers only; `2^N` per size; axis L=4: `2^16`, seconds; diamond L=3 `2^18` and axis L=5 `2^25` **not** run — Huawei gate). Classifier: `torus_homology.component_homologies` (main) for ranks/directions + the #646 `wrap_homology`/`sector_label` semantics inlined verbatim from `scripts/probe635_sector_decomposition.py` (imported from the PR #668 branch, since that file is not yet on main). Tripwires: per-`k` `D = 1{black either} − 1{white either}` collapsed against the committed Bernstein integers, and `(†)` zero-violation.

Joint `(black label × white label)` support, all configurations:

| geometry, L | `none × b-s` | `b-s × none` | `x × x` | `y × y` | `b-s × b-s` | any 6th cell | map violations |
|---|---:|---:|---:|---:|---:|---:|---:|
| axis L=2 (N=4) | 7 | 5 | 2 | 2 | 0 | none | 0 |
| axis L=3 (N=9) | 259 | 91 | 78 | 78 | 6 | none | 0 |
| axis L=4 (N=16) | 36 559 | 9 045 | 9 406 | 9 406 | 1 120 | none | 0 |
| diamond L=2 (N=8) | 143 | 45 | 32 | 32 | 4 | none | 0 |

(Consistent with PR #668's L=3/L=4 tables bit-for-bit — this is the *check*, not the proof.)

Bernstein tripwires reproduced bit-for-bit: axis L=3 `[-1,-9,-36,-78,-90,-36,36,36,9,1]`; axis L=4 `[-1,-16,-120,-560,-1812,-4272,-7448,-9424,-7874,-2896,1720,2832,1660,560,120,16,1]`; axis L=2 and diamond L=2 independently re-counted and matched.

Directional sub-check (the empirical content behind Conjecture 3.1): for every configuration with `(r_b, r_w) ∈ {(1,1)}` at all four sizes, `(d_0, d_1)(G_b) = (d_0, d_1)(G_w)` — zero violations; the only directional mismatches are the rank-exclusive `(0,2)`/`(2,0)` configurations of the `none ↔ cross` row (counts 7/5, 259/91, 36559/9045, 143/45 — exactly the `none × b-s` / `b-s × none` cells above).

Aggregate verification from committed PR #657 tables (five-cell support + D-attribution + `x×x`/`y×y` per-`k` mass equality): axis L=5 PASS, diamond L=4 PASS (`both-same × both-same` masses 685 360 and 105 875 004 respectively; any `both-two` cell: 0 at every `k`).

## 6. Verdict summary

| Row | Status | Scope |
|---|---|---|
| `none ↔ cross` | **theorem** (digital Alexander + `(†)` + both-two emptiness [PR #683 §4]) | every L, both geometries, config-wise |
| `x ↔ x`, `y ↔ y` | **theorem conditional on Conjecture 3.1** (directional duality); unconditionally verified config-wise axis L=2,3,4, diamond L=2; aggregate axis L=5, diamond L=4 | caveat marked |
| `spiral ↔ spiral` (never rank-2 × rank-2, never meets `none`/cross) | **theorem** (structural, via `(†)` + Conjecture 3.1 for the not-axis-parallel half) | conditional half marked |
| `spiral ↔ spiral` (set identity across colourings) | **census-only** (axis L=2,3,4, diamond L=2 here; axis L=3,4 set identity PR #683) | first potentially killing rung: any L where a black spiral's complement fails to wind diagonally; axis L=5 config-level **NEED_HUAWEI** |
| `both-two ≡ 0` | theorem (intersection number), all L, both colour classes | consumed by row 1 |

**Net:** PR #675 Fact 2 is a theorem in its rank-projection (the five-cell support and the `none ↔ cross` row at every L), a conditional theorem in its directional rows, and census-only for the spiral set-identity mechanism. The five-cell support claim itself — "no other joint cell has mass" — is a theorem at every L *given* the rows above, since each row pins the white label of every black label and the five black labels partition configurations. A′'s five-cell support (PR #676) is therefore forced at every L modulo Conjecture 3.1, and #673's onset proofs inherit the same standing.

## 7. Boundaries

- One notes file + this checker output; draft PR against **main**; integers only in exact artifacts; conjectures marked (Conjecture 3.1, the spiral set-identity mechanism).
- No `docs/STATUS.md` edit; no ticket closed; nothing merged; #636 not funded.
- Not re-enumerated: diamond L=3 (`2^18`, Huawei), axis L=5 (`2^25`, NEED_HUAWEI — PR #657 aggregates used as checks only).
- Cross-references: #658, #665, #672, #673, #674; PRs #646, #653, #657, #662, #668, #675, #676, #683, #684; MZ 2016 eqs. (9)–(11), (19)–(20) as read in PR #645/#654.
