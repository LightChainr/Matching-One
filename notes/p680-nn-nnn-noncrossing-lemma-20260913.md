# #680 — Lemma A written out: NN+NNN strip frontier partitions are noncrossing (proof), and a convention correction for the #675 small-width counts

**Ticket:** LightChainr/Matching-One#680. **Date:** 2026-09-13. **Mode:** reasoning +
exact small-width enumeration on the assigned local machine (no Huawei budget used,
no width above 6 enumerated). No `docs/STATUS.md` entry; no ticket closed or merged.

**Answer up front.**

> **Verdict: `proof` for the noncrossing claim, at every width and depth (not just
> w ≤ 6).** The Jordan-curve argument of PR #675 §4 can be completed in full rigour
> once the king graph is *planarized*: the planarized graph has exactly the same
> connectivity as the king graph, and then the classical boundary-terminal
> noncrossing theorem applies. The subtle point — that a path may re-touch the
> frontier between its endpoints — is handled by a joint-minimality choice of the
> crossing tuple, with an explicit case analysis; both cases are written out below
> so a reader who has never enumerated w ≤ 6 can check them.
>
> **The equality of reachable classes (NN vs NN+NNN) remains verified-small-width**
> (here: w = 2..6, R = 1..4, exact integers, brute-force cross-checked), exactly the
> epistemic status it had in PR #675; a general proof is not given and the status is
> marked as such.
>
> **Convention correction (new finding).** The #675 counts "28 states at w=5 over
> 4 rows and 66 at w=6 over 3 rows, identical sets; both strictly inside the
> 42/132 noncrossing partitions" mix two conventions. 28 and 66 are reproduced
> exactly as **NN-side counts of the delayed-closing / fully-occupied-frontier
> convention** (28 at w=5,R=4; 66 at w=6,R=4; the "3 rows" reading gives 65, an
> off-by-one). Under that convention the NN and NN+NNN reachable sets are **not**
> equal (13 vs 28 at w=5,R=4; 24 vs 66 at w=6,R=4) — though all classes are still
> noncrossing. The equality that Lemma A claims IS true — but under the honest
> convention (partition of the occupied frontier sites with *all* bonds, including
> the frontier row's own horizontal bonds), where the reachable class sets coincide
> at every tested (w,R) and the ambient count is not C_w (51 reachable classes at
> w=5 saturation, because vacant frontier positions admit coexistence structures
> that push the count above C_5 = 42).

---

## 1. Setup and statement

Fix w ≥ 1 (strip width) and R ≥ 1 (processed depth). Lattice Λ = {0,…,R−1} ×
{0,…,w−1}: column x = depth index, row y = transverse index, **open** boundary in y.
The frontier is F = {R−1} × {0,…,w−1}; its w sites are the frontier sites. Adjacency:

- **NN**: |Δx| + |Δy| = 1 (von Neumann);
- **NN+NNN (king)**: max(|Δx|, |Δy|) = 1 (adds the two diagonals of each unit face).

For an occupied set S ⊆ Λ, the **frontier occupancy partition** π(S) is the
partition of the occupied frontier sites F ∩ S into connectivity classes of the
graph (S, E). (Vacant frontier sites carry no label.) A partition of a subset of F
is **crossing** if there are positions a < b < c < d and two *different* blocks A, B
with a, c ∈ A and b, d ∈ B; otherwise **noncrossing**.

**Theorem 1 (noncrossing; proved here in full rigour, all w, all R).** For every
S ⊆ Λ and both adjacency rules, π(S) is noncrossing.

**Claim 2 (reachable classes; verified w ≤ 6, R ≤ 4 — not proved here in general).**
Let 𝒞_G(w,R) = {π_G(S) : S ⊆ Λ}. Then 𝒞_king(w,R) = 𝒞_NN(w,R) for all tested
w ≤ 6, R ≤ 4 (§3). Beyond that range the status is unchanged from PR #675:
verified, not proved. (See §6 for what blocks a general argument.)

**Remark 3 (the honest ambient is not C_w).** With vacancies allowed, the honest
ambient (all noncrossing partitions of occupied subsets) *exceeds* C_w: at w = 5
the reachable count saturates at 51 > 42 = C_5. Two structural facts are visible in
the enumeration and trivial to prove: (a) two *adjacent* occupied frontier positions
always lie in a common block (the bond between them is present), so a block's
occupied "interior" frontier positions cannot belong to another block; (b) a block
may leap vacant gaps (e.g. one component touching frontier positions {0,2} with
position 1 vacant — §5). The Catalan comparison in PR #675 therefore only makes
sense in the full-frontier convention of §4.

## 2. Proof of Theorem 1

Throughout, "path" means a simple path in (S, E) for the chosen adjacency, drawn on
the given embedding: NN edges as unit grid segments; king diagonals as curves inside
their unit face (§2.1). Write ∂ = {R−1} × ℝ for the frontier line, H = {x ≤ R−1}
for the closed half-plane containing all of S, and identify the frontier sites with
the *vertices of the graph on ∂* (this is exactly F ∩ S).

### 2.1 Planarization lemma

**Lemma 4.** Let S ⊆ H be finite. There is a graph G* = (S, E*) with E* ⊆ E_king,
the same connectivity relation as the king graph on S, and a plane drawing in which

1. every NN edge is its straight unit segment;
2. every diagonal edge of a unit face is drawn strictly inside that face;
3. at most one diagonal per face is drawn;
4. distinct edges meet only at shared endpoints;
5. all of the drawing lies in H, and the vertices on ∂ are exactly F ∩ S.

*Proof.* Include all NN edges with both endpoints in S. For each unit face f with
corners p₁ = (x,y), p₂ = (x+1,y), p₃ = (x,y+1), p₄ = (x+1,y+1): the two diagonals
are d₁ = p₁p₄ and d₂ = p₂p₃. Include d₁ iff both endpoints are in S and d₂ is *not*
fully occupied-endpointed; include d₂ iff both endpoints are in S and d₁ is not.
*Connectivity is preserved:* if both d₁ and d₂ have occupied endpoints then all four
corners are occupied, so the four NN cycle edges of f are present and either
diagonal is redundant with a 2-edge NN path around the face (e.g. p₂ → p₁ → p₃
replaces d₂); if exactly one diagonal has occupied endpoints it is included verbatim.
*Planarity of the drawing:* grid segments intersect only at shared lattice points;
a diagonal lies in the open interior of its face and meets face-boundary edges only
at corners; two diagonals lie in disjoint face interiors (one per face). All faces
of H lie in H, so the drawing lies in H; a vertex lies on ∂ iff it is a frontier
site. ∎

The same lemma with E_NN in place of E_king (no diagonals at all) is classical. From
now on fix G ∈ {NN, king} and work with its planarization G*.

### 2.2 Choosing a minimal crossing: the paths avoid the frontier

Suppose Theorem 1 fails: there are positions a < b < c < d on ∂ and distinct blocks
A, B of π(S) with a, c ∈ A, b, d ∈ B. Among **all** such crossing data
(a, b, c, d, A, B) together with a simple path P in G* from a to c and a simple path
Q in G* from b to d, choose one minimizing

  N(P, Q) = |P ∩ ∂| + |Q ∩ ∂|   (number of graph vertices of the paths on ∂).

**Lemma 5 (boundary avoidance).** P ∩ ∂ = {a, c} and Q ∩ ∂ = {b, d} (as vertex
sets).

*Proof.* Every vertex of P on ∂ is a frontier site connected to a, hence lies in A;
similarly vertices of Q on ∂ lie in B; and A ∩ B = ∅ because blocks are disjoint as
sets of sites. Let v be a frontier vertex of P other than a, c. Because P is simple
and a, c are its endpoints, v is internal, and P splits at v into simple subpaths
P[a→v] and P[v→c]; replacing the endpoint a (resp. c) of the pair by v drops at
least one frontier vertex, so any resulting valid crossing datum has strictly
smaller N. We exhibit a valid one in each case:

- **v < a.** Then v < a < b < c < d, so (v, c) ⊂ A and (b, d) ⊂ B alternate
  (v < b < c < d); use P[v→c]. N drops (a removed).
- **a < v < b.** (v, c) ⊂ A with (b, d) ⊂ B alternate (v < b < c < d); use P[v→c].
- **b < v < c.** (a, v) ⊂ A with (b, d) ⊂ B alternate (a < b < v < d, since
  v < c < d); use P[a→v].
- **c < v.** (c, v) ⊂ A with (b, d) ⊂ B alternate (b < c < d < v); use P[c→v].

Symmetrically for Q and its frontier vertices v′ ∉ {b, d} (v′ ∈ B; in each case a
subpath of Q has endpoints that alternate against (a, c) ⊂ A):

- **v′ < a.** (v′, b) ⊂ B with (a, c) ⊂ A alternate: sorted, v′ < a < b < c reads
  B A B A. Use Q[v′→b] (drops d).
- **a < v′ < b.** (v′, d) ⊂ B with (a, c) ⊂ A alternate (a < v′ < c < d); use
  Q[v′→d].
- **b < v′ < c.** (v′, d) ⊂ B with (a, c) ⊂ A alternate (a < v′ < c < d); use
  Q[v′→d].
- **c < v′ < d.** (b, v′) ⊂ B with (a, c) ⊂ A alternate (a < b < c < v′); use
  Q[b→v′].
- **d < v′.** (b, v′) ⊂ B with (a, c) ⊂ A alternate (a < b < c < v′); use Q[b→v′].

Each case yields a crossing datum with a strictly smaller N, contradicting the
minimality. Hence no such v (resp. v′) exists. ∎

*(The five Q-cases and four P-cases above are exhaustive: v is either below a,
strictly between a and c — split by comparison with b — or above c; likewise for
v′ against b, d with the intermediate comparison against c. In every case the
alternation check is the displayed strict-order chain, and the subpath used is the
half of P or Q that ends at the *other* original endpoint's side, which is what
makes N strictly smaller.)*

### 2.3 The Jordan contradiction

Let P, Q be the minimal paths of §2.2.

**(a) The curve.** C := P ∪ [a, c], where [a, c] is the straight segment of ∂ from
a to c. P meets ∂ only at a and c (Lemma 5), and P ∩ [a,c] = {a, c}, so C is a
Jordan curve in the sphere (equivalently: compactify H with the point ∞).

**(b) Q avoids C except at b.** Q ∩ P = ∅: a shared vertex w would lie in A and B;
distinct edges of G* cannot cross or overlap (Lemma 4(4)). Q ∩ [a,c] = {b}: any
edge of Q with an endpoint on ∂ has that endpoint in {b, d} (Lemma 5 again — an
edge lying *along* ∂ has both endpoints on ∂, and the only candidate d = b + 1
would leave no integer for c, since b < c < d forces d ≥ b + 2). So
Q ∩ C = {b}, and Q⁺ := Q − {b} is a connected set disjoint from C containing d.

**(c) The left side of the open segment is the bounded side.** Since P ∩ [a,c] =
{a, c} and P, [a,c] are compact with P ∩ [a, c] = {a, c}, there are η, δ > 0 with
dist(P, [a+η, c−η]) ≥ δ. Take p = (−ε, y₀) (coordinates relative to ∂) with
y₀ ∈ (a+η, c−η) and ε < δ. The horizontal ray from p towards ∂ crosses C exactly
once — at (0, y₀) on the segment; it cannot meet P inside the δ-tube, and C has no
points on the far side of ∂. Odd crossing parity puts p in the component I of the
complement that does **not** contain ∞.

**(d) The two ends of Q lie in different components.** The first edge of Q leaves b
into the open half-plane (it cannot run along ∂, by Lemma 5), so Q⁺ contains points
arbitrarily close to b on the left of the open segment (a, c) — at heights within
(a+η, c−η) for small excursions, since b ∈ (a, c). By (c) those points lie in I, so
Q⁺ ∩ I ≠ ∅. On the other hand, the arc ∂* := ∂ ∪ {∞} minus the open segment (a, c)
is connected, meets C only in {a, c}, and contains d and ∞; hence
∂* − {a, c} lies in the single component of the complement that contains ∞, so
**d lies in the ∞-component**, and Q⁺ (connected, avoiding C, containing d) meets it
too. A connected set disjoint from C cannot meet two different components. Contradiction.

**(e)** The contradiction rules out the assumed crossing; π(S) is noncrossing. □

**Remark 6 (what the proof uses).** Only: (i) the occupied graph planarizes with
frontier vertices exactly on ∂ (Lemma 4 — this is where "both diagonals of a face
active ⇒ all four sites occupied ⇒ the NN edges already merge the would-be
crossing" is cashed out); (ii) the frontier line is a straight boundary and the
strip is open in y. **The cyclic (periodic-in-y) strip is *not* covered**: there the
frontier is a circle, linear noncrossing does not imply circular noncrossing, and
Lemma 4's half-plane drawing is unavailable. The periodic convention was not part of
Lemma A's statement and is left untouched here.

**Remark 7 (geometric crossings vs partition crossings).** In the straight-line
drawing the king graph's edges do cross (at face centres, the two diagonals of one
face). Theorem 1 is about the *partition*, not the drawing: the planarization
replaces the crossing pair of diagonals by the NN 2-hop around the face whenever
both are active, which is possible precisely because the four corners are then all
occupied. So "the king graph is non-planar" and "the king frontier partitions are
noncrossing" are consistent; the latter is what any transfer-matrix state space
cares about.

## 3. Enumeration control (honest convention, w = 2..6, R = 1..4)

Method (exact integers only):

- **Transfer BFS.** States are (occupancy mask of the current frontier column,
  canonical restricted-growth string of the partition of its occupied sites).
  Transition = append a column with any of the 2^w masks, union-find update with the
  chosen adjacency (within-column: consecutive bonds, both rules — king adjacency is
  ℓ∞-distance 1, so *no* within-column distance-2 bonds; the distance-2 frontier
  bonds that diagonals create go *through* the previous column, exactly as in the
  row-by-row build of PR #645 Q2). After R columns, the reachable states are exactly
  {π(S) : S ⊆ Λ}: every configuration induces one state, and every state path is
  realized by its prefix configurations.
- **Brute-force cross-check.** Independent per-configuration union-find over all
  2^(wR) configurations for every w·R ≤ 16; agrees with the BFS on every such cell
  (both rules). The R ≥ 3 columns at w = 5, 6 rely on the BFS alone (2^18–2^24
  configurations are beyond the Mac budget for the brute force), but the BFS
  transition itself is the exact connectivity update.
- **Noncrossing check.** Every reachable class in every cell is noncrossing (linear
  test), for both rules — machine confirmation of Theorem 1 on the enumerated box.

Main table (class counts; "equal" = the NN and NN+NNN class *sets* coincide):

| w \ R | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| 2 | 4 = 4 ✓ | 4 = 4 ✓ | 4 = 4 ✓ | 4 = 4 ✓ |
| 3 | 8 = 8 ✓ | 9 = 9 ✓ | 9 = 9 ✓ | 9 = 9 ✓ |
| 4 | 16 = 16 ✓ | 21 = 21 ✓ | 21 = 21 ✓ | 21 = 21 ✓ |
| 5 | 32 = 32 ✓ | 50 = 50 ✓ | 51 = 51 ✓ | 51 = 51 ✓ |
| 6 | 64 = 64 ✓ | 120 = 120 ✓ | 127 = 127 ✓ | 127 = 127 ✓ |

Class-count sequences (w = 2..6):

- R = 1: NN [4, 8, 16, 32, 64]; NN+NNN [4, 8, 16, 32, 64]  (both = 2^w: single
  column, partition = runs of the occupied subset)
- R = 2: NN [4, 9, 21, 50, 120]; NN+NNN [4, 9, 21, 50, 120]
- R = 3: NN [4, 9, 21, 51, 127]; NN+NNN [4, 9, 21, 51, 127]
- R = 4: NN [4, 9, 21, 51, 127]; NN+NNN [4, 9, 21, 51, 127]  (saturated from R = 3
  at every tested w; the counts do not move between R = 3 and R = 4)

Readings: saturation by R = 3 at w ≤ 6; the saturated sequence 4, 9, 21, 51, 127 is
the count of noncrossing partitions of occupied subsets subject to Remark 3(a)
(adjacent occupied positions share a block) — it exceeds C_w precisely because of
the vacancy-leap structures of Remark 3(b). Raw class lists per cell:
`results/p680-noncrossing-lemma-20260913/raw/w{w}_R{R}/{nn,nnn}.txt`.

## 4. Reconciling the #675 numbers (28 / 66) — a convention correction

PR #675 §4 quotes "28 states at w=5 over 4 rows and 66 at w=6 over 3 rows, identical
sets; both strictly inside the 42/132 noncrossing partitions". Enumerating the
natural variants isolates the convention that produces these numbers:

**Delayed-closing / full-frontier convention.** Row-by-row build; when a row is
added, the *previous* row is closed (its own horizontal bonds applied) first, then
the vertical bonds and (for king) the diagonal bonds between prev and new are
applied; the *new* row's own horizontal bonds stay pending. Only configurations
whose final row is fully occupied are accepted, so classes are partitions of the
full position set [w] — this is the only reading in which the ambient C_w = 42/132
comparison is meaningful. Direct per-configuration enumeration validates the BFS
for this convention too (checked at (w,R) ∈ {(3,3),(4,3),(4,4),(5,3),(5,4),(6,3),(6,4)}, both rules).

| w, R | NN | NN+NNN | sets equal? |
|---|---|---|---|
| 5, 4 | **28** | 13 | no |
| 5, 3 | 28 | 13 | no |
| 6, 4 | **66** | 24 | no |
| 6, 3 | 65 | 24 | no |

- 28 at (w=5, R=4) reproduces **exactly**; 66 appears at (w=6, R=4), while the
  "3 rows" reading gives 65 — one state off, most plausibly an off-by-one in the
  #675 note's row count (or a minor bond-timing detail).
- Under this convention the NN and NN+NNN reachable sets **differ in both
  directions**: e.g. at (3,3) the NN class 0|1|0 (sites {0,2} joined, {1} separate)
  is not king-reachable-with-full-frontier, while at (4,4) the king class
  0011 (blocks {0,1},{2,3}) is not NN-reachable — king diagonals both merge states
  NN keeps apart and (via full-frontier relays being unavailable to NN at the same
  depth) fail to reach some NN states. All classes of both rules remain noncrossing
  in every cell.
- Therefore #675's "identical sets" and "28/66" cannot both be literally right. The
  equality claim of Lemma A **is** true — in the honest convention of §1/§3, which
  is what the lemma's own wording ("the set partition induced on the frontier sites
  by occupied-site connectivity") describes, since there the frontier row's own
  bonds are part of the connectivity relation. The 28/66 figures are the NN side of
  the delayed-closing variant.

No claim above requires width > 6; nothing was enumerated beyond w = 6.

## 5. The converse stays false at the cluster level (machine-checked witness)

The recorded #675 fact survives: diagonals create cluster *footprints* that NN
cannot make with the same occupied set. Witness on the 2×3 rectangle (depth 2,
width 3), occupied set {(1,0), (0,1), (1,2)}:

- king: one cluster {(1,0),(0,1),(1,2)}, frontier footprint **{0,2}** (position 1
  vacant — the diagonal relay through (0,1));
- NN with the same three sites: three singletons, footprints {0}, {2}, {1}→∅-on-frontier.

So NN cannot produce that footprint with those three sites. The *partition*
{0,2} (with 1 vacant) is nevertheless NN-reachable at depth ≥ 2 by occupying a
longer relay (e.g. (0,0),(0,1),(0,2)), which is exactly why the *class sets* of §3
coincide while the cluster-level statement stays false. Both halves are checked by
`scripts/p680_frontier_partition_enumeration.py` (witness function).

Per the ticket: no claim is made about *weighted* counts or realizations-per-state;
only the partition classes are compared.

## 6. What this does not do

- Does not prove Claim 2 (class-set equality) in general; it is machine-verified at
  w ≤ 6, R ≤ 4 only, in the honest convention. The obstruction to a cheap general
  proof: king→NN rerouting of a same-depth configuration needs relay sites that can
  collide across blocks, and the delayed-closing counterdirection at §4 shows some
  rerouting statements are genuinely convention-dependent. (The reverse inclusion
  "noncrossing ⇒ NN-reachable at depth O(w)" *is* easy by river routing, but that
  does not give equality at the *same* R.)
- Does not touch the cyclic/periodic-y strip (Remark 6): circular noncrossing is a
  different predicate and the half-plane argument does not apply.
- Does not start #636, does not price any transfer matrix, does not edit
  `docs/STATUS.md` / `docs/ROADMAP.md` / `docs/RESEARCH-FRONTIER.md`, does not close
  tickets, does not merge.

## 7. Outcome taxonomy

```text
noncrossing claim (Theorem 1)   proof — all w, all R (planarization + Jordan)
class-set equality (Claim 2)    verified-small-width — w ≤ 6, R ≤ 4, exact,
                                brute-force-cross-checked (honest convention);
                                unchanged status from #675, cleaner conventions
#675's 28/66                    reproduced as the NN side of the delayed-closing /
                                full-frontier convention (28 exact at (5,4); 66 at
                                (6,4), 65 at (6,3) — off-by-one suspected); under
                                that convention the two rule's class sets differ
converse at cluster level       still false (witness re-verified)
```

## 8. Reproducibility

- Script: `scripts/p680_frontier_partition_enumeration.py` (stdlib only).
- Run: `python scripts/p680_frontier_partition_enumeration.py --out results/p680-noncrossing-lemma-20260913/raw`
- Raw outputs + summary: `results/p680-noncrossing-lemma-20260913/` (see its
  `REPORT.md`, `metadata.json`, `commands.txt`).
- Machine: assigned local Mac; Python 3.x stdlib; no Huawei budget; maximum width
  enumerated: 6.
