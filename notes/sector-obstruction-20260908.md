# #658 sector obstruction: the two D-carrying wrapping amplitudes versus Λ_open − Λ_closed of one operator

**Ticket:** LightChainr/Matching-One#658. **Date:** 2026-09-08. **Mode:** reasoning /
exact-identity on objects that already exist. No transfer matrix was implemented; no
`docs/STATUS.md` entry; no ticket closed or merged. Every exact claim below was
re-verified this session by per-configuration enumeration on the committed geometries
(axis L=2,3,4; diamond L=2,3) against the committed Bernstein integers; the state-space
lemmas are verified up to width 6 and marked as such.

**Answer up front.**

> **Verdict: obstruction — but a sharper one than the ticket anticipated, and with a
> live residue.** The two D-carrying wrapping amplitudes *cannot* be Λ_open − Λ_closed
> of one periodic connectivity operator at finite L, for an exact reason that also kills
> the near-miss variants: **at finite L our two amplitudes are not eigenvalues of
> anything — they are integer configuration counts on a single finite probability
> space, while Jacobsen's Λ_open − Λ_closed is spectral data whose equality criterion
> exists only after the cylinder-length limit is taken first.** The two objects have
> different mathematical types, and no finite-(n,m) identity of the Jacobsen kind exists
> to bridge them. However, the exact finite-L identity that *does* hold is stronger and
> cleaner than the two-graph picture of PR #646: **M(p) = A^cross(p) − A^none(p) is
> already a difference of two homology-class event amplitudes of the black NN graph
> alone** (Fact 1 below). What obstructs the one-operator map is not the sector
> bookkeeping — it is the eigenvalue-versus-probability type mismatch, i.e. precisely
> the limits-first clause of Jacobsen's theorem.

---

## 0. Assets consumed (read first, not repeated)

- `notes/jacobsen-sector-retrieval-20260907.md` (PR #645): Jacobsen 2015 criterion,
  limits order, what is proved; Mertens–Ziff 2016 finite matching identity; the Q2 gap
  note on NN+NNN strip states (partially corrected below, §4).
- `notes/torus-wrapping-retrieval-20260908.md` (PR #654): MZ eq. (20) as the finite-L
  identity behind `D(C)`; MZ pairing (9)–(11), (19); Pinson attribution boundaries.
- `notes/probe635-sector-map-20260907.md` (PR #646): verdict A, both-same collapse.
- `notes/wrapping-type-census-l3l4-20260908.md` (PR #653) and
  `notes/wrapping-type-census-axis-L5-diamond-L4-20260908.md` (PR #657): five-cell
  support, both-two emptiness, axis L=5 joint tables.
- `results/exact_small_matching_polynomials.md`, `scripts/exact_matching_polynomial.py`:
  the Bernstein acceptance gate.

## 1. Exact facts established this session (all machine-checked)

Throughout, label a colouring's wrap class by the displacement-lattice data of its
clusters: `none` (no cluster wraps), `x`/`y` (some cluster wraps exactly one axis and
none wraps both), `spiral` (a cluster wraps both axes with rank-1 displacement lattice —
the #651 `both-same` case), `cross` (a cluster with rank-2 displacement lattice), and
`both-two` (two distinct clusters, one per axis). The rank-2/`cross` bit refines #646's
`both-same`, which lumps `cross` with `spiral`.

**Fact 1 (black-side collapse; exact, per configuration, all five committed geometries).**

```text
D(C) = 1{black has a cross (rank-2) wrapping cluster} − 1{black has no wrapping cluster}.
```

Verified configuration-by-configuration on axis L=2 (2⁴), L=3 (2⁹), L=4 (2¹⁶), diamond
L=2 (2⁸), L=3 (2¹⁸). No exception. Consequently, per occupation number k,

```text
a_k = #cross_k − #none_k          (black side only)
```

and the resulting `a` reproduces the committed Bernstein integers bit-for-bit, e.g.
axis L=3:

```text
a = [-1, -9, -36, -78, -90, -36, 36, 36, 9, 1]      (committed == recomputed)
```

and diamond L=3 (`-1, -18, -153, -816, -3060, -8568, -18438, -30528, -37638, -31640,
-13536, 3816, 9696, 6804, 2844, 804, 153, 18, 1`). This is a strict strengthening of PR
#646's verdict A: there the two amplitudes were `both-same` counts on two lattices
(black NN vs white NN+NNN); here the white side is eliminated entirely and the black
`both-same` count splits as `cross` minus `none` on one graph.

**Fact 2 (white label is a deterministic function of the black label; exact, all five
geometries).**

```text
black none   ⇒ white cross      (white wraps  ⇔  black ≠ cross)
black cross  ⇒ white none
black x      ⇒ white x
black y      ⇒ white y
black spiral ⇒ white spiral
```

No other cell has mass. Verified per configuration (axis L=4: 2¹⁶ configurations; every
black label maps to exactly one white label). The `none ↔ cross` exclusion and the
count-level one-axis pairing are MZ's torus pairing (Mertens–Ziff 2016, eqs. (9)–(11)
and (19); quoted at first hand in `notes/torus-wrapping-retrieval-20260908.md` §Q2);
the labelwise determinism of the x/y/spiral rows is this session's refinement and holds
exactly at every size checked. Fact 1 follows from Fact 2 by pure logic: white wraps
iff black ≠ cross, hence

```text
D = 1{black ≠ none} − 1{black ≠ cross} = 1{black cross} − 1{black none}.
```

**Fact 3 (the five-cell support is exactly these two Facts).** The census cells of PR
#653/#657 (`none×both-same` and `both-same×none` carrying D; the three diagonal
zero-D cells) are Fact 2's graph once `both-same` is split into `cross` and `spiral`:
the D-carrying cell `none×both-same` is black-none/white-cross, and `both-same×none`
is black-cross/white-none. Axis L=5 and diamond L=4 (PR #657) are consistent with both
Facts at every k, including `both-two` mass zero. (No claim beyond the verified sizes:
that `both-two` stays empty at all L remains open, as PR #646 already records.)

## 2. The finite-L dictionary that DOES hold

Combining Facts 1–2 with MZ eq. (20), at every finite L, exactly:

```text
M_L(p) = P_p(black has a rank-2 wrapping cluster) − P_p(black has no wrapping cluster)
       = A^cross_{G_NN}(p) − A^none_{G_NN}(p),
```

a difference of two homology-class event amplitudes **of the black NN graph alone**.
This is a genuine "two sectors of one graph" statement at finite L — but of the
*probability* kind, not the eigenvalue kind. Note what the black-only form does and
does not buy:

- It does **not** make the matching lattice dispensable. The identity `A^cross −
  A^none = R^x − R̂^x` is MZ's pairing (Fact 2); without the white half there is no
  theorem that the difference has a sign change near p_c. The one-graph amplitudes
  A^cross, A^none are individually ~C(N,k)-scale counts whose difference is small;
  the matching function's clean [−1,1] range and monotonicity come from the
  probability-difference reading, which MZ's identity certifies only through the
  two-graph equality.
- It does collapse the sector question: whatever "two sectors" the finite-L identity
  carries, they are `cross` and `none` — homology rank 2 and rank 0 of one graph —
  not "two lattices".

## 3. The obstruction (exact; kills the near-miss variants too)

Jacobsen's criterion, in his own setting (quoted and status-marked in
`notes/jacobsen-sector-retrieval-20260907.md` §1, read at first hand from
arXiv:1507.03027):

> `P_B(q,v)=0 ⇔ Λ_open = Λ_closed`, valid for a basis B of size n×m, with **n finite
> and m→∞**.

with the limit order explicit: the cylinder length m goes to infinity *first*, at
finite circumference n; only then an outer n → ∞ extrapolation. The two Λ's are the
leading eigenvalues of the s=0 block of **one** pTL transfer matrix — one graph, one
edge weight v, one vector space, the sectors being its invariant subspaces.

**Obstruction O1 (type mismatch; exact).** At finite L, the two amplitudes the ticket
names are

```text
a_k = #cross_k − #none_k     ∈ ℤ,
```

configuration counts on the finite probability space {0,1}^N; M_L(p) = Σ_k a_k
p^k(1−p)^{N−k} is a difference of two probabilities of events on *that same space*
(this is what MZ eq. (20) says). Λ_open − Λ_closed at finite (n,m) is a difference of
two leading eigenvalues of a weighted (fugacity-v) operator. No identity of the form
"probability difference = eigenvalue difference" exists at finite (n,m) in the pTL
setting, and none can: the eigenvalue equality in Jacobsen's theorem is a statement
about the *m → ∞* asymptotic slopes of the two sector free energies (exponential growth
rates), made exact by the intermediate-value-theorem crossing. At finite m the two
sides of his criterion are not even defined as equal-or-unequal — the theorem's
hypothesis is the limit. Therefore "map-holds at finite L" in the eigenvalue sense is
dead on arrival, not for a lack of cleverness but because the finite-L object that
would have to appear on the eigenvalue side does not exist.

**Obstruction O2 (one-operator ⇒ one-graph; exact).** Suppose the type mismatch is
waived by asking for the probabilistic shadow: can the two D-carrying amplitudes be
*sector amplitudes* (Pinson-style Z-class weights) of one connectivity operator? A
single connectivity transfer operator's state space carries the connectivity σ-algebra
of **one** graph. The candidate one graph would have to be the union G_U = NN ∪ NNN.
But:

- black-wraps (NN) is not a function of the G_U-connectivity state: two colourings
  with the same G_U cluster partition can differ in whether the *NN-only* subgraph
  wraps (a G_U-wrapping cluster may wrap only through diagonals — e.g. the diamond
  relay path, §4). So the black amplitude is not measurable on G_U's connectivity
  algebra;
- white-wraps (NN+NNN) is trivially measurable on G_U but *not* on G_NN; so G_NN fails
  symmetrically.

Hence no single graph's connectivity algebra carries both amplitudes as sector
observables. Carrying both requires the pair (G_NN, G_NN+NNN) — a two-graph (doubled)
state space on the shared occupation field. That doubled object is not "one operator"
in Jacobsen's sense; his theorem's hypothesis (one pTL operator, sectors = invariant
subspaces of its s=0 block) does not apply, and no published eigenvalue identity covers
the doubled operator. This is the same wall the Q1 search in PR #645 hit from the
literature side ("no source found that writes a signed combination 'primary wrapping −
matching-lattice wrapping' and identifies it with two transfer-matrix sectors"); here
it is shown to be structural, not just unpublished.

**Obstruction O3 (limits-first; exact, and the honest residual).** The only reading of
the ticket's map that survives O1/O2 is: *the two amplitudes become eigenvalue-type
data only after m → ∞, and the map holds in that limit if at all.* Formally: on the
n×m torus (both periods finite), take m → ∞ at fixed n; the sector free energies
f_open(n), f_closed(n) exist; Λ_open − Λ_closed = e^{−m f_open} − e^{−m f_closed}-type
comparisons become slope comparisons. Our M_L is defined with both periods finite; its
root p*_L is a probability crossing. Jacobsen's theorem does not say the probability
crossing equals the eigenvalue crossing at finite m — it cannot, by O1 — and MZ's
identity (which *is* exact at finite m) identifies M_L with a probability difference,
not with anything spectral. So the finite-L defect of the map is not a small
correction: at finite (n,m) the eigenvalue side of the dictionary is simply absent,
and the defect is the whole difference between "difference of two configuration-count
probabilities on {0,1}^N" and "difference of two asymptotic slopes". The near-miss
variant "maybe at finite L the two sector eigenvalues of some cleverly chosen operator
happen to reproduce the Bernstein integers" is killed by O1 (the Bernstein integers are
counts, not eigenvalues, and the count identity is already exact by MZ — there is
nothing spectral left for eigenvalues to reproduce); the variant "maybe open/closed
sectors of the black NN operator" is killed by O2's measurability test applied to the
white half, plus the observation that black-side `cross`/`none` are *events*, not
invariant subspaces: the black operator's sector decomposition by homology weights
(Pinson-style) sums over all homology classes, and `1{some cluster is rank-2}` is a
disjunction over sectors, not a sector.

**Consequence for #636 (pricing, per the ticket's framing).** If anyone builds a
connectivity transfer matrix for square site, its sectors will reproduce *wrapping
probabilities of one graph*; the matching function M(p) will relate to it only through
MZ's two-graph identity, which involves the NN+NNN graph on the same occupation field.
The exact one-graph finite-L dictionary of §2 is the correct place to state what the
map is; the eigenvalue identity of Jacobsen belongs to the m→∞ cylinder and cannot be
imported to finite L. "Limits-first" is not a technicality here; it is the whole gap.

## 4. State-space lemmas (supporting; verified to width 6, proof sketch given)

These settle a question PR #645 left open (its Q2: is the NN+NNN strip state space
still noncrossing / Catalan-counted?).

**Lemma A (no crossing frontier states; verified w ≤ 6, rows ≤ 4; proof sketch).**
For the square strip with NN+NNN site connectivity, the set partition induced on the
frontier sites by occupied-site connectivity is always **noncrossing**. Moreover the
reachable set equals the NN-only reachable set (28 states at w=5 over 4 rows and 66 at
w=6 over 3 rows, identical sets; both strictly inside the 42/132 noncrossing
partitions). Proof sketch: straight-line drawings of king-graph edges cross only at
face centres, between the two diagonals of one face; if both diagonals are active then
all four face sites are occupied, and the NN edges of the face merge the two paths —
so any potential geometric crossing heals into a merge, and a Jordan-curve argument
(occupancy path A from f_i to f_j with i<j−1, interleaved path B from f_{i'} to f_{j'}
with i<i'<j<j', both avoiding each other's endpoints) yields a contradiction. The
converse direction is genuinely false at the cluster level — diagonals DO create
non-NN frontier footprints (e.g. a component touching frontier columns {0,2} with
column 1 vacant: sites (0,r),(1,r−1),(2,r) occupied, verified witness) — they just
never create *crossing* ones. **Status: verified by exhaustive enumeration up to width
6; the proof sketch is not written out in full rigour. Marked: lemma-with-sketch.**

This corrects the strong form of PR #645's Q2 suspicion ("the frontier connectivity
relation admits crossing configurations, e.g. through a K4 on a single face"): the K4
on a face is exactly what *prevents* crossing states, because a K4's diagonals never
occur without the NN edges that merge them. The practical upshot for #636's cost model:
a NN+NNN strip transfer matrix, if built, lives on the same noncrossing class as the
NN one (same reachable partition set, verified w ≤ 6) — the cost blow-up feared in PR
#645 §Q2 is, on this evidence, not caused by crossing states. (The state counts may
still differ from NN through weights/weights-per-state and through the *number of
realizations per state*; the partition *classes* coincide.)

## 5. Quotations relied on (status as marked by the sources)

- Jacobsen 2015 (arXiv:1507.03027), as quoted in `notes/jacobsen-sector-retrieval-20260907.md` §1
  (read at first hand there): "P_B(q,v)=0 ⇔ Λ_open = Λ_closed, valid for a basis B of
  size n×m, with n finite and m→∞"; the proof is by intermediate value theorem
  ("our main result"); the n→∞ step is an extrapolation, not proved. O1/O3 lean on the
  limit order exactly as stated there.
- Mertens–Ziff 2016 (PRE 94, 062152; arXiv:1603.07289), eq. (20), as quoted in
  `notes/torus-wrapping-retrieval-20260908.md` §Q2: `M_L(p) = R^x_L(p) − R̂^x_L(1−p)`
  for x ∈ {c, b, e, h}, exact at every finite L, "the only contribution to the
  right-hand side is the cross-wrapping probabilities". Facts 1–3 and §2 lean on this.

## 6. What this does not do

- Does not implement or price a transfer matrix; does not start #636.
- Does not prove Lemma A in full rigour (enumerated to width 6 only) and does not prove
  `both-two` emptiness at all L (PR #646/#651 boundary untouched).
- Does not claim the root p*_L equals any eigenvalue-crossing root at finite L; per O1,
  no such equality exists.
- Does not edit `docs/STATUS.md`, close tickets, or merge.

## 7. One-line verdicts against the ticket's outcome taxonomy

```text
map-holds      no — killed by O1 (type mismatch: counts vs eigenvalues) at finite L
obstruction    yes — O1 + O2 (one operator ⇒ one graph; neither observable is
               measurable on a single graph's connectivity algebra) + O3 (the
               eigenvalue criterion exists only after m→∞); stated so the
               "finite-L sector eigenvalues reproduce the Bernstein integers"
               near-miss dies (the integers are already exactly reproduced by
               probability counts — there is no spectral residue left to match)
limits-first   the residual live reading: after m→∞ at fixed n, a one-graph
               cylinder has genuine open/closed sector eigenvalues; the map from
               OUR two amplitudes onto them passes through MZ's two-graph identity,
               which is not a one-operator object at any finite size. The finite-L
               defect is O(1)-type (whole-object absence), not a power correction.
```
