# #712 — torus rank-2 onset, theta/wedge cores, diamond 3L−1: independent literature check, 2026-09-12

Retrieval-only literature survey for [issue #712](https://github.com/LightChainr/Matching-One/issues/712) (parent #650, geometry context in draft #710). **Does not enter** `docs/STATUS.md`. No claim-ledger row. No census, no transfer matrix, no Huawei.

**Question asked:** is the onset mass / equality classification proved by #710 for NN site percolation on the honest torus — axis rank-2 onset `2L−1` with `L²` row-column crosses; diamond (periods `(L,L)`,`(L,−L)`) onset `3L−1` with `4L²` minimizers (full `2L`-line plus a straight `L−1` plug), on a two-cycle core (wedge or theta; dumbbell excluded) with systolic length `L(|a|+|b|)` — already a theorem in the wrapping, homology, or combinatorial-topology literature?

**Headline answer: no.** None of the four literature strands named in the ticket classifies smallest wrapping configurations. All four are either probability papers or, where topology appears, record admissibility and coexistence, never onset mass or a count of minimizers. The diamond `3L−1` + `4L²` classification does **not** appear in print (explicit tripwire answer below).

---

## Q1 — Wrapping/crossing minimal occupied sets of homology rank 2 on the square torus

**Verdict: probabilities only. No minimal-set classification anywhere in the strand.**

### Newman–Ziff 2000/2001 — PRIMARY_TEXT_READ

M. E. J. Newman, R. M. Ziff, "Fast Monte Carlo algorithm for site or bond percolation" and companion PRL; fetched full text of the algorithm paper via [ar5iv mirror of cond-mat/0005264](https://ar5iv.labs.arxiv.org/html/cond-mat/0005264).

The paper defines the five wrapping observables verbatim:

> "R_L^(h) and R_L^(v) are the probabilities of wrapping horizontally or vertically around the system respectively; R_L^(b) is the probability of wrapping around both directions simultaneously; R_L^(e) is the probability of wrapping around either direction; and R_L^(1) is the probability of wrapping around one direction but not the other."

with exact inter-relations (`R^(e) = 2R^(h) − R^(b)`, `R^(1) = R^(h) − R^(b)`) and the known Pinson limits (`R_∞^(h)(p_c) = 0.521058290`, `R_∞^(b)(p_c) = 0.351642855`) used as `p_c` markers.

What it does **not** do: no homology language at all; no enumeration or classification of minimal occupied sets. The only "how many sites" object is the *distribution of the first-wrapping occupation number* across random runs — a statistical object, not a combinatorial classification. It explicitly acknowledges topologically distinct wrap geometries without classifying them:

> "(Note that configurations which wrap around both directions are taken to include both those which wrap directly around the boundary conditions and 'spiral' configurations in which a cluster wraps around both directions before joining up.)"

So even the rank-2 (both-directions) wrap is present only as an event probability. The `2L−1` onset and the `L²`-cross count are nowhere.

### Pruessner–Moloney 2004 — PRIMARY_TEXT_READ

G. Pruessner, N. R. Moloney; fetched full text via [ar5iv mirror of cond-mat/0310361](https://ar5iv.labs.arxiv.org/html/cond-mat/0310361).

Computes `P((a,b), n, r)` — probability of `n` clusters with winding numbers `(a,b)` at aspect ratio `r` — and benchmarks against Pinson's CFT formulae. The winding-number bookkeeping is topological but explicitly heuristic:

> "On the torus, it is a topological fact that if a=0 (b=0) then the only path which is not homotopic to a point and does not intersect itself has b=1 (a=1). … If a≠0 and b≠0 then a and b must be relative prime if the path does not intersect itself."

> "Multiple, distinct clusters with the *same* winding number can coexist without intersecting. This, however, does not apply to a cluster with a cross topology: there can only be one such cluster on a torus. … winding clusters with incommensurable winding numbers cannot coexist."

> "The topological considerations in this paper are mainly technically motivated and rather heuristic. For rigorous proofs, we refer to the standard literature."

This is admissibility and coexistence of winding classes — the percolation-language shadow of the intersection form (see Q3) — **not** an onset/equality classification. The words "minimum", "onset", or any minimal-occupied-set enumeration are absent; the only "size" discussions are lattice sizes and aspect-ratio probability decay.

### Langlands–Pichet–Pouliot–Saint-Aubin 1992/1994 — PRIMARY_TEXT_READ

R. P. Langlands, C. Pichet, Ph. Pouliot, Y. Saint-Aubin, *J. Stat. Phys.* **67**, 553 (1992); fetched full PDF from [IAS edition](https://publications.ias.edu/sites/default/files/universality-ps.pdf).

Crossing probabilities on **rectangles and other planar finite geometries** (six 2D models), testing Cardy-style universality. Full-text grep: no "torus", no "periodic", no "winding" anywhere except two bibliography entries. The torus crossing taxonomy that Pinson later uses is *defined* by this paper's programme but the paper itself never touches wrapping. No minimal configurations.

### Pinson 1994 — ABSTRACT_ONLY

H. T. Pinson, "Critical percolation on the torus", *J. Stat. Phys.* **75**, 1167–1177 (1994), [doi:10.1007/BF02186762](https://link.springer.com/article/10.1007/BF02186762). Springer and ADS served bot-challenge pages; only the listing abstract fragment was readable ("We compute the various crossing probabilities defined by R. Langlands, P. Pouliot, and Y. Saint-Aubin for the critical percolation on …"). CFT evaluation of torus crossing probabilities per homology class — the same formulas Pruessner–Moloney benchmark and this repo's `predictions/p156` baseline uses. No abstract-level or quoted indication of minimal-set classification; the probability framing is unambiguous, and Pruessner–Moloney's full text (fetched) confirms her results are probability formulas. **Do not treat as a configuration classification.**

### Supporting: Akhunzhanov–Eserkepov–Tarasevich 2022 — PRIMARY_TEXT_READ

R. K. Akhunzhanov, A. V. Eserkepov, Y. Y. Tarasevich, *J. Phys. A* **55**, 204004 (2022), [arXiv:2204.01517](https://arxiv.org/abs/2204.01517), fetched full text via [ar5iv](https://ar5iv.labs.arxiv.org/html/2204.01517). Closest published object to an *exact wrapping polynomial* on the torus (site, square, one specified direction, L ≤ 12). It constrains coefficients only by an orbit-stabilizer divisibility invariant — "the quantity L²/gcd(i, L²) must be a divisor of c_i" — and never states `c_L = L`, never names the full-row minimizers, and does not treat two-direction wrapping at all ("wrapping probability along one direction", abstract). No `2L−1`, no tilted/diamond cell. So even the exact-enumeration end of the strand stops short of minimal-configuration statements.

## Q2 — Diamond / tilted / 45-degree tori; onset 3L−1 or 4L²

**Tripwire, explicit: NO. The diamond `3L−1` onset and the `4L²` minimizer count do not appear in print**, as far as this retrieval could establish (searched arXiv, Google, publisher indexes on 2026-09-12; fetched texts above all checked). Every torus percolation paper retrieved uses the axis-aligned square cell; the only boundary-condition variations are plane/cylinder/torus (Akhunzhanov, quoted above). No paper states `2L−1` (axis rank-2 onset) either. The #692-derived and #710-proved statements remain unique to this repository.

## Q3 — Theta vs wedge vs dumbbell as cycle-rank-2 cores; intersection form

**Verdict: the needed statements are standard algebraic topology, exactly as #710 uses them, and the percolation strand's coexistence rules are their heuristic echo.**

### Farb–Margalit, *A Primer on Mapping Class Groups*, Ch. 1 — PRIMARY_TEXT_READ

Fetched full PDF ([v5.0 copy](https://pagine.dm.unipi.it/~a019210/Farb%20Magalit_Primer%20on%20Teichmuller%20theory.pdf), Princeton UP 2012), Chapter 1 verified by text extraction.

- **Proposition 1.5** (§1.2.2): "The nontrivial homotopy classes of oriented simple closed curves in T² are in bijective correspondence with the set of primitive elements of π1(T²) ≈ Z²." I.e. an essential single cycle has `gcd(|a|,|b|) = 1` (or `(±1,0)`, `(0,±1)`).
- **Intersection formulas** (§1.2.3): "For two such homotopy classes (p,q) and (p′,q′), we have î((p,q),(p′,q′)) = pq′ − p′q and i((p,q),(p′,q′)) = |pq′ − p′q|," with invariance under the `SL(2,Z)` change-of-coordinates action.

Consequence used by #710 (not printed as a named corollary in the primer, but immediate): two **disjoint** essential simple closed curves on the torus have geometric intersection 0, hence `pq′ − p′q = 0`, hence are homologous. A pair of homology classes generating `H_1(T²) ≅ Z²` (symplectic pairing ±1) must be realized by cycles that **intersect** — so a rank-2 core's two essential cycles cannot live in disjoint pieces of the graph. The **dumbbell** (two disjoint loops joined by a connecting path) has its two cycles with algebraic intersection 0 whatever their classes, and if they are homologous the pair cannot span rank 2; the **wedge** (shared vertex) and **theta** (three shared arcs between two vertices) realize intersecting, spanning pairs. That is the intersection-form exclusion of the dumbbell, in the standard literature.

### Hatcher, *Algebraic Topology*, Ch. 0–1 — PRIMARY_TEXT_READ

Fetched PDF ([pi.math.cornell.edu/~hatcher/AT/ATch1.pdf](https://pi.math.cornell.edu/~hatcher/AT/ATch1.pdf)), text-verified. §1.1/§1.2: every connected graph deformation-retracts onto a maximal tree, and collapsing it gives a wedge of circles, one per edge outside the tree — cycle rank `= |E| − |V| + 1` (worked as Example 1.22 for the cube graph: twelve edges, maximal tree of seven edges ⇒ π1 free of rank five). Wedge, theta, and dumbbell are all homotopy equivalent to the wedge of two circles — so *homotopy alone cannot exclude the dumbbell*; the exclusion needs the embedding in the torus plus the intersection form, which is what the Farb–Margalit facts supply. (Hatcher's Ch. 1/2 text contains no "theta graph"/"dumbbell" graph-classification passage — checked; the ticket's hunch that the citation lives in algebraic topology rather than percolation is right, but it lives in the two sources above, not in a single "graph embeddings forbid disjoint essential cycles" theorem.)

### Echo in the percolation literature — [LIT]

Pruessner–Moloney (fetched, quoted above) state the same obstruction in percolation language: winding clusters with "incommensurable" winding numbers cannot coexist, and there can be only one cross-topology cluster. That is the probability paper's version of "disjoint essential cycles must be homologous". They explicitly defer rigidity: "we refer to the standard literature."

## Q4 — Systolic length `L(|a|+|b|)` for class `(a,b)` on the unit-grid torus

**Verdict: the continuous analogue is classical; the exact discrete `L1` statement was not found in print. It is, however, a two-line covering-space argument, not a gap.**

### Systolic-inequality survey — PRIMARY_TEXT_READ

T. Ikonen, D. Marti, N. Vikman, "Systolic inequalities for metric surfaces via filling minimality", [arXiv:2607.22290v1](https://arxiv.org/html/2607.22290v1) (submitted 2026-07-24), fetched full HTML. States the background exactly:

> "**Loewner's systolic inequality.** Let (𝕋,g) be a 2-dimensional torus with a Riemannian metric. Then, Area(𝕋,g) ≥ (√3/2) sys(𝕋,g)². Equality holds if and only if (𝕋,g) is isometric to the quotient of the Euclidean plane by some hexagonal lattice."

and, for the flat torus, `sys(ℝ²/ℤ², ‖·‖) = min_{z∈ℤ²∖{0}} ‖z‖`, with shortest length in an integer class given by the stable norm `‖z‖_st = min_x d(x, x+z)` (eqs. 4.4–4.5). **No discrete, graph, or Manhattan/L1 version appears** (the extremal metric is the *sup-norm* torus; the paper measures area by Hausdorff measure, so a grid graph is out of scope). [LIT] Loewner's result itself is unpublished (as the survey notes) — cite it through the survey.

### The discrete statement — [LIT] (reconstruction, no print source found)

On the `L×L` unit-grid torus, a closed walk with winding `(a,b)` lifts to a walk in the covering `Z²` from a vertex to its translate by `(aL, bL)`; every grid edge changes L1 displacement by exactly 1, so the walk has at least `|a|L + |b|L` edges. Attainability for the classes #710 needs (axes and `(1,±1)`-type diagonals, which are primitive) is by the straight/staircase path. Caveat worth keeping next to the proof: **simple** cycles on the torus are primitive classes (Farb–Margalit Prop. 1.5 above), so the equality classification for non-primitive `(a,b)` concerns multi-cycle cores, not a single simple cycle. Searched for the exact `L(|a|+|b|)` bound in the wrapping/first-passage/chemical-distance literature: not found in this form; nearest neighbours are the continuous stable-norm statements above. The `2L−1` axis onset *is* this bound at `(1,1)` combined with the row-column minimizer count — which is why it proves nothing less than #710's theorem.

---

## Tripwire summary (required by ticket)

1. **Three fetched primary texts:** exceeded — Newman–Ziff (ar5iv), Pruessner–Moloney (ar5iv), Akhunzhanov et al. (ar5iv), Langlands–Pichet–Pouliot–Saint-Aubin (IAS PDF), Farb–Margalit Ch. 1 (PDF), Hatcher Ch. 1 (PDF), systolic survey (arXiv HTML). All tagged above.
2. **Diamond `3L−1` + `4L²` in print: NO.** Axis `2L−1` in print: also NO. Nothing in the wrapping, homology, or combinatorial-topology literature classifies smallest wrapping configurations; the onset/equality classification of #710 is not a known theorem.
3. **Probability papers are not configuration classifications:** confirmed at the source for Newman–Ziff and Pruessner–Moloney (explicit full-text quotes above); Pinson tagged ABSTRACT_ONLY and used only as a probability source.

## Not established

- that no other literature exists outside the searched strands (this was a bounded retrieval, not an exhaustive survey);
- the exact discrete systolic statement in print (not found; elementary anyway);
- any change to #710's proofs, which this note neither strengthens nor weakens — it establishes that their content is novel as literature.
