# #635 identity write-up — is M(p) the difference of two topological-sector amplitudes?

**Ticket:** LightChainr/Matching-One#660 (write-up deliverable of #635). **Date:** 2026-09-08.
**Mode:** write-up against wrapping-type tables that already exist. No new census. No
`docs/STATUS.md` edit; no ticket closed or merged.

This note consumes, and does not re-enumerate:

- **PR #646** — probe635 sector map, verdict A, degenerate both-same sectors
  (`notes/probe635-sector-map-20260907.md`, `results/probe635-sector-map/latest.json`).
- **PR #653** — #632 bond-lab repair; wrapping-type 4×4 tables at axis L=3,4 and diamond L=3
  (`notes/wrapping-type-census-l3l4-20260908.md`, `results/wrapping-type-census/`).
- **PR #657** — axis L=5 (`2^25`) and diamond L=4 (`2^32`) censuses; Bernstein tripwire vs
  PR #649 pass; **#651 A-continues** (`both-two` mass 0 at every k)
  (`notes/wrapping-type-census-axis-L5-diamond-L4-20260908.md`).
- **PR #654** — wrapping literature retrieval (preferred over #656 for this note)
  (`notes/torus-wrapping-retrieval-20260908.md`).
- **PR #645** — Jacobsen retrieval (#637)
  (`notes/jacobsen-sector-retrieval-20260907.md`).

All exact claims below use integers / `Fraction` only. The tripwire integers were
re-run this session from `scripts/exact_matching_polynomial.py --geometry axis --L 3`
and reproduce bit-for-bit.

---

## 1. Sector labels on finite-L configurations

The #635 question is whether

```text
D(C) = 1{black NN wraps} − 1{white NN+NNN wraps}
M(p) = Σ_k a_k p^k (1−p)^(N−k)
```

decomposes as a difference of two topological-sector amplitudes, **exact at finite L**.
To make that decidable at finite L, the sector labels must be defined combinatorially,
on configurations, not on asymptotic objects.

**Five-name labels (PR #646).** For a configuration C on the L×L torus, compute the
displacement lattice of each occupied cluster (greedy spanning-forest lift into Z²):

| label | definition |
|---|---|
| `none` | no cluster has nontrivial winding |
| `x` | some cluster winds exactly the x-axis (up to sign) |
| `y` | some cluster winds exactly the y-axis (up to sign) |
| `both-same` | a **single** cluster whose displacement lattice has rank 2 in H₁(T²;Q) |
| `both-two` | two distinct clusters, one winding each axis |

Each configuration carries a black label (primal, NN) and a white label (matching,
NN+NNN). The candidate map of #646 is

```text
M(p) = Σ_k [ #(both-same | black, |C|=k) − #(both-same | white, |C|=k) ] · p^k (1−p)^(N−k)
```

i.e. the two sector amplitudes are the black- and white-winding indicators restricted
to the `both-same` homology class.

**Coarse 4×4 labels (PR #653).** The wrapping-type census classifier
(`scripts/wrapping_type_census.py`, `torus_homology.wrapping_channels`) uses
`neither / dir0 / dir1 / both`, where `both` is `dir0 ∧ dir1` and therefore **lumps
rank-2 cross-wrapping with rank-1 spirals**. The correspondence:

| #646 5-name | coarse 4×4 | Mertens–Ziff type (Newman–Ziff (13)–(15) language) |
|---|---|---|
| `none` | `neither` | no wrap |
| `x` | `dir0` (or `dir1`, orientation convention) | R^{(1)} / R^{(h)} single-axis |
| `y` | `dir1` (or `dir0`) | R^{(1)} / R^{(h)} single-axis |
| `both-same` | `both` (cross **and** spiral component) | R^{(b)}: both-axes wrap |
| `both-two` | folded into `both` | two single-axis wraps by distinct clusters |

The fold `both-same ∪ both-two → both` is exactly the distinction MZ do not make in
wrapping-type language; it is the one structural question the coarse 4×4 cannot see
and the 5-name labels can (#651's `both-two` mass).

**Mertens–Ziff cross vs spiral.** The retrieval in PR #654 fixes the identification:
MZ (2016) prove, at every finite L on the torus,

```text
M_L(p) := N_L(p) − N̂_L(1−p) − L² χ(p) = R^x_L(p) − R̂^x_L(1−p),  x ∈ {c, b, e, h}
```

with the **only contribution to the right-hand side coming from the cross-wrapping
probabilities** — all other wrapping types cancel by the MZ pairing (no black wrap ⇒
exactly one white cross; single-direction wraps pair; cross-wrapping is exclusive;
spiraling counts match). The repository observable (x = `either`) therefore equals the
cross difference by the published theorem. In 5-name language this predicts: D mass
only on cells where exactly one side is `both-same` (cross) and the other is `none`;
`x×x`, `y×y`, `both-same×both-same` (spiral pairings) carry D = 0.

## 2. The candidate identity, written so it can fail

The identity under test at finite L, in committed-integer form:

```text
a_k  =  #(black ∈ both-same, |C|=k)  −  #(white ∈ both-same, |C|=k)
```

Equivalently, as a two-amplitude difference:

```text
M(p) = A_black^{both-same}(p) − A_white^{both-same}(p)
A_black^{both-same}(p) = Σ_k #(black ∈ both-same, |C|=k) · p^k (1−p)^(N−k)
```

**Acceptance gate.** The committed Bernstein integers at axis L=3 (N=9), re-run this
session and reproduced bit-for-bit:

```text
bernstein_counts: [-1, -9, -36, -78, -90, -36, 36, 36, 9, 1]
polynomial:       -4p^9 + 18p^8 - 18p^7 + 6p^3 - 1
physical root:     0.5865114551126756356545589766069017348243
```

Any decomposition that does not reproduce these integers exactly (Fraction arithmetic,
no floats) is unused. The run also re-confirms the #628 §8 anchor by exact arithmetic:
M(1/2) = −21/64, so u_c = (1 + M(1/2))/2 = **43/128** at L=3.

PR #646's script, run this session at axis L=3, gives the two amplitude rows and their
sum:

```text
reference bernstein: [-1, -9, -36, -78, -90, -36, 36, 36, 9, 1]
  A_white (m-both-same): [-1, -9, -36, -78, -90, -45,  0,  0, 0, 0]
  A_black (p-both-same): [ 0,  0,   0,   0,   0,   9, 36, 36, 9, 1]
  sum:                   [-1, -9, -36, -78, -90, -36, 36, 36, 9, 1]   ← exact
```

with the independent wrap-label diagnostic path (union-find, not the
displacement-lattice rank path) agreeing position-by-position. **The identity holds
exactly at finite L on every size the enumerations reach — verdict A, in the
wrapping-form.**

## 3. Verdict A / B / C

**Verdict: A at finite L, with a degeneracy caveat that is now evidence-backed, not
speculative.**

- **A (exact at finite L)** — attained in the wrapping-form on all tested sizes:
  axis L=3 (2⁹), L=4 (2¹⁶), L=5 (2²⁵); diamond L=2 (2⁸), L=3, L=4 (2³²). The L=3
  integers are reproduced bit-for-bit (§2); the larger sizes pass the same tripwire
  against the PR #649 committed Bernstein rungs (PR #657, both kernels on diamond L=4).
- **Not B** — no finite-L defect exists at any reachable size; there is nothing to
  exhibit or scale.
- **Not C** — no obstruction: the MZ theorem (PR #654, their (20)) proves the
  difference-of-wrapping-amplitudes form at every finite L on a torus, so the
  wrapping-form identity cannot be false.

**The degeneracy caveat, now closed through axis L=5 / diamond L=4.** At axis L=3, the
per-sector rows above show the D-carrying cells are exclusively
`none × both-same` (D=−1) and `both-same × none` (D=+1); `x`, `y` carry D=0 wherever
they occur, and `both-two` never occurs with D≠0. The structural worry (PR #646 §3)
was that this collapse into a single `both-same` label measured twice might be a
small-L accident. It is not, on the last independently checkable axis rung:

- **Axis L=5** (2²⁵ configs, PR #657): `both-two` mass **0 at every k**; support is
  exactly the five cells `none×both-same` (−1), `both-same×none` (+1), `x×x` (0),
  `y×y` (0), `both-same×both-same` (0). Coarse 4×4 is still the five Mertens–Ziff
  cells; no sixth cell.
- **Diamond L=4** (2³² configs, K1 and K2 kernels, identical 5×5 tables): same five
  cells, `both-two` mass 0. **A-continues.**

Axis L=3 per-k D detail (PR #653 §3, exact): the sign change k=5→6 is two-term
cancellation — k=5 has `neither×both`=45 (D=−1) and `both×neither`=9 (D=+1), net −36;
k=6 has `both×neither`=36 plus a 6-count `both×both` (D=0), net +36. In 5-name terms
this is the same statement: D lives only on the exclusive cross mismatch.

## 4. What axis L=5 / diamond L=4 changed — and what they cannot change

**Changed.** They upgraded the degeneracy caveat from "observation at L ≤ 4" to
"A-continues at the enumeration frontier": `both-two` has zero mass at every k on
2²⁵ and 2³² configurations, and the coarse 4×4 support is unchanged at five MZ cells
on both geometries. They also confirm `dir0 = dir1` at every k on both geometries.

**They cannot change: the Jacobsen map (PR #645).** The published object
(arXiv:1507.03027) is a statement about the **leading eigenvalues of two topological
sectors of one periodic Temperley–Lieb transfer matrix**, with **m → ∞ taken first**
(semi-infinite cylinder), the criterion `P_B = 0 ⇔ Λ_open = Λ_closed` proved by
intermediate value theorem at finite n, and the n → ∞ extrapolation conjectural
(tower Δ_k = 2(k+1) marked as a conjecture in the source). Three specific gaps remain
between our finite-L verdict A and that map, none of which more enumeration can close:

1. **Order of limits.** Our identity is exact at finite L. Jacobsen's criterion is
   proved after taking the cylinder length to infinity. A finite-L amplitude identity
   does not imply equality of leading eigenvalues of the two sector transfer matrices;
   that is an asymptotic statement requiring its own limit.
2. **One operator, two sectors vs two lattices.** Jacobsen's sectors are open/closed
   blocks of **one** pTL operator. Our two amplitudes are black-NN on the primal
   lattice and white-NN+NNN on the **matching** lattice — two different edge sets,
   and the white indicator is not a sector of the primal transfer matrix. Whether the
   matching-lattice content is exactly the open/closed sector distinction is the crux
   #635 named and is untouched by census evidence. (PR #645 Q2: no published
   boundary-state theory for NN+NNN strip connectivity exists; the frontier relation
   admits crossings, so Catalan counting should not be assumed to survive.)
3. **Novelty boundary.** The published nearest neighbor is Mertens–Ziff 2016 itself:
   their identity **is** `M(p)` as a difference of wrapping amplitudes, exact at
   finite L — the pairing is published (PR #654). What is **not** published (PR #645
   Q1) is the sector reading: two sectors of *one* pTL operator vs a signed difference
   on *two* lattices. Any claim from this line must be phrased as the sector
   decomposition, with citation obligations to both Mertens–Ziff 2016 and
   Jacobsen 2015.

So: L=5 / L=4 strengthen A as a finite statement; they are silent on the eigenvalue
identification, which is a limit statement about a different object.

## 5. Boundary: this note does not fund #636

Nothing in this note promotes #636 (transfer-matrix machinery). The honest statement
of position:

- The wrapping-form of verdict A is **published** (Mertens–Ziff 2016). Building
  machinery to re-derive it buys nothing.
- The unpublished piece — the sector reading — is a *conceptual identification*
  (matching-lattice content ↔ open/closed sectors of one pTL operator), and the two
  census facts that would de-risk it (does `both-two` ever carry mass; does a sixth
  4×4 cell ever appear) stayed at zero through the enumeration frontier. Zero mass at
  L=5/L=4 is evidence of degeneracy, not a theorem for all L, and degeneracy, if
  exact, makes the "two sectors" degenerate into one label measured twice — which
  would make #636's sector machinery carry **less** structure than its cost suggests,
  not more.
- If machinery is ever worth building, the reason would be the **crossing frontier
  state space** of the matching side (NN+NNN admits crossing connectivities; PR #645
  Q2 gap, #638 is the correct instrument to measure it) — a measurement question for
  #638 first, a funding question for #636 after. That recommendation lives here, not
  in `docs/STATUS.md`, and does not promote anything.

**Tripwire standing rule:** any future claimed decomposition must reproduce
`[-1, -9, -36, -78, -90, -36, 36, 36, 9, 1]` at axis L=3 bit-for-bit, integers /
Fraction only, or it is unused.

---

**Claim boundary.** Finite identity, not a threshold theorem, not a low-dimensional
TM claim. The 2015 pc interval is superseded in the literature (PR #645 §2: Yang–Zhou
2024, Jacobsen Reply 2024) — recorded, not acted on here. Do not cite #628's 118133
bond `dual_fail` as physics (repaired compound artifact, PR #653/#646; the corrected
bond census has `dual_fail = 0` and exact `r_b + r_w = 2`). `M(p) + M(1−p) ≠ 0` stands
on the site side (complement-transpose of the 4×4 fails as soon as single-direction
wrap appears — the generating fact behind the non-self-complementing F).

Related: #635, #637, #640, #642, #646, #651, #653, #654, #657, #636 (not funded).
