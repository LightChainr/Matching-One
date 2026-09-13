# #552 — R1 (15-state ordered serial algebra): novelty check + scoping decision

**Branch:** `retrieval/p552-r1-novelty-20260913`
**Ticket:** #552  ·  **Source:** issue #552 body; `docs/PUBLICATION-PORTFOLIO.md` (read from PR #548's open branch `docs/publication-portfolio-agent-20260901`, since the file is **not on `main`** — PR #548 is OPEN)

Retrieval mode: novelty judgement only. A negative search is **not** proof of originality (§5).

---

## 1. What the object is (from the ticket + portfolio doc)

The **15-state ordered serial algebra** is a finite semigroup/algebra arising from the project's
connectivity-reduction program: its subsemigroup inclusion lattice (#511), unit-group action orbits
(#512), element centralizers (#513), divisibility-class posets (#514), Green relations /
generalized inverses / power profiles / operator semigroup (#515–#519), partition-monoid
enumeration and idempotent-generated sectors (#493, #495), Burnside orbit certificates (#387),
Möbius inversion for four-terminal partitions (#448), HNF filtration to index 13 (#330). The
portfolio doc's R1 section ("Terminal algebra: reserve until it touches probability") rates it
**high assets, low paper closure** and lists four high-value next tasks — tasks 1–3 = the
probability-connection / obstruction track (Option A); task 4 = "compare the resulting semigroup
with known diagram/partition/planar algebras and finite semigroup classes" (= Option B, pure
finite-semigroup publication).

## 2. Asset-list freshness check (ticket asked: is the merged-asset list outdated?)

Spot-checked every PR number cited in the ticket against the repo:

| PR | title (verified) | state |
|----|------------------|-------|
| #511 | Exact serial subsemigroup inclusion lattice | MERGED |
| #512 | (unit-group action orbits) | MERGED |
| #513 | (element centralizer catalog) | MERGED |
| #514 | (divisibility class posets) | MERGED |
| #515–#519 | Green relations / generalized inverses / power profiles / semigroup structure / commuting census (operator semigroup) | MERGED |
| #493 | Enumerate all local serial partition monoids | MERGED |
| #495 | (partition monoid enumeration / idempotent-generated sectors) | MERGED |
| #387 | Add independent Burnside gadget-orbit certificate | MERGED |
| #448 | Exact: certify four-terminal partition Möbius inversion | MERGED |
| #330 | Exact: extend HNF filtration frontier to index 13 | MERGED |

**Conclusion:** the merged-asset list in the ticket is **current** (all cited PRs are merged). The
only caveat: `docs/PUBLICATION-PORTFOLIO.md` — which the ticket tells me to read for the R1 rating —
does **not** exist on `main`; it lives only on the still-OPEN PR #548. So the R1 "gate" wording is
proposed documentation, not yet merged policy. This does not change the recommendation below.

---

## 3. Retrieval matrix

| # | Query | Source(s) reached | Outcome | Grade |
|---|-------|-------------------|---------|-------|
| N1 | "ordered serial algebra" finite-dimensional Nakayama serial algebra classification | Nakayama/serial-algebra surveys (Grokipedia; AMS *B Proc.* 2018 on homological bounds; arXiv:2403.11359 "Nakayama algebras of small homological dimension"; arXiv:1710.04420) | "Serial algebra" = every fg module is a direct sum of uniserial modules; connected serial algebras = Nakayama algebras, classified by quiver (linear Aₙ or cyclic) + relations (paths of fixed length). Well-classified infinite family over a field — not a 15-element transformation semigroup. | ABSTRACT_ONLY / [LIT] |
| N2 | "15 element semigroup Green relations classification" | Green's relations references (Encyclopedia of Math; Clifford–Preston; Howie) | Green's relations are the standard tool; classification of *order-15* semigroups is not in the abstract references (the famous census reaches order 8). | [LIT] |
| N3 | "idempotent generated partition monoid diagram algebra" | "Enumeration of idempotents in planar diagram monoids" (ScienceDirect/J. Algebra 2019); Jones partition algebra; Halverson–Ram | Diagram/partition monoids are infinite families indexed by n; idempotent-generated sectors studied, but no match to a specific 15-element semigroup. | ABSTRACT_ONLY |
| N4 | "small semigroup enumeration database order 15" | SmallSemi GAP package (known to enumerate semigroups of order ≤ 8); Distler–Mitchell census | Exhaustive semigroup census stops at order 8; **order 15 is beyond the enumerated database**. | [LIT] |
| N5 | "15-state ordered serial algebra" exact phrase | Google Scholar / arXiv / web | **No hit** under this name or a clear synonym. | API_QUERY / WebSearch |

Grade legend: `PRIMARY_TEXT_READ` = full source body; `ABSTRACT_ONLY` = abstract/page only;
`[LIT]` = standard/secondary; `API_QUERY` = database/query.

---

## 4. Novelty conclusion (`cite-or-gap`)

**GAP (with a residual-verification caveat).** The targeted literature search did **not** surface a
prior classification of this exact 15-element ordered serial algebra under a different name. The
closest object classes are:
- *Serial / Nakayama algebras* — but these are algebras over a field, classified by quivers, an
  infinite family, not a single 15-element transformation semigroup. `[LIT]`
- *Diagram / partition monoids* — infinite families in n; idempotent-generated sectors studied but
  no 15-element match. `ABSTRACT_ONLY`
- *Small-semigroup census* — exhaustive only to order 8, so it cannot confirm or deny novelty at
  order 15. `[LIT]`

**Caveat (why this is not a certificate of novelty):** (i) I retrieved only the *asset-list
description* of the object, not its full multiplication table / Green-relation eggbox; (ii)
authoritative novelty confirmation requires either extending a SmallSemi-style enumeration to order
15 (not done) or a direct multiplication-table / Green-relations match against the semigroup census
(not done). So the honest status is **"no prior match found; novelty not yet positively certified."**

---

## 5. Decision memo (Option A vs Option B)

**Recommendation: pursue Option B (pure finite-semigroup publication) as the fast path, in parallel
with — not instead of — Option A.**

Rationale:
1. The novelty check (task 4 / Option B) returned **no prior match** under another name, so R1 is a
   plausible standalone finite-semigroup paper *pending* the residual verification in §4. Option B
   therefore closes fast and is currently unexplored, exactly as the ticket diagnoses.
2. Option A (probability connection / obstruction theorem) remains open and is the higher-impact
   path, but the portfolio doc itself gates R1 on it and explicitly lists task 4 (the Option-B
   comparison) as a legitimate high-value task. The two are not mutually exclusive; Option B can be
   scoped as its own manuscript while Option A continues.
3. The residual verification (GAP SmallSemi at order 15, or a multiplication-table match) is a
   *pre-submission* requirement, not a blocker for starting the Option-B draft. If that verification
   later finds a prior classification, Option B collapses and Option A becomes the only path — which
   the ticket already anticipates.

**Non-claims (per ticket):** this note does not assume the algebra is novel (§4 caveat); does not
claim the probability connection is impossible (Option A stays open); does not merge R1 into P2/P3.

**Next action:** (a) run the order-15 novelty verification (GAP SmallSemi / transformation-semigroup
enumeration or a multiplication-table equivalence check) before any Option-B submission; (b) only
after a clean verification, open the Draft PR skeleton under `docs/manuscripts/r1-serial-algebra/`
using the Option-B outline in the ticket.

---

## 6. Mandatory disclaimer

Novelty judgement only. Negative retrieval ≠ proof of originality. The order-15 verification is the
missing half of the novelty check and must be done before claiming independence. The portfolio-doc
R1 gate is on an open PR (#548) and is not yet repo policy on `main`.

*Full Matching-One repository CI has not been run for this commit.*
