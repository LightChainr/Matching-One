# Probe LIT — Literature retrieval for location-to-rate (standalone)

**Assign to:** the literature officer / any Agent that can fetch arXiv PDFs and published pages. This is **not** one of the three computational/proof probes. Do not run N=725. Do not enumerate tori. Do not implement chart transport.

**The user cannot fetch arXiv.** Every theorem you rely on must appear as a **primary quote** with PDF-stable coordinates (title, authors, year, venue, arXiv id, theorem/lemma/proposition number, page, displayed equation). HTML scrapes that garble numbering are not acceptable. Use the PDF.

---

## Standing

Same epistemic level as #601 / #592 / the #613 note. Theory notes only.

- Does **not** enter `docs/STATUS.md`.
- Does **not** close #276, #321, #613, #566, #592, #601, #606, #615.
- Does **not** mix with PR #602, #604, #606, #607, #615, or with the three probe PRs this packet launches.
- Negative answers are first-class results. “No named theorem for X” is often the deliverable.
- Do not re-summarise sources already quoted in `notes/p613-quantile-convergence-20260907.md` (DCT Thm 1.1 item 3, Grimmett–Li (1.3), DKS abstract, van den Berg counterexample sentence). **Copy those quotes once into a reuse table**, then go looking for what they do *not* contain.
- Do not invent theorems. Do not treat triangular-site theorems as square-site theorems. Do not treat Gate 1 “five widths” or #582 “five transitions” as literature.
- Mark every item **FOUND / NOT FOUND / FOUND-BUT-WRONG-LATTICE / FOUND-BUT-NOT-INVERSE-CDF**.

PR against **frontier** (`claude/matching-one-workspace-pwr5pv`) *or* a fresh `docs/literature-officer-…` branch off frontier. One PR, this topic only.

The rate probe (MID) will emit a list `BLOCKED_ON_LITERATURE_PROBE`. That list, if it exists by the time you start, is **additional** to the questions below, not a replacement. If it does not exist yet, do not wait; this document is self-contained.

---

## Why this exists, and what is already closed

#613 / PR #614 is a scoped corollary:

```text
Q_N(u) → p_c     uniformly on compact subsets of (0,1)
```

for honest periodic square-cell tori with `ℓ_N / log N → ∞`, under H1–H3. It claims **no rate**. DKS is the literature parent of the *location* argument (giant cycles / homological percolation on a torus). DCT supplies subcritical exponential decay, including a site remark. Grimmett–Li supply `p_u + p_c^* = 1`; (H3) `p_c + p_c^* = 1` on amenable Z² uses `p_c = p_u` as well. van den Berg 1981 has counterexamples to matching-sum-to-one outside a restricted class.

The #606 ledger splits the empirical `L^{-4}` root shift:

```text
F1: M'(p_c) ≍ L^{3/4}  ⇔  α_4 ≍ L^{-5/4}, ν = 4/3
    triangular: Smirnov / LSW (named). square: open.
F2: M(p_c) ≍ L^{-13/4}  (Q4, x=21/4) — operator-level conjecture.
Outcome J: F1 and F2 are not currently available in one rigorous model.
```

Your job is **not** to prove F1 on square site (that is a Fields-level problem). Your job is to answer, with quotes, the *rate-shaped* questions that Theorem L and the ledger leave open, including “this does not exist”.

---

## North-star questions (answer each; “not found” is allowed)

The rate probe is forbidden from searching. You are the only search. Every question is of the form: **is there a named theorem that does X, and if so what does it actually say, on which lattice, for which observable (crossing probability, magnetization, inverse-CDF, homological rank, …)?**

### Q1 — Inverse-CDF / quantile finite-size scaling

Is there a theorem that an *inverse* of a sharp-threshold family

```text
Q_N(u) := inf{ p : P_N(event_p) ≥ u }
```

inherits a window width from the event? In particular for percolation crossing events or for `{r_G ≥ 1}` / `{r_G = 2}` on a torus?

- FOUND must quote the exact observable. A theorem about `p_c(N)` (a single quantile, usually u=1/2 of a crossing event) is **not** automatically a theorem about nine interior quantiles.
- Distinguish: existence of *some* window (Friedgut–Kalai, Bourgain, Hatami, DKS sharp-threshold) vs the *percolation* window `N^{-3/4}` / `L^{-1/ν}`.

### Q2 — DKS and descendants: rates, not just location

Duncan–Kahle–Schweinhart, arXiv:2011.11903, AIHP 2025: giant cycles have a *sharp threshold*. Quote:

- what “sharp” means in their paper (window width? 1-arm? which d, i?);
- whether they give a rate for `P(im H_i ≠ 0)` in a `L^{-θ}` window;
- descendants (citations forward) that prove a *power* for homological percolation, if any.

Their d=2, i=1 plaquette model is bond percolation on the square torus. Square-*site* + matching is the #613 lift via H3. Do **not** claim novelty for that lift; do record whether anyone has written the site+matching version in print.

### Q3 — Window-width theorems vs percolation window

For each of the following, quote the window they actually give, the hypothesis class (boolean functions, monotone events, isoperimetry, …), and whether it can ever be as small as `L^{-1/ν}` without RSW/arm inputs:

- Friedgut–Kalai (threshold width `p / log n`);
- Bourgain’s sharp-threshold theorem;
- Hatami’s theorem;
- DKS’s own sharpness (if distinct from Q2);
- any later “sharp threshold ⇒ inverse-CDF window” lemma.

**NOT FOUND** for the last bullet is highly useful.

### Q4 — RSW / box-crossing: who has what lattice

A table, primary quotes, no blending:

| lattice / model | RSW? | source (Thm, page) | implies what for crossing *probabilities* | implies what for *inverse-CDF* (or NOT FOUND) |
|---|---|---|---|---|
| triangular site | | Smirnov / … | | |
| square bond | | | | |
| square site | | | | |
| FK-Ising / critical Ising | | Chelkak–Smirnov / … | | |
| homological torus events | | | | |

Square site is expected **open**. If you find a square-site RSW, it is a big deal: quote it, and flag it as in tension with the #606 ledger’s “square-site conformal invariance open”. Do not hide it.

### Q5 — Four-arm `α_4 = −5/4` and `ν = 4/3`

F1 is this. Quote LSW (which paper, which theorem, which lattice). Quote whatever is known for:

- triangular site (expected: proved);
- square bond;
- square site (expected: open);
- any “four-arm ⇒ `M' ≍ L^{3/4}` on a torus” lemma connecting *arm events* to the derivative of a homological observable. **This connecting lemma is likely NOT FOUND** and that is the actual hole between LSW and F1-on-a-torus.

### Q6 — Finite-size *corrections* on lattices that *do* have conformal invariance

Triangular site and FK-Ising have theorems the square site does not. What is *proved* (not simulated) for:

- `p_c(L) − p_c` (crossing-probability inverse, any u);
- crossing probability at `p_c` minus its limit;
- any quantile of a homological event;
- any `L^{-θ}` with a named `θ`.

Simulations (Ziff, Newman–Ziff, Jacobsen, …) are **not** theorems. If you cite Ziff 2006 (#566 already asked for a first-hand reading), say it is empirical and quote the observable. Do not launder a fit into F1.

### Q7 — Grimmett–Li (1.1) vs (1.3), primary PDF

The #613 note already warns: (1.1) `p_c(G)+p_c(G*)=1` is **not** their theorem; (1.3) `p_u(G)+p_c(G*)=1` is. Verify against the **published PDF** (RSA 65 (2024); arXiv:2205.02734). Quote:

- (1.1), (1.2), (1.3) verbatim;
- the Sykes–Essam “motivation” sentence;
- the amenability / Burton–Keane step the note uses to get (H3);
- companion *Hyperbolic site percolation* arXiv:2203.00981, what it actually proves.

Then van den Berg, JMP 22 (1981) 152–157, DOI 10.1063/1.524747: quote the **abstract** counterexample sentence and the class for which his derivation holds. The note already has a piece of this; you must check the PDF, not the note, for the class.

### Q8 — DCT site adaptation

DCT CMP 343 (2016) Thm 1.1 item 3 is quoted in #613. Verify against PDF v3 (or the journal PDF). Quote the site-percolation paragraph in §1.2 and the `[AB87]` pointer. Confirm Aizenman–Barsky CMP 108 (1987) 489–526 covers finite-range *site* percolation as the note claims. If the note over-claims, say so. This is verification, not new retrieval — still do it, because H2 is load-bearing.

### Q9 — Homological percolation after DKS

Search forward from DKS 2011.11903 for:

- rates / window widths for `im H_i(S) → H_i(T^d)`;
- site vs plaquette vs bond;
- matching / dual statements;
- anything about *quantiles* of a homological height / rank.

NOT FOUND for quantiles is expected and useful.

### Q10 — Does monotonicity + sharpness give a common `ω` for all interior u?

Is there a general lemma (even outside percolation): if `p ↦ P_N(p)` is a monotone family of CDFs with a sharp window of width `w_N` at every u in `[ε,1−ε]`, then the inverse-CDFs share that window, uniformly on `[ε,1−ε]`?

This may live in probability textbooks (lead-lag of quantiles, Bahadur representation, monotone rearrangement) rather than percolation papers. If you find it, quote hypotheses (uniform sharpness, no `u_N→0`). If you find only the single-u case, say so. The rate probe’s R5 (one exponent for the 9-vector is not implied by location) *wants* this lemma or its absence.

### Q11 — Self-dual / self-matching finite-size laws

On triangular site and square bond, `p_L^H = 1/2` exactly (repo). What is proved for the *shape* `Q_L(u)−1/2` or for crossing-probability profiles at self-duality? Cardy’s formula is a *limit* at `p_c`, not a finite-L expansion. Quote any finite-L expansion that is a theorem. NOT FOUND is useful: then F2≡0 does not give a rate for shape either.

### Q12 — Connecting Russo pivotals to four-arm on a torus

In-repo, `M'(p) = pivotal_primal(p) + pivotal_matching(1-p)` is exact. Literature: is there a theorem that the expected number of torus-wrapping pivotals at `p_c` is `≍ L^{3/4}` (or `L^{α_4+2}` / whatever the conversion is)? Quote the conversion between four-arm probability and pivotal count **on Z²**, then say what extra is needed to pass to a **torus** (periodic images, wrapping definition). This is the missing sentence between LSW and F1-as-the-ledger-writes-it.

### Q13 — What Theorem L cannot inherit from sharp-threshold boolean analysis

A negative search: any paper that claims “because `{r_G>0}` is sharp, `Q_N(u)−p_c = O(N^{-α})` for some α from boolean influences”. If found, check whether α is the percolation α. If not found, write NOT FOUND.

### Q14 — X / Twitter, only if primary papers are thin

If Q1–Q13 are all NOT FOUND for inverse-CDFs, a secondary sweep of recent percolation talks/notes is allowed. Do not let a tweet replace a PDF. If you use it, it is a pointer, not a quote.

---

## Method

For each Qi:

1. Search (arXiv, journal PDF, Crossref). Record the query.
2. Open the PDF. Quote ≤ 12 lines that actually answer the question.
3. Tag FOUND / NOT FOUND / FOUND-BUT-WRONG-LATTICE / FOUND-BUT-NOT-INVERSE-CDF.
4. One sentence: what the rate probe is then allowed to write as T vs C vs blocked.

Do not paraphrase a theorem you have not quoted. Do not merge two papers into one “standard fact”. Do not cite Wikipedia.

Already-quoted sources (reuse table, then *stop re-summarising*):

- Duminil-Copin–Tassion, CMP 343 (2016), arXiv:1502.03050
- Aizenman–Barsky, CMP 108 (1987)
- Menshikov 1986
- Grimmett–Li, RSA 65 (2024), arXiv:2205.02734
- Grimmett–Li hyperbolic companion, arXiv:2203.00981
- van den Berg, JMP 22 (1981)
- Duncan–Kahle–Schweinhart, arXiv:2011.11903, AIHP 2025
- Burton–Keane 1989 (for `p_c=p_u` on amenable lattices)
- Smirnov / LSW as named in the #606 ledger (you must pin *which* LSW papers)

---

## What you will not do

- Prove F1 on square site, or “sketch” a proof.
- Fit `ν` or `α_4` to #582 / N=725.
- Discuss affine charts, 1.55, the 4% residual, weighting (HIGH probe).
- Re-enumerate L=3,4 (exact-controls).
- Close #276 / #321 / #613. Location is already a scoped corollary; you are not upgrading it to a rate by quotation density.
- Mix this PR with #602 / #615 / the three probes.
- STATUS.

---

## Deliverables

```text
notes/literature-officer-quantile-rate-YYYYMMDD.md
  — reuse table of already-quoted H1–H3 / DKS / van den Berg
  — Q1–Q14, each with tag, quote or NOT FOUND, and T/C/blocked implication
  — a one-page “rate probe may now write” sheet:
        allowed T sentences
        allowed C sentences (with the extra named)
        still blocked
notes/literature-officer-quantile-rate-biblio-YYYYMMDD.md
  — full bibliographic records, PDF-stable, no broken arXiv HTML
```

No scripts, unless you need a tiny bib helper. No JSON of productions.

---

## Interface

- Rate probe’s R10 questions are *added* to Q1–Q14 if they arrive first. If they duplicate a Qi, answer once and cross-reference.
- HIGH probe should never need you. If they ask for a paper on Fieller / Kriegeskorte / projective inference, that is a **different** literature request (#579/#602 already touched it); do not expand this PR into that topic.
- Exact-controls probe should never need you.

If you find a square-site RSW or a torus four-arm theorem, flag it on #276 and #321 **without closing them**, in one comment, quoting the theorem. That is the only issue-comment this probe is allowed.
