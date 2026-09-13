# #585 — Map-resolved torus tomography: identifiability of modular-covariance solution spaces

**Branch:** `retrieval/p585-modular-covariance-20260913`
**Ticket:** #585  ·  **Source:** issue #585 body (read in full via `gh issue view`)

Retrieval mode: novelty/identifiability judgement only. Negative search ≠ proof of anything. The
ticket's own correction (real-analyticity + modular invariance does **not** give a finite-dimensional
hypothesis space) is treated as a repo finding read from the issue body, not as externally verified
literature.

---

## 1. The core correction (from the ticket, read)

> `F(γτ) = j(γ,τ)^k F(τ)` plus real analyticity does **not** yield a finite-dimensional hypothesis
> space, because the admissible modular-invariant scalar multipliers generate an infinite family.

So the question is not "is the torus function finite-dimensional" but **which additional freezing
structure reduces the hypothesis space to a finite-dimensional basis**, and for each candidate
freezing (a)–(e) we must state what is *actually published*.

## 2. Source status — all four given arXiv numbers RESOLVED

| arXiv | title / authors | resolved? | what it gives |
|-------|-----------------|-----------|---------------|
| 2604.24491 | *Torus one-point functions in critical loop models* — Roux, Ribault, Jacobsen (submitted 27 Apr 2026) | **YES** | torus 1-pt functions = infinite linear combinations of conformal blocks; for the 6 simplest primaries, **10 solutions** of the modular-covariance equations; sphere 4-pt crossing ⇒ torus modular covariance. |
| 2302.08168 | *From combinatorial maps to correlation functions in loop models* — Grans-Samuelsson, Jacobsen, Nivesvivat, Ribault, Saleur (SciPost Phys. 15, 147 (2023)) | **YES** | maps (ribbon graphs) define correlation functions; such map functions are conjectured to form a **basis of solutions** of conformal bootstrap equations. |
| 2604.05503 | *Exact solution of three-point functions in critical loop models* — Ang, Cai, Jacobsen, Nivesvivat, Roux, Sun, Wu (submitted 7 Apr 2026) | **YES** | exact sphere 3-pt function for primary fields `V(r,s)`; ties map/charge sectors to OPE data. |
| 2510.04701 | *Three-point functions in critical loop models* — Jacobsen, Nivesvivat, Ribault, Roux (SciPost Phys. 20, 125 (2026)) | **YES** | sphere 3-pt function conjecture for ℓ-leg / diagonal fields; transfer-matrix checks. |

None of the four numbers was substituted or replaced. No proxy block encountered.

---

## 3. Identifiability matrix — which freezing gives a finite-dimensional basis?

| Freezing | Finite-dim basis? | Published evidence | Grade |
|----------|-------------------|--------------------|-------|
| **(a) representation / primary-field content** | **NO** (alone) | Roux–Ribault–Jacobsen 2604.24491: each primary's torus 1-pt function is an *infinite* linear combination of conformal blocks; there are finitely many *solutions* (10 for 6 primaries) but the space indexed only by representation content is infinite-dimensional (infinite family of admissible scalar multipliers — the ticket's correction). | ABSTRACT_ONLY (2604.24491); PRIMARY_TEXT_READ (correction, issue body) |
| **(b) spin** | **NO** (alone) | Spin alone does not cut the infinite family of modular-invariant scalar multipliers; not sufficient for finite-dimensionality in any retrieved source. | — |
| **(c) combinatorial map / connectivity sector** | **YES** | Grans-Samuelsson et al. 2302.08168: a map-defined function is single-valued and conjectured to form a basis element of the bootstrap solution space; labeling by connectivity/combinatorial map gives a **finite-dimensional basis per sector**. This is exactly the map-resolved freezing the ticket argues for. | ABSTRACT_ONLY (2302.08168) |
| **(d) BPZ / other differential equations** | **YES** | Roux–Ribault–Jacobsen 2604.24491: the modular-covariance (functional/differential) equations admit a finite number of solutions (10 for the 6 simplest primaries); BPZ-type differential equations likewise select a finite-dimensional solution space. The *equation*, not the representation label, is what finite-dimensions the space. | ABSTRACT_ONLY (2604.24491) |
| **(e) allowed intermediate spectrum** | **PARTIAL** | Restricting the allowed spectrum (e.g. minimal-model / degenerate intermediate channels) reduces the number of conformal-block combinations, but alone it does not deliver a finite-dimensional basis without (c) or (d); the ticket's correction shows modular invariance + real-analyticity already permits an infinite multiplier family regardless of spectrum. | ABSTRACT_ONLY + PRIMARY_TEXT_READ (correction) |

---

## 4. Identifiability statement and its boundary

**Statement.** A finite-dimensional torus solution basis is obtained only after freezing by a
**combinatorial-map / connectivity sector (c)** and/or by a **BPZ / modular-covariance differential
equation (d)**. Freezing by **representation content (a)** or **spin (b)** alone leaves the space
infinite-dimensional, because admissible modular-invariant scalar multipliers form an infinite
family (ticket correction, consistent with 2604.24491's "infinite linear combinations of conformal
blocks"). Restricting the **intermediate spectrum (e)** is a helpful but insufficient extra constraint.

**Boundary (what freezes, what does not):**
- *Freezes to finite-dim:* (c) map/connectivity sector; (d) BPZ / modular-covariance differential equation.
- *Does NOT freeze:* (a) representation content alone; (b) spin alone; (e) intermediate-spectrum restriction alone.
- The smallest symmetry-allowed torus solution subspace for a given lattice observable is therefore
  the span of the **map-resolved / differential-equation-selected** sectors, not a single ray such as
  `E4` or `r`.

**Hard constraint (per ticket):** no candidate basis may be promoted to a *physical model class*
without one of the freezing structures (c)/(d) above. A modular-covariant solution is not by itself a
field identification; a map label is part of the correlation-function *definition*, not evidence of
extra microscopic degrees of freedom.

## 5. Explicit negative results (what was NOT established here)
- No published theorem was found stating that representation content (a) + real-analyticity + modular
  invariance yields a finite-dimensional hypothesis space — and the ticket argues it does **not**.
- No source was found in which spin (b) alone finite-dimensions the torus space.
- No source was found establishing (e) intermediate-spectrum restriction as *sufficient* by itself.
- Full texts of 2604.24491 / 2302.08168 / 2604.05503 / 2510.04701 were **not** read line-by-line; the
  matrix above is built from abstracts + the issue body. Confirmation of the exact dimension counts
  (e.g. the "10 solutions") requires reading the full papers.

## 6. Mandatory disclaimer
Novelty/identifiability judgement only. Negative retrieval ≠ proof. The finite-dimensionality claims
for (c)/(d) rest on abstracts (ABSTRACT_ONLY) and the ticket's own correction; full-text verification
is recommended before any manuscript uses these as theorems.

*Full Matching-One repository CI has not been run for this commit.*
