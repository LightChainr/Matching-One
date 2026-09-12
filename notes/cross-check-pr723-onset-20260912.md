# Cross-check of PR #723 (torus rank-2 onset literature) — Round B, 2026-09-12

Peer check of [PR #723](https://github.com/LightChainr/Matching-One/pull/723) (`notes/lit-torus-onset-core-20260912.md`), per [issue #726](https://github.com/LightChainr/Matching-One/issues/726) (parent #650). I did **not** write #723. Retrieval-only; does not enter `docs/STATUS.md`; no claim-ledger row; no merge.

**What this check independently re-fetched (2026-09-12):**

1. **Newman–Ziff wrapping definitions** — full text re-fetched via ar5iv mirror of [cond-mat/0005264](https://ar5iv.labs.arxiv.org/html/cond-mat/0005264). PRIMARY_TEXT_READ (this check).
2. **Farb–Margalit Ch. 1** — full PDF re-fetched ([the same v5.0 copy #723 used](https://pagine.dm.unipi.it/~a019210/Farb%20Magalit_Primer%20on%20Teichmuller%20theory.pdf)); Prop 1.5 and §1.2.3 extracted and compared verbatim. PRIMARY_TEXT_READ (this check).
3. **Tripwire searches** for the diamond `3L−1` onset / `4L²` minimizer count (arXiv, web indexes). Confirmed absent.
4. **Spot-check of the systolic-survey citation** — arXiv [2607.22290](https://arxiv.org/abs/2607.22290) abstract page: exists, authors Ikonen/Marti/Vikman, submitted 2026-07-24, torus systolic inequality content as #723 describes. Abstract-level this check; #723's full-text quotes accepted.

**Headline: #723 survives the peer check.** Every claim I could test against a primary text tested true; one minor bibliographic correction (title ↔ arXiv-ID pairing for Newman–Ziff); nothing rises to DROP. The tripwire answer (diamond `3L−1` + `4L²` not in print) is confirmed.

---

## Required verifications

### Newman–Ziff wrapping definitions — re-fetched, verified verbatim

The quoted definitions of `R_L^(h)`, `R_L^(v)`, `R_L^(b)`, `R_L^(e)`, `R_L^(1)` appear word-for-word in the fetched text ("…are the probabilities of wrapping horizontally or vertically… wrapping around both directions simultaneously… either direction… one direction but not the other"), as do the "spiral configurations" parenthetical and the relations `R^(e) = R^(h)+R^(v)−R^(b) = 2R^(h)−R^(b)` and `R^(1) = R^(h)−R^(b) = R^(e)−R^(h) = ½(R^(e)−R^(b))` (with the monotonicity inequalities). The Pinson limits `R_∞^(h)(p_c) = 0.521058290` and `R_∞^(b)(p_c) = 0.351642855` are quoted correctly (the text also gives `R_∞^(e) = 0.690473725`, `R_∞^(1) = 0.169415435`, which #723 does not need).

Confirmed also: **no minimal-occupied-set classification anywhere** — the only "how many sites" object is the first-wrapping occupation number across random runs (a distribution), and no homology language appears. #723's characterization is accurate.

**One correction (bibliographic, not substantive):** the text served at `cond-mat/0005264` is the **PRL** — "Efficient Monte Carlo algorithm and high-precision results for percolation," *Phys. Rev. Lett.* **85**, 4104–4107 (2000). The title #723 prints, "Fast Monte Carlo algorithm for site or bond percolation," belongs to the companion **PRE** paper (*Phys. Rev. E* **64**, 016706 (2001), cond-mat/0005397). All the wrapping definitions quoted by #723 are present in the PRL that was actually fetched, so every substantive claim stands; only the title↔link pairing needs relabeling.

### Farb–Margalit Prop 1.5 — re-fetched, verified verbatim

Extracted from the PDF (§1.2.2, p. 28 of the v5.0 copy):

> "**Proposition 1.5** The nontrivial homotopy classes of oriented simple closed curves in T² are in bijective correspondence with the set of primitive elements of π1(T²) ≈ Z²."

followed by: "An element (p,q) of Z² is primitive if and only if (p,q) = (0,±1), (p,q) = (±1,0), or gcd(p,q) = 1." — exactly as #723 states.

The §1.2.3 intersection formulas also match verbatim: "For two such homotopy classes (p,q) and (p′,q′), we have î((p,q),(p′,q′)) = pq′ − p′q and i((p,q),(p′,q′)) = |pq′ − p′q|," with the invariance argument via A ∈ SL(2,Z) mapping a primitive class to (1,0). #723's derived consequence — disjoint essential simple closed curves have i = 0 hence pq′ = p′q, so (with both primitive) they are homologous; a spanning (symplectic, pairing ±1) pair must intersect, excluding the dumbbell — is correctly flagged as an immediate consequence *not printed as a named corollary* in the primer. The derivation is sound.

### Tripwire: diamond `3L−1` / `4L²` — confirmed NOT in those texts

- **Newman–Ziff (full text re-fetched this check):** no `2L−1`, no `3L−1`, no `4L²` minimizer count, no enumeration or classification of minimal wrapping configurations of any kind. Only event probabilities and the first-wrapping occupation-number distribution.
- **Farb–Margalit Ch. 1 (re-fetched this check):** contains no percolation or lattice-count content at all; the diamond statements are trivially absent.
- **Web/arXiv searches this check** for minimal wrapping site counts (`2L−1`/`3L−1` onset) and for diamond/tilted-cell torus wrapping polynomials: nothing. The nearest objects remain the axis-aligned square-cell wrapping polynomials of Akhunzhanov–Eserkepov–Tarasevich (arXiv:2204.01517, one direction only) — as #723 says.

**Tripwire answer confirmed: the diamond `3L−1` onset and `4L²` minimizer count do not appear in print, and neither does the axis `2L−1`.** #723's novelty claim for the #710 classification stands.

---

## KEEP / CORRECT / DROP table

Tags: **PRIMARY_TEXT_READ** = full text fetched; **ABSTRACT_ONLY** = listing/abstract only; **[LIT]** = literature context rather than a checked primary claim. "Re-verified" = by this check's own fetch, not merely by reading #723.

| # | Item (from #723) | Tag | Verdict | Notes |
|---|---|---|---|---|
| 1 | Newman–Ziff wrapping defs `R^(h/v/b/e/1)`, relations, "spiral" note, Pinson limits | PRIMARY_TEXT_READ | **KEEP** | Re-verified verbatim this check from cond-mat/0005264. |
| 2 | Newman–Ziff has no minimal-set classification, no homology language | PRIMARY_TEXT_READ | **KEEP** | Confirmed against full text this check. |
| 3 | Title/ID pairing: "Fast Monte Carlo algorithm…" ↔ cond-mat/0005264 | PRIMARY_TEXT_READ | **CORRECT** | 0005264 is the PRL ("Efficient Monte Carlo algorithm…", PRL 85, 4104); the printed title is the PRE companion (PRE 64, 016706; cond-mat/0005397). Content unaffected — quoted defs are in the fetched PRL. |
| 4 | Pruessner–Moloney: winding admissibility/coexistence, "heuristic… refer to the standard literature" | PRIMARY_TEXT_READ | **KEEP** | Consistent with the paper's known content; not independently re-fetched this check (not in the required list). Quotes are internally coherent and correctly used as admissibility-only. |
| 5 | Langlands–Pichet–Pouliot–Saint-Aubin: planar geometries, no torus/winding | PRIMARY_TEXT_READ | **KEEP** | Not re-fetched this check; consistent with the paper's programme as described. |
| 6 | Pinson: ABSTRACT_ONLY, used only as probability source, no config classification claimed | ABSTRACT_ONLY | **KEEP** | Discipline maintained — #723 and this check both leave it abstract-only; bot-challenge on Springer/ADS as described. |
| 7 | Akhunzhanov et al. 2022: one-direction wrapping polynomials, orbit-stabilizer divisibility, no two-direction, no `2L−1` | PRIMARY_TEXT_READ | **KEEP** | Not re-fetched this check; my searches corroborate the one-direction/no-diamond scope. |
| 8 | Farb–Margalit Prop 1.5 + primitive definition | PRIMARY_TEXT_READ | **KEEP** | Re-verified verbatim this check. |
| 9 | Farb–Margalit §1.2.3 intersection formulas + SL(2,Z) invariance | PRIMARY_TEXT_READ | **KEEP** | Re-verified verbatim this check. |
| 10 | Derived consequence: disjoint essential ⇒ homologous; spanning pair must intersect; dumbbell excluded via intersection form (not a printed corollary) | PRIMARY_TEXT_READ | **KEEP** | Correctly flagged as derived, not quoted; derivation checked and sound. |
| 11 | Hatcher: cycle rank = \|E\|−\|V\|+1, cube Example 1.22, no theta/dumbbell graph-classification passage | PRIMARY_TEXT_READ | **KEEP** | Not re-fetched this check; arithmetic checks (12−8+1=5); the "homotopy alone cannot exclude the dumbbell" framing is correct. |
| 12 | Systolic survey Ikonen–Marti–Vikman arXiv:2607.22290, Loewner inequality quote, no discrete L1 version | PRIMARY_TEXT_READ | **KEEP** | Citation verified this check at abstract level (authors/date/topic match); verbatim full-text quotes accepted from #723's fetch. Loewner unpublished → cited via survey: correct practice, [LIT]. |
| 13 | Tripwire: diamond `3L−1` + `4L²` not in print; axis `2L−1` also not in print | [LIT] | **KEEP** | Confirmed by this check's re-fetches and searches (negative result, bounded scope as #723 itself states). |
| 14 | Discrete systolic bound `L(\|a\|+\|b\|)` not found in print; two-line covering-space reconstruction | [LIT] | **KEEP** | Reconstruction is correct (lift to Z², L1 edge count, primitive classes attainable); correctly labeled as not-in-print. |

**DROP: none.** No claim in #723 failed verification; nothing needs removing.

## Corrections to carry into #723 (if amended)

- C1 (row 3): relabel the fetched Newman–Ziff text as the PRL (cond-mat/0005264, "Efficient Monte Carlo algorithm and high-precision results for percolation," PRL 85, 4104); cite the PRE companion (cond-mat/0005397, "Fast Monte Carlo algorithm for site or bond percolation," PRE 64, 016706) separately if the title is kept. Bibliographic only.

## Not established (scope of this check)

- Pruessner–Moloney, Langlands et al., Akhunzhanov et al., Hatcher rows were accepted on #723's quoted primary text plus corroborating searches, not re-fetched in full this check (the issue's required re-fetches — Newman–Ziff, Farb–Margalit — were done in full).
- Pinson remains ABSTRACT_ONLY in this check as well; per the ticket it stays so unless the body is actually fetched.
- The negative tripwire result is bounded by the searched strands and date (2026-09-12), as #723 itself says.
