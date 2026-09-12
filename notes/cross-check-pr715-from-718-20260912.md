# Round B peer check: #715 against #718 Φ / periodic-IC TL

**Date:** 2026-09-12. **Ticket:** #719 (parent #650). Peer check of
`notes/lit-double-pulse-observability-20260912.md` (PR #715, which I did **not** write).
**Type:** retrieval/verification only. No Huawei. No new production, no enumeration campaign,
no transfer matrix, no hardware. No `docs/STATUS.md` edit; no merge; no close of #713/#715/#718.
**Do not re-derive #709 ranks** — they are out of scope here and are not re-scored below.

The trigger: draft **#718** (`notes/p398-is-periodic-tl-20260912.md`) constructs
Φ: NC-partitions(w) ↔ NC-matchings(2w) with `Φ D_i = e_{2i} Φ`, `Φ J_i = e_{2i+1} Φ`.
If correct, that falsifies #715's central Q3 claim ("TL is join-only on pairings, P398 is
join+detach on partitions, therefore different generators"). This note re-fetches the required
primary sources, independently re-verifies Φ, and issues a KEEP/CORRECT/DROP verdict on every
load-bearing sentence of #715.

## 1. What was fetched this round (all PRIMARY_TEXT_READ unless marked)

| Source | Scope read | How |
|---|---|---|
| #715 note `notes/lit-double-pulse-observability-20260912.md` (pr715 head 308f66c7) | full | `git show pr715:…` |
| #718 note `notes/p398-is-periodic-tl-20260912.md` (pr718 branch, local) | full | `git show pr718:…` |
| Pearce–Rittenberg–de Gier–Nienhuis, *Temperley-Lieb Stochastic Processes*, arXiv:math-ph/0209017v2, J. Phys. A 35 (2002) L661 | §2 Eqs. (2.1)–(2.18), §3 Eq. (3.8) | arXiv HTML |
| Cantini–Sportiello, arXiv:1003.3376v1, JCTA 118 (2011) 1549–1574 | §2.2 Eq. (4) and rules (5a)–(5d); §2.4 Eqs. (22)–(24); §3 proof statement | arXiv HTML |
| Petreczky, *Realization theory … Part II*, ESAIM: COCV **17** (2011) 446–471 | §2.1; Def. 2.8; Thm 2.3(iii) p. 452; Thm 2.7; Rem. 2.5 | Numdam PDF, `pdftotext -layout` |

Not re-fetched this round (retained on #715's own marks): Lucarini arXiv:2502.07908v2,
Petreczky–Wisniewski–Leth arXiv:1605.04414, Bertoin arXiv:0704.3122, Arbib–Manes 1980
(ABSTRACT_ONLY), Isidori 1973 (ABSTRACT_ONLY), Fliess 1974 [LIT], Ho–Kalman 1966 [LIT],
Mueller et al. 2020 (ABSTRACT_ONLY). Rows below that rest only on these are marked "not re-fetched".

## 2. Independent verification of #718's Φ (my own script, not #718's)

I implemented the fattening, the Cantini–Sportiello Eq. (4) reconnection map, and P398's
`detach_i` (point i becomes a singleton) / `join_i` (merge blocks of i and i+1) from scratch and
checked exhaustively for w = 2,…,7:

- Φ maps NC-partitions(w) bijectively onto NC-matchings(2w) (Catalan counts
  2, 5, 14, 42, 132, 429);
- `Φ(detach_i π) = e_{2i} Φ(π)` and `Φ(join_i π) = e_{2i+1} Φ(π)` for every state and every i,
  including the cyclic seam.

**The identities hold.** #718's correction is real: the periodic TL chain has **2w** local
generators, not w, and the even-indexed ones induce the detach half on partitions. #715's Q3
separation argument compared w terms against w terms and missed the factor of two.
(Independent negative finding about my own first attempt: `detach_i` must make i a *singleton*;
deleting the point i outright is not a map on a fixed-w state space and is not what Φ
conjugates. The check above uses the fixed-w semantics, consistent with #718 Eq. (3).)

## 3. Verdict table on #715's load-bearing sentences

Legend: **KEEP** = stands as written. **CORRECT** = quotation/attribution right, conclusion
needs amendment. **DROP** = falsified, must be withdrawn. "Quote" is verbatim from the cited
primary page; where the supporting source was not re-fetched, the quote is #715's own verbatim
quotation and the row says so.

| # | #715 sentence (abridged) | Verdict | Quote (source, location) |
|---|---|---|---|
| T1 | "The only bilinear-Hankel text I read in full is Petreczky, ESAIM: COCV 17 (2011) 446–471, Part II"; Fliess/Isidori/Ho–Kalman not upgraded to PRIMARY | **KEEP** | "The dimension of Σ_min equals the rank of the Hankel-matrix H_Φ of Φ, i.e. dim Σ_min = rank H_Φ." — Petreczky Thm 2.3(iii), p. 452 (Numdam PDF; re-verified this round, verbatim) |
| T2 | "Nothing retrieved matches the specific scalar observer K̂(z)=2(z²+11z+27)/[(z²+10z+23)(z²+11z+26)], the 4×4 determinant −16, or the width-4/5 odd-sector rank attainment … Absence of a retrieved match is **not** an originality certificate" | **KEEP** | Retrieval-gap statement about #709-specific objects; nothing in #718 or in the new fetches matches them, and #718 explicitly does not re-score the ranks. (No single-quote row possible for a negative claim; scoped correctly in #715.) |
| Q1a | "P398's K(τ)=S H exp(τG) H F is an instance of this language, not a new response formalism" (two-time second-order/Volterra kernels for Markov chains) | **KEEP** (not re-fetched) | #715 quotes Lucarini Eq. (18): "G^(2)_{m,Ψ}(k,p) = Θ(k)Θ(p) ⟨m^T (M^T)^p m^T (M^T)^k Ψ, ν_inv⟩", §IV.1.3: "The nonlinear Green functions can be formally seen as **Volterra kernels**". Lucarini arXiv:2502.07908v2, PRIMARY_TEXT_READ in #715; owner's round-1 comment on #715 independently confirms the paper's existence and type. Nothing in #718 touches this. |
| Q1b | "A primary source that uses a two-pulse kernel to recover a *parity-odd* hidden sector after an even quotient … was **not found**. Treat that as a gap in my retrieval, not as a no-go." | **KEEP** | Negative claim, correctly self-scoped in #715; unchanged by #718. |
| Q2a | Word coefficients (2.6), Hankel matrix, finite-rank ⇔ realizable, Thm 2.3(iii) "dim Σ_min = rank H_Φ", Thm 2.7, Rem. 2.5 attributed to Petreczky Part II | **KEEP** | "It turns out that the generating series c has a representation of the form (2.6) if and only if the column rank of H_f is finite." (§2.1, Numdam PDF); "the rank of H_Φ … is the dimension of the linear space spanned by the columns of H_Φ" (Def. 2.8); "constructing a bilinear switched system realization of Φ from the columns of the Hankel-matrix H_Φ" (Rem. 2.5, p. 457). All re-verified verbatim this round. |
| Q2b | Petreczky–Wisniewski–Leth 2016 "does not state a Hankel-rank criterion … not the place to attribute the rank theorem to; the 2011 Part II is" | **KEEP** (not re-fetched) | #715's in-body negative reading (its §2a quotes the 2016 paper's span-reachable/observable minimality test instead). Consistent with T1; the 2011 Part II carries the rank theorem, as re-verified above. |
| Q2c | "A finite 4×4 or 16×16 nonzero Hankel minor is an **instance** of Theorem 2.3(iii)/Theorem 2.7 … Fliess-series/word-coefficient representation, Hankel construction, min order = Hankel rank, partial realization are all prior art and must not be presented as new" | **KEEP** | Direct consequence of Thm 2.3(iii)/Thm 2.7 (p. 452, 456–457), re-verified. The required "confirm #715 kept Petreczky Thm 2.3(iii) correctly" check: **confirmed, verbatim.** |
| Q3a-1 | PRdGn abstract generator Eq. (2.1), line generator Eq. (2.6), cylindrical Eq. (2.15); state space link patterns, dims (2.9)/(2.13) | **KEEP** | "H = Σ_a c_a(1 − w_a) c_a ≥ 0" (2.1); "H = Σ_{j=1}^{L−1}(1 − e_j)" (2.6); "H = Σ_{i=1}^{L}(1 − e_i)" (2.15) — arXiv:math-ph/0209017v2 §2, all re-verified verbatim. |
| Q3a-2 | "Crucially, verbatim: 'the terms in the Hamiltonian may connect disconnected lines but it is not possible to have the reverse process'. So this generator is a **join-only** (irreversible) dynamics on planar connectivities." | **CORRECT** | The quotation is verbatim (PRdGn §2, discussion after Eq. (2.6), in the context of Jordan cells and the right-ideal/word filtration). But it is a statement about the algebra's line/defect filtration in the strip representation, not a representation-independent invariant. In the periodic **IC** quotient on disk link patterns each e_a is the CS Eq. (4) reconnection map, and under Φ the even generators induce the detach half on partitions. The word "So this generator is…" overreaches. |
| Q3a-3 | "TL acts on link patterns (noncrossing **matchings**), P398 on noncrossing **partitions** of w cyclic points; TL's Σ(1−e_j) is **join-only**, while P398's G = Σ_j(J_j + D_j) contains the detach half D_j as well. The two are therefore different generators, and the TL result does not transfer as a statement about G." | **DROP** | Falsified by #718's Φ (`Φ D_i = e_{2i} Φ`, `Φ J_i = e_{2i+1} Φ`), independently re-verified here at w=2..7 (§2 above). On partitions the TL maps split into two interleaved families: odd endpoints join, even endpoints detach. There are 2w local TL maps, not w (CS §2.2: "Define the 2n maps {e_j}_{1≤j≤2n} acting over LP(n)"). |
| Q3a-4 | "The *cylindrical* TL generator Σ_{i=1}^{L}(1−e_i) is structurally the nearest published object to a sum of local generators on cyclic planar connectivities, which is why it is worth naming. But it is not P398's G." | **CORRECT** | It is not merely the nearest object: under Φ and the periodic **IC** quotient (PRdGn Eq. (2.16)→(2.17): "The cylinder can be considered as closed at the top and becomes a disk … We call this case periodic boundaries with identified connectivities (IC) … its dimension is reduced to C_{L/2}") it **is** P398's generator at η=0, with L=2w, on the function-side transpose convention. |
| Q3b | Bertoin 2007 EFC "does have both join and detach, and can be reversible. But it is **exchangeable** (non-local, all partitions …), so it is not the planar/local object P398 uses." | **KEEP** (not re-fetched) | #715's §3b quotation of Bertoin §2.2. Untouched by #718; EFC remains a non-match. The "not an EFC process" half of Q3c survives for exactly this reason. |
| Q3c-1 | "**No named process matching P398's generator was found.**" | **DROP** | The named process exists: the periodic **identified-connectivity**, zero-defect, loop-weight-one TL stochastic process — PRdGn Eq. (2.15) with the IC quotient (2.17), dimension C_{L/2} = Catalan(w); reconnection map exactly CS §2.2 Eq. (4). This is an in-retrieval reading miss, not a retrieval gap: #715 had already retrieved PRdGn as PRIMARY_TEXT_READ but read Eq. (2.15) through the join-only sentence and missed the IC discussion and the factor of two. |
| Q3c-2 | "**Do not call `G` 'the Temperley–Lieb process' or 'an EFC process.'**" | **CORRECT** | Withdraw the TL half: at η=0, G is exactly (the transpose of) the periodic-IC TL intensity generator Σ_a(E_a − I). Keep the EFC half (exchangeable/nonlocal ≠ planar local). For η≠0 the alternating-rate TL chain remains unnamed in what was retrieved; #718 claims no name and no integrability for that family. |
| Q3c-3 | "an unretrieved name is not evidence that the process is unnamed, and it is not an originality claim" | **KEEP** | Methodological caveat; vindicated by the outcome — the name was inside an already-retrieved source. |
| Q4a | "K(0) = 0 with nonzero Hankel rank is the ordinary **strictly proper / relative degree ≥ 1** case … 'zero at zero delay, nonzero rank afterwards' is **generic**, not a named counterexample" | **KEEP** | Consistent with the rank criterion of Thm 2.3(iii) (p. 452), which involves no condition on the first block C B. Re-verified source this round. |
| Q4b | "I found **no published example** of exactly the P398 shape — a two-pulse Markov-chain delay kernel with K(0) = 0 and Hankel rank 4 on an odd sector" | **KEEP** | Scoped negative retrieval on #709-specific objects; #718 leaves the #709 certificates standing ("The previously certified double-pulse claims remain scoped as before"). |
| Q4c | "The *linear* shape is standard and prior; the *P398-specific* two-pulse certificate is not in print." | **KEEP** | As Q4a/Q4b. |

Sources table in #715 (rows 1–11): marks and retrieval URLs are accurate as far as this round
could re-verify (rows 3 and 5 re-verified directly; others not re-fetched, retained). Row 5
(P RdGn) needs an annotation: the §2 IC/DC quotient discussion Eqs. (2.16)–(2.18) and the §3
periodic-IC FPL paragraph (Eq. (3.8)) were retrieved this round and are load-bearing for the
correction; #715 read only the join-only parts of that paper.

## 4. Which of the four pillars survive (explicit, as required)

| Pillar | Verdict after #718 |
|---|---|
| **Hankel** (Petreczky Part II; finite instance, not a new theorem) | **Survives intact.** Thm 2.3(iii) re-verified verbatim; #718 does not touch realization theory. |
| **Volterra** (Lucarini two-time kernel framing of K(τ)) | **Survives.** The TL identification adds an exact model dictionary but does not change the response-theoretic status of K(τ). |
| **EFC** (Bertoin non-match) | **Survives as a non-match** (exchangeable/nonlocal ≠ planar local). |
| **TL-process-id** (#715: "no named process", "do not call it TL") | **Does not survive.** P398's G at η=0 *is* the periodic identified-connectivity zero-defect O(1) TL stochastic process; H = J_0 − J_{w−2} maps to E_1 − E_{2w−3} inside the TL operator algebra (#718 Eq. (5), consistent with Φ J_i = e_{2i+1} Φ re-verified here). |

## 5. Consequences and boundary conditions carried over from #718

- #715's realization-theory citations, its Volterra framing, and the model/observer-specific
  character of #709 remain useful and stand; the process-identification paragraph must be
  withdrawn and replaced by the #718 dictionary.
- The identification is the **IC disk** quotient (front/back identified, dimension C_w). It must
  NOT be substituted for #708's homology-preserving lifted torus closure; P398 remains distinct
  from microscopic square-site percolation, and η is not occupation probability p.
- The CS/RS stationary law (proved, not conjectured: CS §2.4 Eqs. (22)–(24), "H_n |s_n⟩ = 2n|s_n⟩",
  with §3 titled "Proof of the conjecture"; PRdGn Eq. (3.8) "A(n) = Π (3k+1)!/(n+k)! = 1, 2, 7, 42, …"
  for the periodic-IC normalization) becomes applicable to P398 at η=0 as an import. This check
  confirms the CS generator definition matches the TL side of the dictionary (2n maps, no-op when
  (j, j+1) already paired, deleted-loop multiplier −q−q⁻¹ = 1 at q = e^{2iπ/3}).
- No continuum-field claim for #275 and no Jordan-visibility claim for any P398 channel follows
  from the name or from the stationary correspondence.

## 6. What this note does not do

- No novelty claim in either direction; no re-derivation or re-scoring of #709 ranks.
- No merge, no close, no `docs/STATUS.md` edit; #715 and #718 stay open as drafts.
- No re-run of #718's own `scripts/p398_tl_fattening.py`; the check in §2 is an independent
  implementation, which is the point of a Round B peer check.

### Primary sources read this round (for the record)

- Pearce–Rittenberg–de Gier–Nienhuis, arXiv:math-ph/0209017v2 (HTML): §2 Eqs. (2.1)–(2.18),
  join-only sentence, IC/DC quotients, §3 Eq. (3.8). **PRIMARY_TEXT_READ.**
- Cantini–Sportiello, arXiv:1003.3376v1 (HTML): §2.2 Eq. (4) + rules (5a)–(5d) + loop-factor
  remark "−q−q⁻¹ = 1 for q = e^{2iπ/3}"; §2.4 Eqs. (22)–(24); §3 proof opening. **PRIMARY_TEXT_READ.**
- Petreczky, ESAIM: COCV 17 (2011) 446–471 (Numdam PDF, pdftotext): §2.1, Def. 2.8,
  Thm 2.3(iii) p. 452, Thm 2.7, Rem. 2.5 p. 457. **PRIMARY_TEXT_READ.**
