# #637 retrieval note — critical polynomials and sector crossings: what is already published

**Ticket:** LightChainr/Matching-One#637. **Date:** 2026-09-07. **Mode:** literature
retrieval, theory input only. No compute entered `docs/STATUS.md`; no ticket closed.
Every answer below is cite-or-gap, with the exact status (theorem vs conjecture)
marked as the source itself marks it.

Headline answer up front, because the dispatch pack says Q1 decides whether #636
is funded:

> **Q1 verdict: a weaker sibling of our identity is published (Mertens–Ziff 2016),
> but our identity as stated — matching function = difference of two *sector
> amplitudes* with different connectivity classes on the two sides — is NOT in the
> literature.** Mertens–Ziff give an exact finite-size identity equating a cluster
> count difference to a wrapping-probability difference, both taken at the SAME
> bond class per side (primary NN for primal, matching NN+NNN for white). That is
> #628's structure, not #635's. A targeted search of Scullard–Jacobsen's own
> writing, Pinson, and the sector-decomposition literature finds nothing that
> decomposes a single-lattice matching function into two topological sectors of
> one transfer matrix. The gap is real but narrow, and the Mertens–Ziff identity
> must be cited as the nearest neighbor.

---

## Primary source read at first hand: Jacobsen, arXiv:1507.03027

Read from the arXiv HTML (v1, 33 pages, 7 figures, 5 tables). This is the paper
`references.bib` already carries as `Jacobsen2015EigenvalueIdentity`.

### 1. The exact form of the criterion

- **Algebra:** periodic Temperley–Lieb (pTL), with the augmented
  Jones–Temperley–Lieb extension used to close the cylinder in the horizontal
  direction (extra generator E_{2n-1}).
- **Sectors:** the transfer matrix T, restricted to the s=0 (zero through-strings)
  block, decomposes into two topologically distinct sub-sectors labeled
  **open** vs **closed** (FK clusters vs dual FK clusters that penetrate the
  cylinder). The critical criterion is
  `Λ_open(p) = Λ_closed(p)` at the root of the graph polynomial `P_B(q,v)`,
  quoted by the paper as: "P_B(q,v)=0 ⇔ Λ_open = Λ_closed, valid for a basis B of
  size n×m, with n finite and m→∞."
- **Order of limits — the part #636 must not get wrong:** **m → ∞ FIRST**, so the
  basis becomes a semi-infinite cylinder of circumference n, and the criterion is
  applied at finite n; only afterwards is the outer n → ∞ extrapolation taken.
  This is stated in the abstract ("bases B are semi-infinite cylinders of
  circumference n") and again in the main result.
- Context: this is the eigenvalue reformulation of the Scullard–Jacobsen critical
  polynomial `P_B ≡ R_2 − R_0` (probability of wrapping in both directions minus
  probability of wrapping in neither) defined on finite periodic bases
  (Scullard–Jacobsen, arXiv:1209.1451; Scullard 2012 deletion–contraction version,
  arXiv:1207.3340).

### 2. What is proved versus what is numerical

- The **eigenvalue-equality criterion itself is proved**, by the intermediate
  value theorem: for v≫1 the paper argues `Λ_open > Λ_closed`, for v≪1 the
  reverse, both grow exponentially in m, so a root must exist where the two are
  equal; the paper says "We can prove this statement as follows" and calls it
  "our main result". **This part is a theorem within its setting.**
- What is NOT proved: that the n → ∞ limit of the root p_c(n) equals the true
  p_c. The paper's own language for the exactly-solved cases is "observed to be
  independent of the size" / "was observed to converge quickly" — the n-independence
  on solvable lattices is empirical, supported afterwards by the CFT argument
  (both sectors carry the same magnetic exponent x_m = 5/48, so the free-energy
  difference is o(n^{-2})). **If we consume the extrapolated number as exact, we
  are inheriting an unproved extrapolation step — flag it as such.**
- **Post-publication status of the number itself (new fact for the ledger).**
  Yang & Zhou, J. Phys. A 57, 258001 (2024), doi:10.1088/1751-8121/ad4d2c, computed
  exact cylindrical p_c(n) up to n = 24 and showed the 2015 extrapolated value
  `0.59274605079210(2)` is **incorrect**; their corrected value is
  `0.5927460507896(1)`. Jacobsen's Reply (J. Phys. A 57, 258002, 2024,
  doi:10.1088/1751-8121/ad4d33) published an updated extrapolation giving
  `0.59274605079016(1)` (the Reply's number as recorded in the HandWiki
  percolation-threshold table and OEIS A377420). **Implication for this
  repository:** the "Jacobsen 2015 interval [0.59274605079208, 0.59274605079212]"
  quoted in `docs/astra/Q4-why-square-site-resists.md` and used by
  `scripts/threshold_claim_intake.py` as a published interval is **not the current
  state of the literature**; the Yang–Zhou and Jacobsen-Reply values sit ~2e-12 to
  ~2e-11 BELOW it, outside the 2015 interval. The census results referenced
  against the 2015 interval should be re-referenced against both newer values
  before any claim of non-membership is reused. This is a status-relevant finding
  but per ticket boundaries it is recorded here, not in docs/STATUS.md.
  (In 2022, Mertens' transfer-matrix/cluster-counting route gave
  0.592746050786(3) — consistent with the corrected values, not with the 2015
  digits.)

### 3. Finite-size structure (the correction part, Q4-adjacent)

- The paper's FSS ansatz: `p_c(n) = p_c + Σ_k A_k n^{-Δ_k}` with measured
  Δ_1 = 4.0001(2), Δ_2 = 6.00(1) for square-site, conjecturing Δ_k = 2(k+1), i.e.
  a clean tower 1/n^4, 1/n^6, 1/n^8, .... For kagome bond, A_1 = 0 and the leading
  correction is Δ_2 ≈ 6. The CFT backing: both sectors carry the same critical
  exponent, and the difference `f_open − f_closed = o(n^{-2})`, with the
  pseudo-critical shift `p_c(n) − p_c = O(n^{-4})` (eq. 50 of the paper).
- **So the leading correction exponent is 4 for square-site (n^{-4}), NOT the
  trivial 1/n^2; this is exactly the structure #636 would consume.** The paper
  marks the full tower Δ_k = 2(k+1) as a conjecture ("conjecture Δ_k = 2(k+1)"),
  not a theorem.
- Scullard–Jacobsen arXiv:1209.1451 independently discusses the same convergence
  behavior for the critical polynomial (the "probabilistic definition" paper) and
  is where the R_2−R_0 form and its very fast convergence are documented; Xu et
  al. 2021 (patchy particles on Archimedean lattices) use the FSS form
  `P_B ≃ a_1 (p−p_c) L^{y_t} + b_1 L^{y_1} + b_2 L^{y_2}` and report the critical
  polynomial's finite-size corrections are "much smaller" than for other wrapping
  quantities. None of these works match the correction structure against an
  independently measured finite-size observable of OUR kind (a signed difference
  of wrapping indicators on two lattices, D(C)); they all match against a known
  answer or against MC estimates of the same threshold. **The specific comparison
  in Q4 — criterion root vs independent finite-size observable, no known answer
  assumed — has no published instance found by this search.** That is where a
  contribution is still plausible.

### 4. Square-site specifically (feeds Q3)

- Handled inside pTL via the FK (q=1) formulation with the site-specific R-matrix
  `check R_i = E_{i+2}E_i + v E_{i+1}`, p = v/(1+v). Site percolation is treated
  on the same cylinder geometry, s=0 sector, reduced state space.
- Scale reached: n×∞ bases up to **n_max = 21** in 2015 (Table 2 gives p_c(21) =
  0.592744551481371482002735520463 unextrapolated); extended to **n = 24** by
  Yang–Zhou 2024. Best current published value: **0.59274605079016(1)** (Jacobsen
  Reply 2024) with Yang–Zhou at 0.5927460507896(1) and Mertens 2022 at
  0.592746050786(3) — all three agree pairwise at their own quoted precisions,
  and all three sit outside the 2015 quoted interval.
- The paper's own statement on solvability: "Site percolation on the square
  lattice. This renowned problem is unsolved essentially because the four-regular
  square-lattice hypergraph is not selfdual." Site square is listed among the
  three significant *unsolved* cases (with kagome bond and square SAP). Ziff,
  in his 2021 retrospective (Physica A 564, 125557, arXiv:2101.00550, writing
  about the Stauffer challenge), states it more bluntly: "the square lattice with
  site percolation cannot be put in the self-dual form that is required to use
  those methods. I believe it is insoluble."

---

## Q1 — is "matching function = difference of two sector amplitudes" published?

**Cited sibling, gap for the identity itself.**

What IS published, at first hand from Mertens & Ziff, "Percolation in Finite
Matching Lattices", J. Stat. Mech. (2016) / arXiv:1603.07289 (v2, Nov 2016):

- The **finite-size matching identity**, exact on any L×L torus:
  `N_L(p) − N̂_L(1−p) − L² χ(p) = R^x_L(p) − R̂^x_L(1−p)`,
  where N_L is the black cluster count on the primary lattice, N̂_L the white
  cluster count on the matching lattice, χ(p) the matching polynomial (for square:
  p − 2p² + p⁴), and R^x the probability of wrapping in one or both torus
  directions. This is an exact theorem, generalizing Sykes–Essam, and the paper
  explicitly notes "The criterion that follows is related to the criterion
  Scullard and Jacobsen use to find precise approximate thresholds, and our work
  provides a new perspective on their approach."
- So: the **difference of a primal quantity and a matching-lattice quantity being
  exactly equal to a difference of wrapping probabilities** is published and exact
  at finite size. Note the shape, though: in Mertens–Ziff, each side's indicator
  uses ONE connectivity class per lattice (black on NN primary; white on NN+NNN
  matching), and the difference is between the lattices' separate observables at
  dual weights p and 1−p.

What is NOT published (this search): the decomposition of a **single**
single-lattice matching-function observable into a **difference of two
topological sectors of one transfer matrix** (open vs closed amplitudes of the
same pTL operator, as in Jacobsen's Λ_open − Λ_closed). In particular:

- Scullard–Jacobsen (1209.1451, 2012; JPA 2015; PRR 2, 012050 (2020)) define the
  critical polynomial as R_2 − R_0 (two-direction wrapping minus no-wrapping) —
  that is a difference of two *wrapping-type* probabilities of the same lattice,
  not of two connectivity classes, and not a matching function.
- Pinson (J. Stat. Phys. 2011, arXiv:1006.1307 context) gives exact torus
  wrapping probabilities by homology subgroup (π_G = Z_G/Z), later extended by
  Arguin to Q-Potts FK clusters; these are sector-resolved amplitudes of the SAME
  model, with no matching-lattice side. The Ising/Potts-spin-cluster versions
  (arXiv:1402.0879 and follow-ups) likewise stay within one connectivity.
- No source found that writes a signed combination "primary wrapping −
  matching-lattice wrapping" and identifies it with two transfer-matrix sectors.
  Scullard–Jacobsen's own probabilistic definition of P_B is the closest published
  object to a "sector difference" and it is a two-direction/no-wrap difference on
  ONE lattice.

**Verdict for Q1: the identity is not published. The bond-vs-site clause in the
ticket matters too:** Mertens–Ziff's theorem covers site percolation directly
(their example lattice is exactly our square-site/square-with-diagonals matching
pair, χ(p) = p − 2p² + p⁴). So the published sibling is in fact a SITE result —
what remains unpublished is the SECTOR decomposition, not the bond/site swap. Any
#636 claim of novelty must therefore be phrased narrowly: not "first to connect
matching functions to wrapping differences" (Mertens–Ziff did that, exactly, in
2016) but "first to realize the matching function as a difference of two
topological sectors of one transfer matrix." If #635's map reproduces the
Mertens–Ziff identity plus a sector reading on top of it, the citation obligation
is to both Mertens–Ziff 2016 and Jacobsen 2015.

Related repo context honored: #628's finding (r_b + r_w = 2 fails in
118133/262144 square-bond configurations) is precisely the statement that
`R^x_L(p) + R̂^x_L(1−p) = 1` — the p↔1−p-symmetric reading of the Mertens–Ziff
identity — is NOT an identity at finite size; Mertens–Ziff's exact statement
instead keeps the −L²χ(p) term, consistent with #628. No conflict.

## Q2 — the matching lattice's state space (NN+NNN strip connectivity)

**No published boundary-state theory for the matching strip was found; gap with
strong circumstantial evidence.**

- The matching lattice of the square lattice is the square lattice with both
  diagonals (NN+NNN). As a strip-connectivity problem on w frontier sites, a
  single NNN bond along the frontier (i, i+2 is connected through a diagonal of
  the NEXT row — in a row-by-row build, the diagonals of the row just closed are
  horizontal-distance-2 bonds on the frontier) creates a connectivity graph on
  the frontier that is NOT a noncrossing partition: e.g. sites i and i+2
  connected while i and i+1 are not forces a crossing structure in the
  nesting-nested-arc picture. The standard noncrossing boundary states (Catalan
  C_w for pure NN connectivity, as in the repo's
  `scripts/noncrossing_connectivity_codec.py`) therefore do NOT bound the
  matching side.
- Prior art found, and what it settles:
  - **Pinson's homology/sector work and its descendants** (Langlands–Pinson–
    Saint-Aubin; Arguin 2006; arXiv:0812.2925 etc.) handle wrapping probabilities
    on the torus for percolation, i.e. plane triangulation-era connectivity —
    they do not provide strip boundary-state counts.
  - **Mertens 2022 (arXiv:2109.12102)**, "Exact site-percolation probability on
    the square lattice," is the closest technical neighbor: it builds a
    transfer matrix for exactly our kind of object (spanning-configuration
    generating function on n×m square-site) using "signatures" = connectivity
    states of a row, and Appendix A discusses the count of signatures via
    balanced-parenthesis / Catalan-triangle arguments. **But it is NN-only for
    the black side** (spanning of primary clusters); it does NOT carry NN+NNN
    connectivity for the white/matching side on the same footing. Its signature
    machinery would need nontrivial extension for the matching side, and the
    paper's Catalan-based state counting confirms the NN side but says nothing
    about the NN+NNN side. So the Q2 question — correct boundary-state space for
    a matching (NN+NNN) connectivity on a strip, and whether it is still counted
    by Catalan numbers — is NOT answered by the closest prior art.
  - Transfer matrices for lattice strips with diagonals exist for BOND problems
    on triangulated strips (e.g. Shrock–Tsai-type cyclic/Möbius strip transfer
    matrices, which treat the triangular lattice as a square strip plus
    diagonals — arXiv:cond-mat/0407070 family; Potts transfer matrices with
    Catalan configuration spaces, Alcaraz-type family-tree papers
    S0010465514003464). Those handle diagonal BONDS on a planar strip within
    noncrossing states, because a single diagonal per face keeps the connectivity
    planar. **The matching side is different: it has BOTH diagonals per face
    plus the NN edges, and the union is non-planar as a connectivity relation**
    (the square-with-both-diagonals graph contains K4 minors on every face, and
    crossing connectivity states on the frontier are generated by row additions).
    None of the strip-literature entries handle this union.
- **Circumstantial evidence on cost:** Jacobsen 2015 reached n = 21 on the
  square-SITE problem using pTL reduced states ~4^n (the paper says the
  eigenvalue formulation "allows us to double the size n for which T_c(n) can be
  obtained, using the same computational effort", and reduced-state hashing). If
  the matching side required only the same noncrossing class, the NN-only
  machinery at the same widths would have been cheap; instead, the one
  published exact treatment of the matching pair Mertens–Ziff 2016 works
  enumeratively/analytically at small L and the transfer-matrix treatments
  (Mertens 2022, signatures) explicitly bound themselves to the primary side.
  This is consistent with the matching side being genuinely harder, but it is
  not a proof.
- **Verdict for Q2: gap.** No published boundary-state space for NN+NNN strip
  connectivity; Catalan counting is not expected to survive (the frontier
  connectivity relation admits crossing configurations, e.g. through a K4 on a
  single face), so #636's cost model, if built on noncrossing states, is suspect
  in the direction of being an UNDERestimate. #638's pure enumeration (open
  ticket, correct instrument) should be the arbiter; this retrieval could not
  find literature doing #638's job.

## Q3 — state of the art for square-site, and why the site case resists

**Cite, and the resistance is stated by the sources themselves.**

Best published treatments, chronologically:

1. **Scullard 2012 (arXiv:1207.3340):** generalized critical polynomials by
   deletion–contraction, bases to 36 bonds, all Archimedean lattices except
   kagome; the (4,8²) prediction "though not exact, is not ruled out by
   simulations". Not site-square.
2. **Scullard–Jacobsen 2012/2013 (arXiv:1209.1451):** transfer-matrix
   computation of the critical polynomial, bases up to 96/162/243 edges for
   (4,8²), kagome, (3,12²) bond problems; explicitly notes "the alternative
   definition of P_B(p) can be applied to study site percolation problems" —
   this is the door Jacobsen 2015 then walked through for square-site.
3. **Jacobsen 2015 (arXiv:1507.03027):** eigenvalue method, square-site to
   n = 21, extrapolated 0.59274605079210(2) — later corrected, see §2 above.
4. **Mertens 2022 (arXiv:2109.12102):** transfer matrix over "signatures" for
   the number of spanning configurations, F_{n,m}(z), exact for small systems,
   extrapolated p_c = 0.592746050786(3). **This is the single best independent
   route to the square-site number and its state-counting is via
   Catalan/balanced-parentheses structures.**
5. **Yang–Zhou 2024 (arXiv-equivalent, J. Phys. A 57, 258001):** exact
   cylindrical p_c(n) up to n = 24, c ≈ 2.7459 state growth, value
   0.5927460507896(1); triggered Jacobsen's Reply (258002) re-extrapolating the
   2015 data to 0.59274605079016(1).

Why the site case resists, in the sources' own words:

- Jacobsen 2015: "unsolved essentially because the four-regular square-lattice
  hypergraph is not selfdual."
- Ziff 2021 (arXiv:2101.00550): "the square lattice with site percolation cannot
  be put in the self-dual form that is required to use those methods. I believe
  it is insoluble, but you never know — new things come along."
- The Scullard–Jacobsen line solves what it can by *building the criterion into
  an approximation scheme* rather than by solving; the site case is not special
  there — the method is exactly as applicable, and the approximation converges
  (fast) but has no exact special case to snap to, because square-site is
  neither self-dual nor self-matching (matching lattice ≠ itself).
- `docs/astra/Q4-why-square-site-resists.md` asks whether the degree/height
  regularity of exactly-solved thresholds is mechanism-bound and whether
  square-site provably falls outside the reachable set. **The literature does
  NOT answer that question.** What the literature DOES establish (and this is
  the one clause the Q4 file should absorb rather than leave as an open ask):
  the three negative facts the file lists (not self-dual, not self-matching, no
  star-triangle) are each stated or implicitly assumed in the sources above, and
  no author since Sykes–Essam has proposed a fourth mechanism for the site case.
  The file's open question ("is that obstruction provable?") remains genuinely
  open — no published theorem bounds the reachable set of exact-threshold
  mechanisms. **Retire nothing; the file's Q remains live, but its premise
  (that the negative facts are agreed) is confirmed by Jacobsen 2015, Ziff 2021,
  and the whole Scullard–Jacobsen line.**

## Q4 — corrections/correction structure between our observables and that spectrum

**Gap, and it is the only live direction of real contribution.**

What exists:

- Jacobsen 2015 §8: correction tower for p_c(n), Δ_k = 2(k+1) conjectured,
  leading 1/n^4, CFT-motivated via x_m = 5/48 and c = 0 cylinder scaling
  (eqs. 46–50). This is a correction analysis of the CRITERION against a known
  pc — the "known answer" mode.
- Xu–Chen– Cummings-style FSS of P_B in MC contexts (arXiv:2111.12568 and
  references): `P_B ≃ a_1(p−p_c)L^{y_t} + b_1 L^{y_1} + b_2 L^{y_2}`, i.e. the
  crossing probability itself has an FSS correction structure; again consumed
  against known thresholds.
- Mertens–Ziff 2016: exact finite-size identity (their eq. 4) — the one place
  where a wrapping-difference observable is tied to an exact statement at every
  finite L, not asymptotically.

What does NOT exist (searched, not found anywhere):

- A correction/closing analysis of a sector-crossing criterion matched against
  an **independently measured finite-size observable that is NOT the criterion
  itself and is not validated against a known pc**. Every published use of
  R_2−R_0 or Λ_open=Λ_closed validates by (a) recovering a known exact answer
  on solvable lattices, or (b) reproducing MC thresholds. No one has asked "does
  the root of the sector-crossing criterion, at finite size, track another
  measurable with a predictable error structure, in a regime where no exact
  answer exists?" That is the exact shape of our D(C)-vs-criterion comparison.
- Also absent: any published correction tower for a criterion applied to a
  DIFFERENCE of wrapping indicators on TWO lattices (our D(C) is signed and
  changes sign near p_c; the published towers are all for nonneg wrap
  probabilities of one lattice).

**Verdict: Q4 is a genuine gap.** If #635's map holds, the publishable content
is exactly this correction structure, and the ticket's warning is confirmed:
novelty belongs there, not to the transfer matrix.

## Q5 — symmetry-forced vanishing of first-order response for transfer matrices

**Gap at the level of a named theorem for stochastic/transfer-matrix operators;
the Hamiltonian and master-equation cases are classical.**

- Hamiltonian case: Wigner–Eckart / selection-rule machinery (any standard
  reference; matrix elements vanish unless Γ_f* ⊗ Γ_T ⊗ Γ_i contains the
  identity) — that is the Wigner–Eckart route, not a transfer-matrix result.
- Master equations: Hänggi (1978) type results — established. (Not re-verified
  here per ticket boundary; #601 settled the named-rule landscape for
  equivariant Markov generators, i.e. the gap #601 found is real and this
  search found nothing to fill it.)
- For **transfer matrices / Perron–Frobenius operators with a reflection
  symmetry**: the search surfaces only the generic Perron–Frobenius structure
  (nonnegative matrices, unique dominant eigenvalue, positivity of the leading
  eigenvector — e.g. the PF-theorem exposition in Phys. Rev. E 105, 044105
  (2022)) and FDT-type relations for equilibrium Markov/transfer operators
  (Kubo–Ruelle theory; detailed-balance reversibility). The relevant folklore
  — that for an equilibrated transfer operator T with a reflection R
  commuting with T, and a perturbation/response channel of odd parity under R,
  the first-order response ∂Λ/∂F vanishes because the leading eigenvector is
  even and the adjoint channel is odd — exists in the literature ONLY in
  application-specific forms: parity-sector decompositions of transfer matrices
  (Ising/loop models) are standard, and "the dominant eigenvalue is even-parity
  so odd-sector observables have zero linear response in equilibrium" is used
  implicitly in equilibrium-statistical-mechanics computations (it is the
  reason equilibrium susceptibilities to odd fields are computed from
  eigenvalue differences, not from ∂Λ_0). But no source was found that states
  and proves it as a theorem in the general transfer-matrix setting with a
  name attached. **#610's transfer-matrix selection lemma is therefore not
  anticipated by a named published theorem as far as this search can tell;
  the gap #601 identified stands.** Caveat honestly stated: a proof like this
  is elementary enough that it may be folklore embedded in homework/exercise
  culture or in specialized papers this search's queries did not reach; the
  recommendation stands nonetheless — keep #610's lemma, cite it as "we could
  not find this stated elsewhere."

---

## Per-ticket deliverable checklist

- Q1: cite-or-gap → **gap for the identity, cite Mertens–Ziff 2016 as the
  published sibling (site-exact, torus-exact) and Jacobsen 2015 for the sector
  machinery.** The claim for #636, if it proceeds, must be phrased as novelty in
  the sector decomposition, not in connecting matching to wrapping.
- Theorem status of Jacobsen's criterion: **proved (IVT), within the m→∞-first
  setting; the n→∞ step and the correction tower are conjectural/empirical.**
- New ledger-relevant fact: **the 2015 published interval for square-site pc has
  been superseded** (Yang–Zhou 2024; Jacobsen Reply 2024; consistent with
  Mertens 2022). Any census/exclusion conclusion currently keyed to
  [0.59274605079208, 0.59274605079212] needs re-referencing before reuse.
  Flagged here only; docs/STATUS.md untouched per boundaries.
- Q2: gap, with the cost-model implication spelled out (noncrossing assumption
  likely wrong on the matching side; #638 is the right instrument).
- Q3: cited (five sources above); `docs/astra/Q4-why-square-site-resists.md`
  premise confirmed, its open question remains open; retire nothing.
- Q4: gap — and the only gap where a real contribution is plausible.
- Q5: gap at the named-theorem level; #610's lemma stands.

## What I could not do

- I could not access the full text of Yang–Zhou 2024 / Jacobsen Reply 2024
  (paywalled at IOP; bot-managed at HAL). The Comment's abstract, the Reply's
  published corrected value 0.59274605079016(1) as recorded in HandWiki's
  percolation-threshold table and OEIS A377420, and the two DOIs are verified;
  the Reply's internal argument is not, and the note above quotes the secondary
  sources for the Reply's number.
- I did not fetch Pinson 2011 at first hand; the sector structure is quoted
  from the Langlands/Arguin descendants (arXiv:0812.2925, arXiv:1402.0879) and
  from secondary citations of Pinson. If the sector claim in Q1 is load-bearing
  for a decision, Pinson's paper itself should be pulled before #636 is funded.
- This was a retrieval ticket; no enumeration or computation was attempted, so
  no Q2 count can be produced here. #638 remains the right instrument for the
  state-space question.
