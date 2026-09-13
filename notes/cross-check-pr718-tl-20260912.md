# Round B cross-check of #718: P398 as periodic identified-connectivity TL

2026-09-12. Peer check of the Astra package #718 (TL half) per #721. The note
and script under review are the author's, not mine; this file records an
independent verification, not a new claim. Scope: read #718, hand-check Phi at
small widths, confirm the two literature anchors, and audit the stated
non-claims. Parent: #650.

## 0. Source marks

| Source | Access | Mark |
| --- | --- | --- |
| #718 `notes/p398-is-periodic-tl-20260912.md`, `scripts/p398_tl_fattening.py` | full read from branch `pr718`; script executed locally at w=2..8 | PRIMARY_TEXT_READ |
| Pearce, Rittenberg, de Gier, Nienhuis, arXiv:math-ph/0209017v2 (J. Phys. A 35 L661) | full HTML text fetched, section 2 (Eqs. (2.7)-(2.18)) | PRIMARY_TEXT_READ |
| Cantini, Sportiello, arXiv:1003.3376v1 (JCTA 118 (2011) 1549) | HTML text fetched through section 4.1; section 4.5 gyration proof and later sections beyond the fetched text | PRIMARY_TEXT_READ (partial; see section 6) |
| Levy/Martin/Saleur CTL algebra, Kreweras complement, #708/#709/#715 internals | cited here only as reported by #718 or Pearce et al.; not independently fetched this round | [LIT] |

Verdict summary: every checked claim of #718 is **KEEP**. No CORRECT, no DROP.
Residual caveats are listed in section 6 and none affects the identification.

## 1. Independent hand-check of Phi at w=2 and w=3

Definitions as in #718: L_i=2i, R_i=2i+1 (zero-based); a block
B={i_1,...,i_k} in cyclic order contributes pairs (R_(i_j), L_(i_(j+1)));
a singleton contributes the adjacent pair (2i, 2i+1).

**w=2.** C_2=2; the two noncrossing partitions give the two noncrossing
matchings on 4 endpoints, both distinct:

    Phi({0}{1})   = (0,1)(2,3)
    Phi({0,1})    = (0,3)(1,2)

Both generator identities, both families, including the no-op cases:

    Phi(detach_0 {0,1}) = Phi({0}{1}) = (0,1)(2,3);
        e_0 Phi({0,1}): (0,3)(1,2) -> remove (0,3),(1,2), add (0,1),(2,3). OK.
    Phi(join_0 {0}{1})  = Phi({0,1})  = (0,3)(1,2);
        e_1 Phi({0}{1}): (0,1)(2,3) -> remove (0,1),(2,3), add (1,2),(0,3). OK.
    detach_0 {0}{1} = {0}{1};  e_0 on (0,1)(2,3): 0,1 already paired -> no-op. OK.
    join_0  {0,1}  = {0,1};   e_1 on (0,3)(1,2): 1,2 already paired -> no-op. OK.

**w=3.** C_3=5 = Bell(3); all five partitions of {0,1,2} are noncrossing, and
the five matchings are distinct and noncrossing, so Phi is a bijection onto
the Catalan set (not merely a cardinality match):

    {0}{1}{2} -> (0,1)(2,3)(4,5)
    {0,1}{2}  -> (0,3)(1,2)(4,5)
    {0,2}{1}  -> (0,5)(1,4)(2,3)
    {0}{1,2}  -> (0,1)(1,5)(3,4)
    {0,1,2}   -> (0,5)(1,2)(3,4)

Identities checked, including the cyclic seam (i=w-1):

    Phi(detach_1 {0,1,2}) = Phi({0,2}{1}) = (0,5)(1,4)(2,3);
        e_2 Phi({0,1,2}): (0,5)(1,2)(3,4) -> remove (1,2),(3,4),
        add (2,3),(1,4) = (0,5)(1,4)(2,3). OK.
    Phi(join_2 {0}{1}{2}) = Phi({0,2}{1});
        e_5 Phi({0}{1}{2}): (0,1)(2,3)(4,5) -> seam, remove (4,5),(0,1),
        add (5,0),(4,1) = (0,5)(1,4)(2,3). OK.

**K^2 = one-site rotation**, hand-checked at w=3 on two orbits:

    K({0}{1}{2}) = {0,1,2}  and  K({0,1,2}) = {0}{1}{2}, so K^2 = rotation (trivial here);
    K({0,1}{2}) = {0,2}{1}  and  K({0,2}{1}) = {0}{1,2} = rotation of {0,1}{2}. OK.

K is not an involution at w=3 (orbit length 3 on the singleton block states is
consistent with K^(2w)=id), matching #718's warning.

**Inverse direction.** In a noncrossing perfect matching of 2w cyclic
endpoints every arc spans an even number of endpoints, so pairs are
(even, odd); mate(2i+1)/2 is a well-defined permutation whose cycles are the
blocks, and noncrossing forces each block convex in cyclic order. This
reconstructs (1) of #718, so the bijection is genuine at every width.

## 2. Pearce et al. math-ph/0209017: IC vs DC

PRIMARY_TEXT_READ. The paper distinguishes exactly what #718 says it does:

- Eq. (2.15) is the periodic cylinder Hamiltonian H = sum_(i=1)^L (1 - e_i)
  in the cylindrical TL (CTL) algebra.
- Eq. (2.16), "periodic (DC)": the quotient keeping at most one
  non-contractible loop, where front/back half-loops are *distinct*;
  dimension (1 + L/2) C_(L/2).
- Eq. (2.17) and the surrounding text, "periodic (IC)": one further quotient
  closing the cylinder into a disk, front/back half-loops isotopic;
  dimension C_(L/2).

With Pearce's L = 2w (2w endpoints), the IC dimension is C_w. P398's state
space of w-point noncrossing partitions has exactly C_w states, so it matches
the IC quotient and cannot match DC ((1+w)C_w for w >= 2). #718 claims IC and
explicitly disclaims DC ("not the cylinder's distinct-connectivity (DC)
representation") -- **confirmed**.

The join-only remark #718 quotes is present ("the terms in the Hamiltonian
may connect disconnected lines but it is not possible to have the reverse
process"). #718's reading that this concerns the line/defect filtration of
the faithful/DC representations and does not forbid detaches on the IC link
patterns is consistent with Cantini-Sportiello's e_j (section 3 below): on
link patterns, e_j, viewed back on the w-point partitions, does detach a
point whenever the arc (j, pi(j)) is not the adjacent pair. I verified this
concretely at w=3 in section 1 (e_2 Phi({0,1,2}) detaches point 1). The
factor-of-two warning ("there are 2w local TL maps, not w") is also correct:
Pearce/CS index e_1..e_L resp. e_1..e_2n on the 2w endpoints.

## 3. Cantini-Sportiello 1003.3376: stationary law imported, not reproved

PRIMARY_TEXT_READ (partial; caveat in section 6). Confirmed against the text:

- Section 2.2, Eq. (4) is exactly the reconnection map #718 uses: e_j acts as
  the identity if (j, j+1) is a pair, otherwise removes (j, pi(j)),
  (j+1, pi(j+1)) and adds (j, j+1), (pi(j), pi(j+1)). The paper states this
  in words, verbatim, and derives the affine TL relations from it.
- Section 2.4, Eqs. (22)-(24): H_n = sum_(k=1)^(2n) e_k, and the
  Razumov-Stroganov statement H_n |s_n> = 2n |s_n> where |s_n> sums square
  FPL configurations refined by boundary link pattern. Stated as Conjecture
  2.1 there; the abstract and the proof skeleton (Lemma 3.1 through Eq. (67),
  "thus completing the proof") show the paper *proves* it by Wieland
  gyration. Square grid, alternating boundary conditions, in bijection with
  ASMs; q = e^(2i pi/3) so -q - q^(-1) = 1, i.e. loop weight 1, so the no-op
  case carries no scalar -- exactly the "loop weight 1 is essential" point of
  #718 section 2.
- #718 section 5 says: "This is an import of their theorem, not a new proof
  or discovery" and explicitly disclaims a Markov dynamics on FPL
  configurations and any eta != 0 extension. **Confirmed**: nothing in #718
  reproduces or extends the CS proof; the stationary law pi_0(pi) =
  FPL_w(Phi(pi))/A_w with the ASM product is a correctly-scoped import.

## 4. Script reproducibility (executed check)

`git show pr718:scripts/p398_tl_fattening.py` executed locally:

- w=2..8: 31,040 join/detach conjugacy equalities -- matches #718 section 6's
  "all 31,040" exactly; bijection, cyclic seam, TL relations, half-step
  square, and readout identities all assert-clean.
- State counts 2, 5, 14, 42, 132, 429, 1430 = Catalan(w), and the w<=6 runs
  agree with the independent restricted-growth-string enumeration.
- Stationary primitive weight sums 2, 7, 42, 429 at w=2..5 equal the ASM
  product values, and the eta = +/-1/4 conjugacy assertions pass.

## 5. Non-claims audit

- **Square-site percolation**: #718 section 3 states "The P398 process
  remains different from microscopic square-site percolation, and its
  parameter eta is not thereby identified with occupation probability p."
  No identification is made or implied anywhere in the note. **Confirmed.**
- **#708 annular rank / #709 certificates**: #718 section 3 states the IC
  disk quotient "must NOT be substituted for #708's homology-preserving
  lifted torus closure"; section 6 keeps "the previously certified
  double-pulse claims ... scoped as before; we do not re-score or extend
  their rank certificates here." The opening even states #709's certificates
  "remain valid". **Confirmed** -- no replacement of the annular-rank work is
  claimed or performed.

## 6. KEEP / CORRECT / DROP

| #718 claim | Verdict | Basis |
| --- | --- | --- |
| Phi(detach_i pi) = e_(2i) Phi(pi), Phi(join_(i,i+1) pi) = e_(2i+1) Phi(pi), all w, incl. no-op and seam cases | **KEEP** | hand-checked w=2,3 (section 1); CS Eq. (4) is the same map with the same no-op convention; script's 31,040 equalities at w=2..8 |
| relevant representation is periodic IC, dimension Catalan(w), not DC | **KEEP** | Pearce Eqs. (2.16)-(2.17): DC dim (1+L/2)C_(L/2) vs IC dim C_(L/2); C_w matches IC only |
| detach half exists on IC link patterns despite the source's join-only remark | **KEEP** | join-only remark concerns the line/defect filtration; CS e_j concretely detaches in partition language; verified at w=3 |
| K = Phi^-1 rho Phi satisfies K^2 = one-site cyclic rotation, K^(2w)=id, K not an involution | **KEEP** | hand-checked at w=3 on two orbits; script asserts half_step^2 = rotation for every state, w=2..8 |
| A_w = ASM numbers 2, 7, 42, 429 at w=2..5; pi_0 = FPL_w(Phi(pi))/A_w | **KEEP** | standard ASM product; stationary solves reproduce the primitive sums; CS theorem (proved) H_n\|s_n> = 2n\|s_n> grounds the import |
| RS stationary description is imported, not reproved; eta=0 only; no FPL dynamics constructed | **KEEP** | #718 section 5 disclaims all three explicitly; nothing in #718 is a proof |
| no square-site percolation identification; eta != p | **KEEP** | #718 section 3 explicit |
| no replacement of #708 annular rank or #709 certificates | **KEEP** | #718 sections 3 and 6 explicit |

Nothing to CORRECT; nothing to DROP.

## 7. Residual caveats (observations, not corrections)

1. My Cantini-Sportiello read is truncated after section 4.1; the gyration
   proof (Prop. 3.4 / section 4.5) was not read in this pass. The theorem
   status (proved) rests on the abstract plus the visible proof completion at
   Eq. (67). This does not touch the import, which cites their result.
2. Notation mapping Pearce L <-> #718 w (L = 2w) is implicit in #718; the
   dimension check C_(L/2) = C_w pins it down, so no correction is needed.
3. The eta != 0 staggered chain and the Kreweras-type naming of K remain, by
   #718's own hedging ("Kreweras-type"), analogies; nothing downstream should
   cite them as published theorems.
4. Pearce et al. is a letter: it *names* IC/DC and dimensions but does not
   prove representation-theoretic fine structure; #718 does not claim
   otherwise.
