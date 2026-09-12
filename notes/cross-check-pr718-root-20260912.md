# Round B cross-check of #718's aspect-uniform root note (root-vs-law half)

2026-09-12. Ticket #722 (parent #650). ASSIGNED_MACHINE: local. Model:
glm-5.3-flash (deepseek-v4.1-flash returned 429). Peer check of
`notes/aspect-uniform-balance-consistency-20260912.md` on branch `pr718`
(commit 72d2ee6). Peer check, not original retrieval: branches were read with
`git show prNNN:path`; nothing merged; STATUS untouched; #613 is not weakened
by anything below.

## 0. Scope

The root-vs-law half only: #718 §1–§7 (theorem (1), H profile, thin-sequence
contrast, Fisher identity, executed checks). The P398 periodic identified-TL
half of the same PR is out of scope for this round.

## 1. The five required checks

### 1.1 The #718 note itself — PRIMARY_TEXT_READ (branch pr718, full text)

Claim under review: on axis rectangles m>=w>=2, the unique zero p_(w,m) of
M=P_2-P_0 satisfies sup_(m>=w) |p_(w,m)-p_c^site(NN)| -> 0 as w->infinity,
with NO m/w bound and NO w/log(wm) condition; H=P_2/(P_0+P_2) is a conditional
rank profile, not the birth-time mixture F; #613's full-law theorem stays.

### 1.2 #716's two-birth note — PRIMARY_TEXT_READ (branch pr716, full text)

#716 proves law(T) => (delta_0+delta_1)/2 for fixed w as m->infinity (its §2)
and extends it to super-thin sequences w_m = o(log m) (its §6). #718 §5 uses
w_j=j, m_j=ceil(exp(j^2)); w_j/log m_j = 1/j -> 0, inside #716 §6's hypothesis,
so the citation is legitimate. The three-way statement on that sequence —
F(p)->1/2 while p_(w_j,m_j)->p_c and H(p)->1{p>p_c} — is exactly the root-vs-law
distinction both notes claim, and the two notes' limits are orthogonal
(#716's "mean -> 1/2, not q_w" remark concerns the fixed-w cylinder limit).
KEEP.

### 1.3 Duminil-Copin–Tassion site adaptation — PRIMARY_TEXT_READ [LIT]

Read arXiv:1502.03050v3 (HTML): Theorem 1.1 and §1.2. Theorem 1.1(3) (verbatim):
"If (J_{x,y}) is finite-range, then for any β<β_c, there exists c=c(β)>0 such
that P_β[0<->Λ_n^c] <= e^{-cn} for all n>=0." Printed in bond language.

§1.2 "Comments and consequences", paragraph "Site percolation" (verbatim):
"As in [AB87], the proof may be adapted to site percolation on transitive
graphs. In this context, one can obtain the inequality P_p[0<->∞] >=
(1/(d-1))·(p-p_c)/(1-p_c) (d is the degree of G) for p>=p_c by introducing
φ_p(S)=Σ_{x∈S}Σ_{y∉S,{x,y}∈E} P_p[0<->^S x]."

Precision finding: §1.2 is an adaptation remark whose PRINTED site inequality
is the supercritical mean-field bound with the modified φ_p; the subcritical
one-arm exponential decay for site is obtained by running the same adapted
proof, not printed as a numbered site theorem. It is not a bond result copied
to site — the φ_p modification is genuinely site-specific — and [AB87]
(Aizenman–Barsky 1987) inside that sentence, plus Menshikov (1986), cover site
subcritical exponential decay on finite-range transitive graphs directly.
#718 already flags "the main theorem there is printed in bond language; §1.2 is
essential here" and "actual site sharpness on the degree-eight graph, not a
bond result". KEEP, with a precision upgrade suggested in §5.

### 1.4 Matching duality for the p_c relation — PRIMARY_TEXT_READ [LIT]

This is the load-bearing citation and the direction must be stated correctly.
The elementary crossing duality for a planar matching pair gives only
p_c(G)+p_c(G*) <= 1 (closed matching crossings are subcritical exactly when
1-p < p_c(G*)). The SUPERCRITICAL half of #718 §4 needs the opposite inclusion:
p > p_c(NN) must force 1-p < p_c(NN+NNN), i.e. p_c(NN)+p_c(NN+NNN) >= 1 — the
Sykes–Essam equality.

Verified chain (arXiv:2205.02734v3, HTML intro, = RSA 65 (2024) 832–856,
DOI 10.1002/rsa.21226): Eq (1.1) is Sykes–Essam motivation ("presented
motivation ... verified in a number of cases when G is amenable"), NOT that
paper's theorem; Eq (1.3) p_u^site(G)+p_c^site(G*)=1 IS proved, in the
companion Grimmett–Li "Hyperbolic site percolation" (arXiv:2203.00981, accepted
RSA) for quasi-transitive one-ended plane matching pairs derived from a
mosaic; the intro states "When G is amenable, we have p_c^site(G)=p_u^site(G)";
G = square lattice Z^2 (amenable, transitive, one-ended) with G* = Z^2 plus
face diagonals is the paper's Figure 1.1 pair. Hence
p_c^site(NN)+p_c^site(NN+NNN)=1 holds for exactly this pair. This matches #613's
§(H3) reading of the published PDF, including its warning (van den Berg 1981
counterexample in wider generality — do not carry (H3) beyond the amenable
case). #718 cites "Eqs. (1.1),(1.3) and the amenable p_u=p_c discussion", which
reproduces the provenance chain without overclaiming (1.1) as a theorem of the
RSA paper. KEEP. Upgrade suggested in §5: cite arXiv:2203.00981 directly.

### 1.5 Harris/FKG for the P_0 lower bound — verified

#718 §2 supplies the induction proof of association for decreasing indicators
and explicitly does not pretend the overlapping boxes A_v are independent; the
overlap is handled by (3): intersection_v A_v^c ⊆ {r=0}, via the winding-lift
escape argument under the embedding condition 2R+2<=w<=m (both graphs have all
edges of sup-norm length <=1, so the same condition prevents edge
identification and spurious boundary adjacency for the NN+NNN box). FKG then
gives P_0 >= (1-a_R)^(wm). KEEP. The exact-rational FKG lower bound appears in
the executed script output (see §3).

## 2. Independent re-derivation of the proof skeleton

Items checked by hand, no gap found:

- Floor bound: floor(m/(R+1)) >= m/(2(R+1)) holds iff m >= 2(R+1), which the
  hypotheses give (m >= w >= 2R+2). With R=floor(w/8), 2R+2 <= w for w >= 8,
  and the L=floor(m/(R+1)) slabs k(R+1)..k(R+1)+R fit the m rows.
- Rate split (7): exponent/m <= floor(m/(R+1))/m · log(w a_R) - w log(1-a_R),
  so gamma_w = -log(w a_R)/(2(R+1)) + w log(1-a_R) with the second term
  NEGATIVE and retained; (8) is exactly the exponentiation of gamma_w > 0.
- §4 limits: a_R <= C e^{-cR} gives -log(w a_R)/(2(R+1)) -> c/2 (R ~ w/8) and
  w log(1-a_R) ~ -w C e^{-cw/8} -> 0, so gamma_w >= c/4 eventually, uniformly
  in m >= w. Same on the matching graph (finite-range transitive; diagonal
  steps change y by at most 1, so the slab containment (5) survives).
- Supercritical half: complementation interchanges P_0 and P_2 through the
  rank-sum identity r_NN(ω)+r_match(ω^c)=2 (repository digital Alexander,
  assumed as in #613's (H1)); the matching graph at 1-p is subcritical exactly
  when p > p_c(NN) by the equality of §1.4.
- Uniqueness: Margulis–Russo pivotality gives P_0' < 0, P_2' > 0, M' > 0 on
  (0,1); M(0)=-1, M(1)=+1; unique zero trapped for every m >= w.
- Fisher identity (11): re-derived symbolically from
  P_0=E(1-H), P_2=EH, P_1=1-E; the cross terms 2E'H' cancel, giving
  (E')^2/[E(1-E)] + E(H')^2/[H(1-H)]; at H=1/2 the second term is (M')^2/E
  since M'=2EH' there. Confirmed.
- H quantiles: H <= P_2/P_0 below p_c and 1-H <= P_0/P_2 above p_c trap all
  interior quantiles uniformly in m >= w; H is not F (F=(1+M)/2, different
  function; §5 shows the limits decouple on the thin sequence).

## 3. Executed controls (reproduced on this machine)

`git worktree` of pr718; pytest unavailable locally, run via unittest:

- `tests/test_rectangular_rank_odds.py`: 8 tests, OK.
- `tests/test_p398_tl_fattening.py`: 8 tests, OK (out of scope, ran anyway).
- `scripts/rectangular_rank_odds.py`: reproduced §7 — 65,536 configurations
  per graph (131,072 graph/configuration evaluations) on the 4x4 torus,
  complementary rank sum checked on all 65,536 pairings, 236 structural
  graph cases, exact probability versions of (4), (6), (11) at p=1/10, 1/3,
  1/2 with `inequalities_exact: true` and `exact_decomposition: true`,
  elapsed ~1.1 s.

These are finite sign checks; the note correctly does not claim the asymptotic
theorem from them (§7 of the note says the same).

## 4. KEEP / CORRECT / DROP

| # | Claim / citation in #718 | Verdict | Basis |
|---|---|---|---|
| 1 | uniform-in-m root consistency as w->infinity (theorem (1)) | KEEP | skeleton re-derived end-to-end; the needed matching-pair equality is rigorously sourced (§1.4); no gap found |
| 2 | H is not the birth-time mixture F | KEEP | different functions ((1+M)/2 vs P_2/(P_0+P_2)); #718 §5 shows the limits decouple (F->1/2, H->step); not proposed as an independent observation |
| 3 | #613 remains valid for the full law | KEEP | #613 Theorem L is scoped under (H0); #718 claims neither necessity of (H0) nor a sharp full-law condition; nothing here weakens #613 |
| 4 | no tilted-torus claim | KEEP | skew quotients explicitly excluded ("does not extend the geometry to arbitrary skew quotients without an additional strip construction"); nothing in the proof implies one |
| 5 | bond sharpness cited for site? | NO FLAG — KEEP | the note flags the bond language itself and relies on the §1.2 site adaptation plus [AB87]; see the precision upgrade below |

No DROP entries. No CORRECT entries of substance; three non-blocking
precision upgrades follow.

## 5. Non-blocking precision upgrades

1. Quote the §1.2 "Site percolation" passage (done in §1.3 above) when
   relying on the site one-arm decay, or cite Aizenman–Barsky (1987) /
   Menshikov (1986) directly for subcritical site exponential decay on
   finite-range transitive graphs — §1.2's printed site inequality is the
   supercritical mean-field bound, the subcritical decay being an adaptation
   the reader runs.
2. Cite Grimmett–Li arXiv:2203.00981 alongside rsa.21226 for (1.3); the RSA
   paper's (1.1) is motivation, not theorem — exactly as #613's (H3) warning
   already records.
3. In §5, name #716 §6's hypothesis (w_m = o(log m)) explicitly when invoking
   the super-thin sequence w_j=j, m_j=exp(j^2); "the previous full/empty-row
   bounds" currently leaves the hypothesis implicit.

## 6. Sources and read marks

- pr718 `notes/aspect-uniform-balance-consistency-20260912.md`:
  PRIMARY_TEXT_READ (branch pr718, full text).
- pr716 `notes/thin-torus-two-birth-limits-20260912.md`:
  PRIMARY_TEXT_READ (branch pr716, full text).
- main `notes/p613-quantile-convergence-20260907.md`: PRIMARY_TEXT_READ.
- Duminil-Copin & Tassion, arXiv:1502.03050v3 (HTML), Theorem 1.1 and §1.2:
  PRIMARY_TEXT_READ. [LIT]
- Grimmett & Li, "Percolation critical probabilities of matching
  lattice-pairs", RSA 65 (2024) 832–856, DOI 10.1002/rsa.21226; read as
  arXiv:2205.02734v3 (HTML introduction, Eqs (1.1)–(1.4), amenability
  passage): PRIMARY_TEXT_READ; the Wiley page itself: ABSTRACT_ONLY.
  [LIT]
- Grimmett & Li, "Hyperbolic site percolation", arXiv:2203.00981v3:
  ABSTRACT_ONLY — the p_u+p_c(G*)=1 theorem verified via its quotation as
  Eq (1.3) in the primary text and via #613's published-PDF reading, not from
  the companion's full text. [LIT]
- Background, not read in this round: [LIT] Aizenman–Barsky, CMP 108 (1987)
  489–526; Menshikov (1986); van den Berg, JMP 22 (1981) 152–157;
  Burton–Keane, CMP 121 (1989) 501–505; Sykes–Essam, J. Math. Phys. 5 (1964)
  1117–1127.

## 7. Verdict

Not falsified. KEEP x5, CORRECT 0, DROP 0, three non-blocking precision
upgrades. The aspect-uniform root theorem's citation base (site one-arm decay
via the DC–T adaptation; matching-pair equality via the Grimmett–Li pair;
Harris/FKG with overlapping boxes handled honestly) supports the claim as
stated; #613's full-law theorem and the root-vs-law distinction are preserved.
