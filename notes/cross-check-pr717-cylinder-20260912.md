# #720 Round B peer check of #717: Jacobsen 2015 / Mertens–Ziff 2016 cylinder claims

2026-09-12. Retrieval only. Peer check of `notes/lit-cylinder-tm-width4-20260912.md`
(PR #717); I did not write #717. Parent #650. No census, no transfer engine, no
Huawei. Additive on `main`; independent of #708–#716 and #718.

Every check below was performed by re-fetching the primary sources in this session
and comparing #717's quotations and equation numbers against them. Local copies used:
arXiv PDFs (pdftotext) and arXiv/ar5iv HTML, dated 2026-09-12.

## 0. What this check read, and its marking

| # | Source | Access route this session | Mark |
|---|---|---|---|
| S1 | Jacobsen 2015, *J. Phys. A* **48** 454003, arXiv:1507.03027v1 | arXiv PDF v1 (pdftotext) **and** arXiv HTML v1 (LaTeXML) | **PRIMARY_TEXT_READ** (§2, §4, §4.3, §6.1, §7.1, §8; Eqs (4),(8)–(13),(26),(28),(31)–(34),(42)–(44); Table 2) |
| S2 | Scullard & Jacobsen 2012, arXiv:1209.1451v1 | arXiv PDF v1 | **PRIMARY_TEXT_READ** (abstract; Eqs (4),(5) and surrounding text) |
| S4 | Scullard & Jacobsen 2015, arXiv:1511.04374v1 | arXiv PDF v1 | **PRIMARY_TEXT_READ** (§3; Eq. (2)) |
| S5 | Mertens & Ziff 2016, arXiv:1603.07289v2 | arXiv PDF v2 **and** ar5iv HTML | **PRIMARY_TEXT_READ** (abstract; Eqs (12),(15),(17)–(24),(31)–(41) region) |
| S6 | Akhunzhaev, Eserkepov & Tarasevich 2022, arXiv:2204.01517v1 | arXiv PDF v1 | **PRIMARY_TEXT_READ** (§1; Eqs (4)–(9)) |
| S7 | Yang & Zhou 2024 Comment, DOI 10.1088/1751-8121/ad4d2c | Crossref API abstract | **ABSTRACT_ONLY** |
| S8 | Jacobsen 2024 Reply, DOI 10.1088/1751-8121/ad4d33 | Crossref API abstract | **ABSTRACT_ONLY** (body not fetched this session; not re-tried) |

**Equation-numbering record for S1 (ticket item 2).** The ticket says "HTML (13)=(PDF
13); HTML (50) may be PDF (44)". Both resolved by direct comparison:

* HTML (13) = **PDF (13)** (main result `P_B(q,v) = 0 ⇔ Λ_open = Λ_closed`).
* HTML (50) = **PDF (44)** (`p_c(n) − p_c = O(n^{−4})`), exactly; not merely "may be".

The offset is a constant **+6** from the R-matrix section onward: HTML (32) = PDF (26);
HTML (34) = PDF (28); HTML (40) = PDF (34); HTML (48) = PDF (42); HTML (49) = PDF (43);
HTML (50) = PDF (44). Eqs (4), (8), (9), (10), (13) carry the same number in both
renderings. #717's §0 caveat is accurate.

**New: equation-numbering caveat for S5 that #717 did not record.** The ar5iv HTML and
the arXiv PDF v2 of 1603.07289 agree on (20)–(24) and (31) but **diverge by +2 in the
scaling section**: ar5iv (39) `p*_L − p_c ~ L^{2−x−1/ν}` = **PDF (37)**; ar5iv (40)
`L^{2−y−3/ν}` = **PDF (38)**; ar5iv (41) integral estimator `L^{−1.65}` = **PDF (39)**.
#717 cites the ar5iv numbers (39)/(40)/(41) while claiming PDF v2 as an access route.
This does not change any content verdict below, but #717 §0 should add an S5 caveat
mirroring its S1 one. **This is the only correction found.**

## 1. KEEP/CORRECT/DROP table (ticket-required)

| # | Claim in #717 | Verdict | Evidence checked this session |
|---|---|---|---|
| C1 | "(13) is m→∞" | **KEEP** | PDF (13) followed verbatim by "valid for a basis B of size n × m, with n finite and m → ∞". The proof (monotonicity + intermediate value theorem) is quoted accurately; Eq. (4) `P_B = Z_2D − qZ_0D` and the `Z_2D ~ (Λ_open)^m`, `Z_0D ~ (Λ_closed)^m` steps (PDF (11)–(12)) confirmed. |
| C2 | "(50) observed not proved" | **KEEP** | PDF (44) is introduced by "What we have observed in sections 4 and 7 is…"; PDF then says "It is clear that more work would be required to establish whether (43) [HTML (49)] can be shown—obviously using more ingredients—to actually imply (44) [HTML (50)]." #717's reading is exactly right: (50) is an observed law, not deduced from (49). |
| C3 | "no 5+15+16 in print" | **KEEP** | S1's only cylinder decomposition is the string-graded (9) `T~ = ⊕_{k=1}^n T^{(s=2k)} ⊕ T_open ⊕ T_closed` — six sectors at n=4, not 5+15+16. Greps of S1 and S5 full texts find no `B_5/B_15/B_16`, no `tr B_15^m − tr B_5^m`, no trace-difference identity. |
| C4 | "no c ρ^m/m in print" | **KEEP** | S1 §4.3: "the rate of convergence being exponential in m" — no rate, no closed form. S5: all displacements are power laws in L (`w ≈ 4`, primary-lattice `L^{−2.75}`, `L^{2−x−1/ν}` with w = −4.17, `L^{2−y−3/ν} ≈ L^{−1.55}`, integral `L^{−1.65}`). S6 Eq. (9): power-law series. #717 correctly marks this a negative retrieval result, not a no-go theorem. |
| C5 | "MZ matching-lattice ≠ digital-Alexander M" | **KEEP** | MZ's matching lattice is the Sykes–Essam matching lattice (square lattice plus face diagonals); their `M_L(p) = N_L(p) − N̂_L(1−p) − L^2 χ(p)` (PDF (15)) `= R^x_L(p) − R̂^x_L(1−p)` (PDF (20)), x ∈ {c,b,e,h}. It is a cluster-count/wrapping-probability difference on a *different graph*, not the repository's configuration-rank observable `M = P_2 − P_0` of #705/#710/#718. #717 §8 keeps the two uses distinct. |
| C6 | Invented quotations | **none found — KEEP all** | Every quotation in #717 §§1–5 was located verbatim in the fetched texts (log in §2). No reconstructed or fabricated quotation detected. |

**CORRECT (2 items, both editorial):**

* **N1.** Add the S5 ar5iv-vs-PDF-v2 numbering caveat to #717 §0 (detail in §0 above;
  ar5iv (39)/(40)/(41) = PDF v2 (37)/(38)/(39)).
* **N2.** Inside verbatim S1 quotations, #717 silently renders the equation references
  in HTML numbering — e.g. the "more work would be required" quote prints "(49)…(50)"
  where the PDF prints "(43)…(44)". This is consistent with #717's declared convention
  ("the HTML/LaTeXML numbers … Content, not number, is what is quoted") but it means
  those quotations are not byte-verbatim against the PDF. Recommend bracketing such
  substitutions, e.g. "(49) [PDF: (43)]", so no later reader counts this as a
  transcription error.

**DROP:** nothing. No claim of #717 required deletion.

## 2. Quotation verification log (spot map)

* S1 quote before (13): "For v >> 1 … intermediate value theorem implies our main
  result" — PDF verbatim (≪/≫ rendered `<<`/`>>` in #717).
* S1 (9)-context quote: "This is so precisely because contributions to Z_1D are
  excluded from (4), implying that loops winding around the cylinder carry the weight
  n_wind = 0." — PDF verbatim.
* S1 §4.3 quote: "As expected, the results converge rapidly to the m = ∞ limit…" —
  PDF verbatim; §4.3 attribution correct (Kagome n=2 subsection, where the sentence
  appears).
* S1 §8 quotes around (42)–(44) — PDF verbatim modulo the disclosed numbering
  substitution (N2).
* S1 §6.1 R-matrix: "Ř-matrix contains only two out of fourteen possible terms" —
  PDF verbatim; HTML (32) = PDF (26); §6.1 attribution correct.
* S1 Table 2 n=4 `0.5914171708531384817988341017359231779642` — read from PDF v1
  Table 2, row n=4; matches. `n_max = 21` sentence — PDF verbatim.
* S2: (4) normalization `P(0D;B)+P(1D;B)+P(2D;B)=1`; (5) `P(2D;B)=P(0D;B)`; "Despite
  its apparent simplicity, eq. (5) is the main result of this paper."; unique-root
  sentence; contraction–deletion equivalence left open — all PDF verbatim.
* S4: Eq. (2) `P_B = P_2D − q P_0D`; "…and we could then set P_2D = P_0D. However,
  this is wasteful because P_1D is never used for anything." — PDF verbatim.
* S5: (20) main result + "This is the main result of this paper."; (21) SJ condition
  and the two following sentences; (22) `M_L(p_c) = 0` (self-matching) "for all values
  of L"; (24) (self-dual); unique-root passage with `w ≈ 4` and `L^{−2.75}` — all
  verbatim in PDF v2 and ar5iv.
* S6: "up to L = 16 for a cylinder, and up to L = 12 for a torus"; Eqs (4)–(6)
  estimators; Eq. (9) series — PDF verbatim.
* S7/S8: Crossref abstracts match #717's verbatim blocks, including the disputed
  digits `0.592 746 050 792 10(2)` vs `0.592 746 050 7896(1)`.

Not re-checked this session: #717 §5's publisher search-index snippet for S8 (marked
`[LIT]` there); it remains a search-index artifact and no argument depends on it.

## 3. Tripwire check (ticket item 4)

#717 must not imply that `B_5`, `B_15`, `B_16`, `tr B_15^m − tr B_5^m`, or the
aspect-uniform `p_{w,m} → p_c` statement of #710/#718 are Jacobsen (13) or any
published result. Result: **pass**.

* #717 §0 "Tripwire honoured" and §1.2 explicitly state the objects do not appear in
  S1 and are not attributed to any source; §1.2, §4 and §7 say so again in
  substance ("Neither a 5+15+16 split nor a `tr B_15^m − tr B_5^m` identity appears
  anywhere in S1").
* #717 correctly reports (13) as an `m → ∞` eigenvalue identity over `s = 0` reduced
  states (count (8) = 35 at n = 4), not as a finite-m trace identity.
* #717 does not mention #718's `p_{w,m} → p_c`-without-aspect-bound theorem anywhere;
  no misattribution is possible from its text.

## 4. Verdict

#717 is accurate against the primary texts on every substantive point checked. Two
editorial corrections (N1: add the S5 numbering caveat; N2: bracket the numbering
substitutions inside S1 quotes). Nothing to DROP. No invented quotation found.
