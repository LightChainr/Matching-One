# Round B peer check of #724: thermal `m λ^m` vs Jordan/LCFT diagnostics

**Date:** 2026-09-12
**Ticket:** #727 (parent #650); peer check of PR #724 / note
`lit-thermal-jet-jordan-20260912.md` (ticket #714)
**Claim level:** peer-check status only; no repository claim is upgraded,
downgraded or moved
**This note does not enter** `docs/STATUS.md`
**Scope guard:** this note does **not** adjudicate #275, does not declare
#275 solved (or refuted), does not merge #724, and adds no new retrieval
beyond re-verification of #724's sources.

## Verdict in one line

**#724 survives Round B peer check.** Every load-bearing verbatim quote was
re-fetched on 2026-09-12 and matches the primary texts; the central
distinction — the published Jordan diagnostic is *weight degeneracy + field
mixing + a logarithm in a correlation function*, while `m λ^{m-1}` from a
linear split of a semisimple pair is ordinary calculus on analytic branches
and is **not** that diagnostic — is correctly drawn. Two minor
paraphrase/citation items are CORRECTed below; nothing is DROPped.

## What was (re)read for this check

| # | Source | Where read | Tag this round |
|---|---|---|---|
| 1 | Vasseur, Jacobsen, Saleur, [arXiv:1206.2312](https://arxiv.org/abs/1206.2312), J. Stat. Mech. **L07001** (2012) | ar5iv HTML re-fetched 2026-09-12: eqs. (8), (9), (11) verbatim + conclusion paragraph | PRIMARY_TEXT_READ |
| 2 | Gurarie, [arXiv:hep-th/9303160](https://arxiv.org/abs/hep-th/9303160), Nucl. Phys. B **410**, 535 (1993) | ar5iv HTML re-fetched 2026-09-12: eqs. (13)/(19), (22), Jordan-cell sentences | PRIMARY_TEXT_READ |
| 3 | Cardy, [arXiv:cond-mat/9911024](https://arxiv.org/abs/cond-mat/9911024) | ar5iv HTML re-fetched 2026-09-12: degeneracy mechanism, unitarity exclusion, `q`-derivative passage | PRIMARY_TEXT_READ |
| 4 | Bamieh, [arXiv:2002.05001](https://arxiv.org/abs/2002.05001) (v2, 2022) | arXiv HTML re-fetched 2026-09-12: semisimple hypothesis, §3.2, Example 1, eq. (25) | PRIMARY_TEXT_READ |
| 5 | #724 note in full (`git show pr724:notes/lit-thermal-jet-jordan-20260912.md`) | local branch `pr724`, commit `ef6bdcf` | PRIMARY_TEXT_READ |
| 6 | Creutzig, Ridout, [arXiv:1303.0847](https://arxiv.org/abs/1303.0847) | **not re-fetched this round**; accepted from #724's full-text reading | [LIT: not re-verified this round] |
| 7 | Jacobsen, [arXiv:1507.03027](https://arxiv.org/abs/1507.03027) | **not re-fetched this round**; accepted from #724's full-text reading | [LIT: not re-verified this round] |
| 8 | Mertens, Ziff, [arXiv:1603.07289](https://arxiv.org/abs/1603.07289) | **not re-fetched this round**; accepted from #724's full-text reading | [LIT: not re-verified this round] |
| 9 | Kato, *Perturbation Theory for Linear Operators*, Grundlehren **132** | not accessed this round; [LIT] status **inherited from #724** (body never verified) | [LIT: primary text not verified] |
| 10 | Pinson, J. Stat. Phys. **75**, 1167 (1994) | not accessed this round; [LIT] status inherited from #724 | [LIT: primary text not verified] |

## 1. Mandated re-fetch 1: Vasseur–Jacobsen–Saleur eq. (11) and the
Boltzmann-weight-derivative paragraph

Re-fetched from ar5iv on 2026-09-12. **Matches #724 verbatim.**

Eq. (11) and the Jordan-cell statement (Section 4):

> Logarithms in LCFTs can be associated with the non-diagonalizability of the
> scale transformation generator—i.e., the Hamiltonian in the usual radial
> quantization—of the theory. … at $Q=1$
>
> $$\tilde{\psi}_{ab}(\Lambda r)=\Lambda^{-5/4}\left(\tilde{\psi}_{ab}(r)+\frac{2\sqrt{3}}{\pi}\log\Lambda\ \varepsilon(r)\right). \tag{11}$$
>
> The field $\tilde{\psi}_{ab}$ is therefore mixed with the energy operator
> $\varepsilon(r)$ after a scale transformation. … In other words, the scale
> transformation generator (or Hamiltonian) is non-diagonalizable, with a
> rank-2 Jordan cell mixing the two fields $\tilde{\psi}_{ab}$ and
> $\varepsilon$.

Conclusion paragraph, verbatim (the exact paragraph #724 leans on for the
DROP side of the inference):

> In conclusion, it is important to stress that logarithmic terms such as
> those we have identified would not be present for generic $Q$, and occur
> solely because of the special degeneracies present at $Q=1$. This is of
> course quite different from logarithmic dependencies in other non-local
> quantities—see e.g. [27, 28]—**which are obtained as derivatives of
> correlation functions with respect to the Boltzmann weights (such as
> $Q$)**.

Supporting equations also re-verified:

```text
<ψ̃_ab(r) ψ̃_cd(0)> = 2A(1) r^{-5/2}[ … + (4√3/π) log r ]        (their eq. (8))
lim_{Q→1} (Δ_ψ̂ − Δ_ε)/(Q−1) = √3/π                              (their eq. (9))
```

**Assessment:** #724's characterization is exact. The published chain is
(i) degeneracy of two scaling weights at `Q=1` (`Δ_ε = Δ_ψ̂ = 5/4`),
(ii) **mixing** of the two fields under scale transformation (eq. (11)),
(iii) the logarithm in a **two-point function** (eq. (8)). At no point does
the published diagnostic run through a factor `m` in a parameter derivative
of a transfer-matrix trace.

## 2. Mandated re-fetch 2: Gurarie 1993 L0 Jordan cell

Re-fetched from ar5iv on 2026-09-12. **Matches #724.** His eqs. (13)/(19)
(the general form is (19); (13) is the `c=−2` instance):

```text
L0|C,n>   = (h_C+n)|C,n>
L0|C_1,n> = |C,n> + (h_C+n)|C_1,n>                              (his eq. (19))
```

> "Ordinary primary operators are known to be the eigen vectors of the `L0`
> operators, and their eigen values are the dimensions of these operators.
> It will be shown that those "new" operators, which I will call
> pseudo-operators, are the basis of the Jordan cell for `L0`."

Two-point functions (his eq. (22)):

```text
<C_1(z)C_1(w)> = −2/(z−w)^{2h_C} [log(z−w) + λ']
<C(z)C_1(w)>   = 1/(z−w)^{2h_C}
```

Criterion (end of his §2):

> "…we must include logarithmic operators in the theory if it possesses at
> least two operators the product of which when expanded according to the
> fusion rules … contains the contribution of at least two operators with
> the same dimension."

**Assessment:** the defining object is the Jordan cell of `L0` with the log
appearing in correlation functions — exactly as #724 records.

## 3. Additional spot-checks performed this round (beyond the mandate)

**Cardy 1999 — re-verified verbatim.** All three #724 quotes match:
the degeneracy/cancellation mechanism ("…the leading terms will cancel
leaving a logarithmic term proportional to `r^{−2x_i} ln r`"), the unitarity
exclusion ("Such operators should not occur in unitary conformal field
theories, such as correspond to pure critical systems with positive
Boltzmann weights…"), and the percolation passage ("The connectivities of
the percolation problem are given by the derivatives with respect to `q` at
`q=1`, and are finite. But in this case **there are no logarithmic terms of
the above form**."). **One nuance #724 does not mention:** Cardy's footnote
[17] to that same sentence reads "However, at least in `d=2`, logarithms do
arise at a higher level because the coefficient of `TT̄` in the OPE is
`O(1/c²)`." This does not weaken #724 — the higher-level logs are again
correlator/OPE logs inside LCFT, not Boltzmann-weight-derivative Jordan
evidence — but it belongs on the record. Logged as a CORRECT (annotation).

**Bamieh 2020/2022 — re-verified verbatim.** The semisimple hypothesis is
word-for-word as #724 quotes it ("Throughout this note, we will assume the
semi-simple case, i.e. that `A_ε` has a full set of eigenvectors (i.e.
diagonalizable) for each `ε` in some neighborhood of zero."), eq. (25)'s
scalar form `λ_1i = w*_0i A_1 v_0i` and Example 1
(`λ_ε1 = 1+εα, λ_ε2 = 1+εβ`) are exact, and — confirmed this round — **the
word "Jordan" does not appear anywhere in the document**. The degenerate
case is handled via `(Λ_1)_11 = (V*_0 A_1 V_0)_11`, whose diagonal entries
are the eigenvalues of the perturbation restricted to the degenerate
eigenspace. No Jordan block enters.

## 4. KEEP / CORRECT / DROP

| #724 item | Verdict | Evidence status |
|---|---|---|
| VJS eq. (11) + "rank-2 Jordan cell mixing the two fields" as the published diagnostic | **KEEP** | re-fetched 2026-09-12, verbatim match — PRIMARY_TEXT_READ |
| VJS conclusion paragraph: Boltzmann-weight-derivative logs are "quite different" from the Jordan-cell logs | **KEEP** | re-fetched 2026-09-12, verbatim match — PRIMARY_TEXT_READ |
| Gurarie L0 Jordan cell, eqs. (13)/(19), (22), pseudo-operator sentence, fusion criterion | **KEEP** | re-fetched 2026-09-12, verbatim match — PRIMARY_TEXT_READ |
| Cardy 1999: unitarity exclusion + percolation `q`-derivatives have "no logarithmic terms of the above form" | **KEEP** (with annotation) | re-fetched, verbatim match — PRIMARY_TEXT_READ; CORRECT: footnote [17] higher-level `TT̄` caveat should be noted |
| Core distinction: `m λ^{m-1}` from a linear split of a semisimple pair is **not** the published Jordan diagnostic | **KEEP** | the entire quoted diagnostic chain (VJS (8)/(9)/(11), Gurarie (19)/(22)) runs through degeneracy + mixing + log in correlators; Bamieh's semisimple split chain produces `2mλ^{m-1}` with no Jordan block — both sides re-verified this round |
| Bamieh semisimple hypothesis quote; "Jordan" absent from the document | **KEEP** | re-fetched, verbatim match — PRIMARY_TEXT_READ |
| "analytic branches `λ̄ + ε μ_j`" attributed to Bamieh §3.2 | **CORRECT (paraphrase, not quotation)** | Bamieh never uses "analytic branches" or `μ_j`; his notation is `(Λ_1)_11 = (V*_0 A_1 V_0)_11` with diagonal entries the eigenvalues of the restricted perturbation. Mathematical content identical; #724 should mark this as paraphrase, not verbatim. Re-verified — PRIMARY_TEXT_READ |
| Bamieh dated "2020" | **CORRECT (citation)** | v1 2020, v2 2022-04-24; cite as (2020, v2 2022) — metadata re-verified this round |
| Kato Ch. II §1.2/§1.6 contrast: linear split ⇒ semisimple, fractional-power split ⇒ defective | **KEEP as [LIT]** | not re-verified this round; #724 correctly did not quote the unaccessed body — [LIT: primary text not verified] |
| Jacobsen 2015: crossing of Perron eigenvalues used, never read as Jordan; `q4 = 0.59141…` is a finite-size estimate, not a new `pc` | **KEEP** | not re-fetched this round; accepted from #724's full-text reading — [LIT: not re-verified this round] |
| Mertens–Ziff 2016: no eigenvalue/Jordan reading of `ρ^m/m`; main result is `M_L(p) = R^x_L(p) − R̂^x_L(1−p)` | **KEEP** | not re-fetched this round; accepted from #724's full-text reading — [LIT: not re-verified this round] |
| Creutzig–Ridout modern Jordan-cell statement, eq. (1.10) | **KEEP** | not re-fetched this round — [LIT: not re-verified this round] |
| Pinson 1994 listed with abstract unverified | **KEEP as [LIT]** | status honestly declared in #724; unchanged — [LIT] |
| #710 toy `diag(λ+ε, λ−ε)` gives `2mλ^{m-1}` with no Jordan block | **KEEP** | elementary algebra, independently re-derived; consistent with #710's own wording ("Semisimple crossing, not a Jordan block") |
| Final table: "#275's `m λ^m` thermal jet is thereby resolved — **not claimed**" | **KEEP** | correct restraint; nothing in the re-verified literature settles the empirical identity of the #275 jet |

**DROP count: 0.** No claim in #724 required removal.

## 5. Not established (unchanged from #724, endorsed)

- The empirical identity of #275's thermal jet remains open; this note and
  #724 record only what the published diagnostics do and do not license.
- #275 is **not** solved by #724 and is not declared solved here.
- Kato's exact theorem wording ([LIT], body never accessed by either round).
- The repository's own `thermal-jordan-spin4-descendant.md` construction is
  neither confirmed nor refuted by any source quoted in either note.

## Sources (this round's own reads tagged)

1. arXiv:1206.2312 — Vasseur, Jacobsen, Saleur (2012) — PRIMARY_TEXT_READ (re-fetch 2026-09-12).
2. arXiv:hep-th/9303160 — Gurarie (1993) — PRIMARY_TEXT_READ (re-fetch 2026-09-12).
3. arXiv:cond-mat/9911024 — Cardy (1999) — PRIMARY_TEXT_READ (re-fetch 2026-09-12).
4. arXiv:2002.05001 — Bamieh (2020, v2 2022) — PRIMARY_TEXT_READ (re-fetch 2026-09-12).
5. `pr724` branch, commit `ef6bdcf`, file `notes/lit-thermal-jet-jordan-20260912.md` — PRIMARY_TEXT_READ.
6. arXiv:1303.0847, arXiv:1507.03027, arXiv:1603.07289 — [LIT: not re-verified this round].
7. Kato, Grundlehren 132 — [LIT: primary text not verified, both rounds].
8. Pinson, J. Stat. Phys. 75, 1167 (1994) — [LIT: primary text not verified, both rounds].
