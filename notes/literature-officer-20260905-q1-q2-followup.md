# Literature officer — Q1/Q2 follow-up, 2026-09-05

Repo search for `2508.16047`, `0905.1727`, `Diamantis` returned **zero**. Two primary sources that change what Q1 and Q2 should spend output on. Theory input; **does not enter** `docs/STATUS.md`.

X this pass: no posts from Ziff / Deng / Jacobsen / Nolin / Garban on exact thresholds, wrapping, or LCFT. Quanta/Diskin already noted.

---

## Q1 — the only constructed log pair is still energy–hull

**Camia and Feng,** [arXiv:2508.16047](https://arxiv.org/abs/2508.16047) v2 (1 June 2026). *The percolation energy field and its logarithmic partner.* Triangular **site**. No square-site statement. No `x=21/4`. No spin-4.

They define two lattice fields and prove that the scaling limits of their two- and three-point functions exist and have the structure of an LCFT logarithmic pair (abstract):

> One of the two fields can be identified with the percolation analog of the Ising energy field, while the other is related to the percolation four-arm event.

The energy field is normalized by `a^{5/4} log |a|`; the two-point of the log partner scales as `|z|^{-5/2}` (twice the four-arm exponent `5/4`). That is the same pair He 2024 computed in the continuum (`γ^perco = −5/4`, energy–hull) and the same pair the q1-he note already said does **not** answer Astra Q1 (level-4, spin-4, `x=21/4`).

What is new is the *status*: this is now a theorem about lattice fields, not a CFT calculation. Q1 remains open. The thing a paid query would be buying — a level-4 spin-4 pairing fixed by `μ=−5/4` — is still not in the literature, and the one log pair that *is* constructed sits at `x=5/4`, scalar four-arm, triangular site.

**Do not** fit matching-odd to `21/4` on the back of this paper. Tan 2019 already has `21/4` as spin-0 (`P_{4s}`); Camia–Feng does not move that.

---

## Q2 — second-order modular forms are weight 0, not `E4`

**Diamantis and Kleban,** [arXiv:0905.1727](https://arxiv.org/abs/0905.1727) (2009). *New percolation crossing formulas and second-order modular forms.* Extends Kleban–Zagier [math-ph/0209023](https://arxiv.org/abs/math-ph/0209023).

Three crossing probabilities of Simmons–Kleban–Ziff:

1. may be written in terms of Cardy’s `Π_h` and Watts’ `Π_{h\bar v}` only;
2. are weakly holomorphic **second-order modular forms of weight 0** (type `(1, χ)`) on `Γ(2)`;
3. under boundedness-at-cusps assumptions, are **completely determined by their transformation laws**, with Cardy as the only physical input.

“Second-order” here means the transformation law involves a first-order form (an iterated integral / Rankin–Cohen, not a holomorphic Eisenstein series). Kleban–Zagier already had `Π_h'` first-order and `Π_{h\bar v}'` second-order. This paper puts the three new crossings in the same space.

**Why this is the Q2 alternative.** The N=290 fingerprint tested the conjunction “weight-4 holomorphic shape `g₂(τ)` AND a normalization that leaves `Ê4(2i)/Ê4(i)=11/4`”. Measured `1.880 ± 0.177`, excluded `11/4` at `4.9σ`. If the additive shape `A~(τ)` of a crossing-type (or wrapping-type, or Cardy-descendant) observable lives in the Diamantis–Kleban space, it is **weight 0, second-order**, not weight 4 holomorphic. Then excluding `11/4` is the expected outcome, not a surprise, and the next freeze should name a `Γ(2)` second-order competitor rather than another Eisenstein ratio.

This does **not** identify matching-odd as a crossing probability. It names the function space that wrapping / Cardy-type observables actually occupy, which the fingerprint did not include. Combine with the #576 non-claim: matching-odd is not Pinson `π({1,0})`; if a later readout *is* wrapping-flavoured, score it in this space, not against `E4`.

Kleban–Zagier was named in the 2026-09-05 brief and never opened. `0905.1727` is the paper that makes the function space precise.

---

## P2, one line

Akhunzhanov cites Jacobsen–Scullard graph-polynomial `p_c = 0.59274605079210(2)` as the most accurate square-site *estimate*. That is an estimator root, not an exact algebraic relation. Keep it out of the exact-`p_c` table. Unchanged from the first brief; restated because the 2022 paper is now in the wrapping stack.

---

## Opinions for subsequent analysis

1. **Q1.** Do not pay a query to re-derive energy–hull. Camia–Feng 2026 and He 2024 already have it. The query is still only worth running if it can say whether the level-4 spin-4 pair is fixed by `μ=−5/4` — and the 2026 lattice construction still does not.
2. **Colour, before any `21/4` fit.** Unchanged (Tan scalar; Tassion mono `>` poly).
3. **Q2 / next freeze.** Add a `Γ(2)` weight-0 second-order competitor, or write that matching-odd is not claimed to live in that space. Leaving `E4` as the only modular shape is how N=290 excluded everything on a list that did not contain the actual crossing-type functions.
4. **#576 / #567.** `pinson_pi10_ratio` (2.969 at `r=2`) and a Diamantis–Kleban-type shape are different competitors. Name both or non-claim both. Do not fold them into `11/4`.
5. **P2 / #566.** Unchanged. Height 6 is Wierman; square site has no exact form in Ziff 2006.

## Not established

- that matching-odd is a second-order modular form;
- a square-site analog of Camia–Feng’s lattice fields;
- anything in the claim ledger.
