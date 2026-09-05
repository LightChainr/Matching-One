# #574 quote check — 2026-09-05

Primary-source check of the three items in [issue #574](https://github.com/LightChainr/Matching-One/issues/574). Theory input; does not enter `docs/STATUS.md`.

Sources opened: [cond-mat/0510245](https://arxiv.org/abs/cond-mat/0510245) (HTML + abs + APS landing), [cond-mat/9811416](https://arxiv.org/abs/cond-mat/9811416) (HTML), [cond-mat/0610813](https://arxiv.org/abs/cond-mat/0610813) (HTML + abs).

---

## Q1 — Ziff, PRE 73, 016134 (2006) — MATCH

**Citation is correct.** APS: *Phys. Rev. E* **73**, 016134 (published 27 Jan 2006). arXiv: `cond-mat/0510245`. Title on APS is *Generalized cell–dual-cell transformation and exact thresholds for percolation*.

**Polynomial: yes.** Abstract and body both print

```text
p^5 - 4 p^4 + 3 p^3 + 2 p^2 - 1 = 0,   p_c = 0.625457…
```

for the “A” lattice (bond). Same string in the APS abstract.

**Limitation sentence: MATCH, verbatim.**

> Unfortunately, the method does not appear to work for some of the more notorious unsolved systems: site percolation on the square and honeycomb lattices, and bond percolation on the kagomé lattice.

Safe to attribute to Ziff as his own statement that square site is outside the 2006 cell/dual-cell method.

---

## Q2 — Suding & Ziff, PRE 60, 275 (1999) — CLOSE, not character-exact

Printed sentence (HTML of `cond-mat/9811416`):

> We did not consider the square, triangular, and Kagomé lattices, as p_c (site) is either known exactly (triangular and Kagomé [14]), or has already been measured to a high degree of precision (square [15]) for these cases.

Ticket draft omitted `[14]`, `[15]`, and the trailing “for these cases.” Sense is unchanged: they skip square site because it is already a high-precision *measurement*, not because they have an exact algebraic form. **No exact algebraic claim for square-site `p_c` in this paper.** Quote the printed sentence, not the shortened draft.

---

## Q3 — Ziff & Scullard, J. Phys. A 39, 15083 (2006) — one height-6, not two height-4s

Authors on the arXiv abs: Robert M. Ziff, Christian R. Scullard. Journal ref: *J. Phys. A: Math. Gen.* **39** (2006) 15083–15090. `cond-mat/0610813`.

**Wierman bow-tie: yes, printed.**

> Setting the single bonds, r_1 and r_2, to p gives the condition 1 − p − 6p² + 6p³ − p^5 = 0, with solution in [0,1] p_c = 0.404518 … as found by Wierman.

Coefficient height 6 is on the page. Table 3 row (a) reprints the same equation, citing Wierman [13].

**bow-tie (d) vs martini-A: two drawings, one relation.**

The paper lists both:

| table | name | p_c | printed equation |
|---|---|---:|---|
| Table 1 | martini-A | 0.625457… | `p^5 − 4p^4 + 3p^3 + 2p^2 − 1 = 0` |
| Table 3 (d) | bow-tie (d) | 0.625457 | `1 − 2p^2 − 3p^3 + 4p^4 − p^5 = 0` |

These are the **same polynomial up to an overall sign**. The text does not treat them as independent solutions:

> …make the basic triangular cell the “A” generator (Fig. 3(b)) which produces the lattice shown in Fig. 4(d), with a transition point identical to that of the regular “A” lattice (given in Table 1) since no double bond arises in this case.

So: one generator, two embeddings, one threshold relation. **Do not print a P2 table with two independent height-4 data points.** Cite the A / martini-A relation once; mention Fig. 4(d) as another lattice in the same class if a figure is wanted.

The two *other* generalised bow-ties in Table 3 are independent:

- (b) deg 11, `p_c = 0.533213`, height 36
- (c) deg 11, `p_c = 0.672929`, height 35

Rows (e)–(h) are duals (`1 − p_c` of (a)–(d)), not new polynomials.

---

## Manuscript consequence

| Claim | Use |
|---|---|
| Square site is outside the 2006 cell/dual-cell method | Ziff PRE 73, 016134, the “notorious unsolved” sentence, verbatim |
| Suding–Ziff does not give an exact square-site form | printed sentence with the two citations |
| Published exact-bond height reached 6 in 1984 | Wierman, reprinted as Table 3(a) |
| Height-4 exact relation | **one** row: A / martini-A (bow-tie (d) is the same) |
| Deg-11 examples | Table 3(b) and 3(c) only |

## X this pass

No current posts on algebraic / exact square-site `p_c`, nor on Pinson–Arguin wrapping. Older hits are Tarasevich 2022 wrapping *polynomials* (finite-L exact probabilities, already in the deeper note), not a closed form for `p_c`.
