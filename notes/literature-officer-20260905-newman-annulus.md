# Literature officer — Newman–Ziff four wrappings, annulus complex roots, 2026-09-05

Repo search for `0.521058290` / `2410.04767` returned **zero** (Pinson `R_h` was computed in the #576 note; the other three Newman–Ziff numbers were not frozen). Theory input; **does not enter** `docs/STATUS.md`.

X since 2026-06: empty for Deng / wrapping / exact threshold / LCFT.

Do not re-read as new: Pinson `π({1,0})` at `r=1,2,4`, Duncan–Kahle A vs S, Nolin backbone value, Picco no-even-spin, Jacobsen 2015 estimator.

---

## #576 — Newman–Ziff, four published numbers, one unspecified engine

**Newman and Ziff,** *Phys. Rev. E* **64**, 016706 (2001); the 2000 preprint [cond-mat/0005264](https://arxiv.org/abs/cond-mat/0005264). They quote Pinson’s infinite-volume wrapping probabilities at criticality:

| symbol | event (their words) | Pinson value |
|---|---|---|
| `R_h` | wraps horizontally (specified direction) | **0.521 058 290** |
| `R_e` | wraps either horizontally or vertically, or both | **0.690 473 725** |
| `R_b` | wraps in both directions | **0.351 642 855** |
| `R_1` | wraps around one specified axis but not the other | **0.169 415 435** |

Identities: `R_b = R_e − 2 R_1`, `R_h = R_e − R_1`. The #576 note’s `π({1,0})(r=1) = 0.521058290` is `R_h`. Duncan–Kahle **A** is “at least one giant 1-cycle” ≈ `R_e`; **S** (homology basis) ≈ `R_b`. Specified-direction A is `R_h`.

They then invert `R_L(p) = R_∞(p_c)` to estimate square-site `p_c = 0.59274621(13)` from `R_h`. That is an estimator, not an exact form (Jacobsen 2015 already classified square site as the finite-size-dependent class).

**#576 freeze.** The wrapping channel on the next ladder is one of these four, or it is matching-odd and not wrapping. Write which. Comparing an unnamed mixture to `E₄` `11/4` or to a single Pinson number is how the last fingerprint excluded a list that did not contain the observable.

---

## Q1 — the backbone CFT spectrum contains **complex** dimensions

**Sun, Xu, Zhuang,** [arXiv:2410.04767](https://arxiv.org/abs/2410.04767) v2 (9 Feb 2025). Exact annulus-crossing formulae. Third event: two disjoint open paths connecting the two boundaries. Leading exponent = backbone `x_B`, the unique real root in `(1/4, 2/3)` of

> `√(36x+3)/4 + sin(2π √(12x+1)/3) = 0`

other than `−1/12` and `1/4`. The same elementary equation has **countably many complex roots**. Those roots appear as the exponents of the subleading terms in the crossing formula. They conclude:

> This suggests that the backbone exponent is part of a CFT whose bulk spectrum contains this set of roots.

**What this does to a matching-odd log-slope.** Polychromatic arms are rational. Backbone (monochromatic two-arm) is transcendental, and its CFT companions are complex. Picco–Ribault’s connectivity four-point has no even spin. A single real Kac weight (`21/4`, `5/4`, `2`) is the wrong default for a leftover that has not been colour-decomposed. Do not fit `21/4` on the back of “some exponent in percolation is 21/4” (Tan’s `P_{4s}` is scalar, different channel).

Nolin–Qian–Sun–Zhuang PRL **134**, 117101 (2025) is the letter version of the same claim. Already named; the complex-root sentence is in the companion [2410.04767](https://arxiv.org/abs/2410.04767), which was not opened.

---

## One line, generic-Q Jordan

**Liu, Jacobsen, Saleur,** [arXiv:2403.19830](https://arxiv.org/abs/2403.19830). Lattice “emerging Jordan blocks” confirm rank-2 logarithmic structure of Potts at generic `Q`. At `c=0` they recall that the `T`–`t` block was already seen on the lattice (the Ridout `t` of the earlier note). They do **not** construct a spin-4 / rank-4 block at percolation.

---

## Opinions

1. **#576 / #567.** Freeze `R_h` / `R_e` / `R_b` / `R_1`. The numbers are in Newman–Ziff. Aspect `r=2,4` changes them (Pinson `π({1,0})` ratio 2.969 at `r=2`); the square-aspect table above is `r=1` only.
2. **Q1.** Colour first. The monochromatic-two-arm CFT is not a real Kac table. A paid query is still only the level-4 spin-4 question He left open.
3. **P2.** Newman–Ziff `p_c` from wrapping is an estimator. Keep out of the exact table.
4. **X.** Dead since June. arXiv remains the channel.

## Not established

- which of the four `R` the engine wrapping channel is;
- matching-odd ∈ backbone CFT spectrum;
- anything in the claim ledger.
