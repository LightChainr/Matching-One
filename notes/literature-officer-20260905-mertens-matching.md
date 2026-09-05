# Literature officer — Mertens–Ziff matching identity, 2026-09-05

**Role:** theory input. Companion to the other `notes/literature-officer-*` on this branch and #571. **Does not enter** `docs/STATUS.md`.

Repo search for `1603.07289`, `Mertens matching function`, `2410.23250` returned **zero**. This is the project's lattice pair (square site / NN+NNN) on the torus, as a published identity.

---

## Ticket map

| Ticket | What this pass changes |
|---|---|
| matching-odd / engine | **Mertens–Ziff, J. Phys. A,** [1603.07289](https://arxiv.org/abs/1603.07289): exact finite-size Sykes–Essam. Cluster-count difference on the pair equals a wrapping-probability difference. Checkable at every `L` with no new continuum input. |
| #567 / two `11/4`s | Same paper: the unique root of the matching function converges empirically `~ L^{-4}`, **faster** than Newman–Ziff wrapping estimators `L^{-2.75} = L^{-11/4}`. A third `11/4` (estimator exponent) is still not the modular ratio. |
| colour / Q1 | **Radhakrishnan–Tassion,** [2410.23250](https://arxiv.org/abs/2410.23250) (v2 Aug 2025): monochromatic arm exponents are strictly larger than polychromatic; two-arm strictly larger than twice one-arm (Garban–Steif question). Polychromatic 8-arm is `21/4`. Same-colour 8-arm is **larger**. |

X this round: Quanta on Diskin–Tassion (already in the He/Q3 note). No new percolation-research posts.

---

## 1. Finite-size Sykes–Essam is an identity on our pair

**Mertens and Ziff,** [arXiv:1603.07289](https://arxiv.org/abs/1603.07289). Square **site** percolation; matching lattice = square with NN+NNN (Sq8). Matching polynomial, eq. (2):

```text
χ_□(p) = p − 2 p² + p⁴
```

Finite-size identity on an `L × L` torus, eqs. (4), (12):

```text
N_L(p) − N̂_L(1−p) − L² χ(p)  =  R_L^x(p) − R̂_L^x(1−p)
```

`N` is mean cluster count, hat is matching, `x` runs over wrapping types (both axes, either, horizontal, cross). The left-hand side is the matching function `M_L(p)`. The right-hand side is exactly a wrapping difference on the pair the engine already samples.

They record that the **both-axes** case is Scullard–Jacobsen's critical polynomial (an estimator, not an identity for unsolved `p_c`). The identity above is not an estimator: it holds at every `p` and every `L`.

**Plumbing.** At `L ≤ 12` the left-hand side can be evaluated from Akhunzhanov wrapping polynomials (deeper note) plus exact cluster enumeration, or from the transfer engine. If it fails, the wrapping channel is wrong. If it holds, later Pinson comparisons sit on a checked channel. Cheaper than another 200 M-sample rung.

**The unique root.** `M_L` has a unique root `p*_L ∈ (0,1)` with empirical

```text
p*_L − p_c  ∼  L^{-w},   w ≈ 4
```

against wrapping-on-primary-alone `L^{-2.75}`. Do not treat this `w ≈ 4` as a modular weight. Do not put `χ_□` in the P2 exact-`p_c` table: it is the matching polynomial, not a threshold relation.

**Matching-odd.** The RHS is a wrapping difference on (square, matching). Matching-odd in this repository is orientation-sensitive, not `R − R̂`. The identity still constrains the *even* wrapping pair. Name it as a non-claim or a control on #567, not as the odd slope.

---

## 2. Colour — monochromatic is strictly heavier

**Radhakrishnan and Tassion,** [arXiv:2410.23250](https://arxiv.org/abs/2410.23250) v2 (26 Aug 2025). Square **bond**. Two theorems:

1. two-arm exponent `> 2 ×` one-arm (quantitative Harris–FKG). Answers a question of Garban–Steif on exceptional times.
2. monochromatic `j`-arm exponent `>` polychromatic `j`-arm. New proof of Beffara–Nolin; framed as quantitative Reimer.

Polychromatic 8-arm is `(8²−1)/12 = 21/4`. If matching-odd wrapping is same-colour, the continuum exponent is **strictly larger than** `21/4`, and there is still no formula (Nolin–Sun: mono `k>2` open). Fitting `21/4` is a polychromatic identification. Colour-decompose before the fit (q1-noise note).

---

## 3. Opinions for subsequent analysis

1. **Engine / #567.** Evaluate `M_L` at a few committed sizes against the Mertens–Ziff identity. Zero-new-theory plumbing of the wrapping pair.
2. **#567 competitors.** Add `p*_L` (matching-function root) as a *threshold estimator*, not as a modular shape. Keep Newman–Ziff `L^{-11/4}` labelled as an estimator exponent. Pinson/Arguin or an explicit non-claim, unchanged.
3. **Colour, before N=580 or any `x=21/4` fit.** Mono `>` poly is now a theorem on square bond; expected on square site. Same-colour 8-arm is not `21/4`.
4. **P2 / Q4 / #566.** Unchanged. `χ_□` is not a `p_c` polynomial.
5. **Q1.** Unchanged (He does not compute the level-4 pairing; Tan `21/4` is spin 0).
6. **#227.** Tassion's two-arm `> 2×` one-arm is the exceptional-times input GPS asked for. Still: GPS on the even wrapping first, then noise on the odd channel.

## Not established

- that matching-odd *is* `R − R̂`;
- a published numerical Pinson value at `r=2,4` (still compute from the formula);
- a monochromatic 8-arm exponent;
- anything in the claim ledger.
