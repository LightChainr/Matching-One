# Literature officer — Picco–Ribault spectrum, X θ(p_c), 2026-09-05

Two items not in earlier notes. Theory input; **does not enter** `docs/STATUS.md`.

Do not re-read as new: Ziff 2006, Camia–Feng, He, SKZ, Diamantis–Kleban, Ridout `t`, Vasseur `F(r)`, Grimmett–Li, Jacobsen 2015, Nolin backbone, Liu–Sun 2410.12724, Diskin–Tassion sharpness.

---

## Q1 — Picco–Ribault: the bootstrapped connectivities have **no even spin**

**Picco, Ribault, Santachiara,** [arXiv:1607.07224](https://arxiv.org/abs/1607.07224), *SciPost Phys.* **1**, 009 (2016). Numerical conformal bootstrap for percolation / Potts four-point functions. Spectrum ansatz: infinite, discrete, **non-diagonal** Virasoro.

They compute the cluster-connectivity four-points `R_σ` (linear combinations of `P_0,P_1,P_2,P_3`) with spectrum

> `ℳ_{2ℤ, ℤ+1/2} = {(Δ_{(r,s)}, Δ_{(r,−s)})}_{r∈ 2ℤ, s∈ℤ+1/2}`

which they state has **no even spin**. Leading state `(Δ_{(0,1/2)}, Δ_{(0,1/2)}) = (5/96, 5/96)` — the density/spin field. First few states at `c=0`: spin 0 `(5/96,5/96)`, spin `±1` `(39/32, 7/32)` and conjugate, scalar `(77/96,77/96)`, …. No `x=21/4`. No spin-4.

A second ansatz `ℳ_{2ℤ, ℤ}` *does* contain even spin; it is **not** the spectrum they use for the connectivity four-points that match Monte Carlo.

**Why this is the colour/spin split, from the other side.** Tan 2019 already has `x=21/4` as the scalar `P_{4s}` (spin 0), not spin-4. Picco–Ribault says the four-point *connectivities* that have been bootstrapped live in a channel with no even spin at all. A matching-odd leftover that looked like weight-4 holomorphic (`E₄`, spin 4) would not be this connectivity four-point. Colour-decompose the engine readout before any `21/4` or spin-4 fit. The even-spin ansatz exists on paper and has not been fitted to matching-odd.

Delfino–Viti three-point connectivity ([arXiv:1009.1314](https://arxiv.org/abs/1009.1314), imaginary DOZZ) was proved for triangular site by Ang–Cai–Sun–Wu (LQG/CLE; UCSD talk 3 Oct 2024). Three-point of the density field. Not matching-odd, not torus wrapping, not `21/4`.

---

## X, 2026-09-03 — Komargodski on `θ(p_c)=0`

[Zohar Komargodski, 3 Sep 2026](https://x.com/ZoharKo/status/2095513556942721254) (~8k views): physicists’ `θ(p_c)=0` is “now solved”; the 2-d exponent `β=5/36` is still not a theorem; AI may close physics/math gaps.

**2-d is not news.** Harris 1960 / Kesten 1980: square *bond* `p_c=1/2` and `θ(1/2)=0`. Square *site* has `θ(p_c)=0` by the same 2-d theory (RSW + Russo–Seymour–Welsh); it is not an algebraic statement about the value of `p_c`.

**3 ≤ d ≤ 10 is not our ticket.** [arXiv:2511.01851](https://arxiv.org/abs/2511.01851) v2 (1 Mar 2026) still lists continuity at `p_c` as open in those dimensions. A Lean 4 formalisation circulating as Anthropic `formal-math/percolation` claims all `d≥2` via Kozma–Nitzan [arXiv:2401.12397](https://arxiv.org/abs/2401.12397) (a *reduction to a conjectured inequality*). Until that inequality is a refereed theorem, do not cite `θ(p_c)=0` in 3-d as proved, and do not put it in STATUS. It does not algebraise square-site `p_c` and does not touch Q1–Q4 / P2 / #566 / #576.

X otherwise: still no researcher thread on algebraic square-site `p_c`, wrapping polynomials, or LCFT spin-4.

---

## Opinions

1. **Q1.** Picco–Ribault is the spectrum of the connectivities that exist. It has no even spin. Do not fit matching-odd to `21/4` or to spin-4 without a colour decomposition that puts the leftover in `ℳ_{2ℤ,ℤ}` rather than `ℳ_{2ℤ,ℤ+1/2}`.
2. **Q2.** Unchanged: Ridout-`t` / Cardy–Watts / `pinson_pi10_ratio`, or non-claim. Picco’s no-even-spin 4-point is another reason the holomorphic-`E₄` row was the wrong first test.
3. **P2 / #566.** Unchanged. Komargodski’s thread is 3-d continuity, not a square-site algebraic form.
4. **X.** Worth watching Komargodski / Deng / Sun. This week’s hit is not ours.

## Not established

- matching-odd ∈ `ℳ_{2ℤ,ℤ}`;
- `θ(p_c)=0` in 3-d as a refereed theorem;
- anything in the claim ledger.
