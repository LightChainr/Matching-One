# Literature officer — Q1 He primary, Q3, neighbourhoods, 2026-09-05

**Role:** theory input. Companion to the other `notes/literature-officer-*` on this branch and #571. **Does not enter** `docs/STATUS.md`.

Repo search for `Camia Feng`, `Diskin`, `Cieplucha`, `a0=-25/48` returned **zero**. He 2024 is cited in #571 as owner of `μ = -5/4`; this pass is the primary reading of what that paper actually computes.

---

## Ticket map

| Ticket | What this pass changes |
|---|---|
| Q1 | He, SciPost Phys. **19**, 008 (2025) = [2411.18696](https://arxiv.org/abs/2411.18696) computes `γ^perco = -5/4` for the **energy–hull** rank-2 pair and `b = -5` for `(t, T)`. It does **not** compute a coupling for `x=21/4` or spin 4. The Astra hope that `μ = -5/4` fixes the level-4 pairing is not a published result. |
| Q1 rank | Same paper, and He–Saleur [2109.05050](https://arxiv.org/abs/2109.05050): rank-3 at `T\bar T`, `a_0 = -25/48`. He 2025: Liouville-norm comparison “suggests the potential existence of **arbitrarily high rank** Jordan blocks.” That is a warning against treating the energy-level `γ` as determining descendants. |
| Q3 | Cardy [0907.4070](https://arxiv.org/abs/0907.4070): lattice parafermions on the square grid are **C4** discrete-holomorphic. Multiplication by `i` is the Cauchy–Riemann rotation. A local C3 character is not in that toolkit; the repo no-go is consistent with the literature, not a surprise. |
| matching / P2 | Malarz [2303.10423](https://arxiv.org/abs/2303.10423) and Ciepłucha–Utnicki–Wołoszyn–Malarz [2503.16703](https://arxiv.org/abs/2503.16703): square **site** with complex neighbourhoods (up to 6th/7th zone, including NN+NNN). Thresholds are Newman–Ziff MC, `p_c(ζ) ∝ ζ^{-0.545}`. Not algebraic. Do not put them in an exact-threshold table. |
| community | Diskin–Easo–Radhakrishnan–Sudakov–Tassion [2603.03257](https://arxiv.org/abs/2603.03257), Quanta 2026-08-31: supercritical sharpness on **every infinite transitive graph**. For `ℤ²` this is 1980s. Does not move square-site algebraicity, wrapping, or matching-odd. |

X, 2026-09-04: [Quanta](https://www.quantamagazine.org/stunning-percolation-proof-solves-decades-old-puzzle-about-phase-transitions-20260831/) on Diskin–Tassion. A 2026-09-01 post repeats a **rumour** that Levant Alpoge has a major percolation result (tba). Rumour is not a source.

---

## 1. Q1 — He 2025 does not answer the Astra question

Abstract, verbatim scope:

```text
We give a generic construction of logarithmic operators based on Kac
operators and focus on the rank-2 pair of the energy operator mixing
with the hull operator.
```

Computed couplings (primary quotes):

```text
γ^perco = 1 / (2 (h'_{2,1} − h'_{0,2})) = −5/4
b = −1 / (2 h'_{1,2}) = −5
```

`γ` is energy `ε = Φ̂_{2,1}` mixing with 2-hull `Φ̂_{0,2}`. `b` is the `(t, T)` pair. Neither is the level-4 `x=21/4`, `s=4` pairing Q1 asks about. Search of the paper for `21/4` and spin-4: **absent**.

He–Saleur, [2109.05050](https://arxiv.org/abs/2109.05050): identity module at `c=0` has a rank-3 Jordan cell at `h = \bar h = 2` involving `T\bar T`, with non-chiral coupling `a_0 = -25/48`, identical for percolation and polymers.

He 2025 on higher rank:

```text
comparison with c < 1 Liouville CFT suggests the potential existence
of arbitrarily high rank Jordan blocks
```

and, if cluster Kac norms match Liouville, “rank-4, rank-5, and even arbitrarily higher rank Jordan blocks whose bottom fields are identified with these Kac operators.”

**What this does to Q1.** Paying an external model to ask whether `μ = -5/4` *fixes* the level-4 spin-4 coupling is still the right question, because the published computation stops at energy/hull and at `T\bar T`. Do not treat He 2025 as already answering it. Do not fit matching-odd slopes to `−5/4` or to `−25/48`.

**Camia–Feng, JHEP 08 (2024) 103.** Lattice logs as a sum of similar connectivity events at many scales. Four-point of the density field. Support for an LCFT description; no spin-4 number.

---

## 2. Q3 — discrete holomorphy on `ℤ[i]` is C4

**Cardy, J. Stat. Phys. (2009),** [0907.4070](https://arxiv.org/abs/0907.4070). Discretely holomorphic lattice observables satisfy a Cauchy–Riemann relation on neighbouring faces. They have **fractional spin by construction**, exist only on the integrable critical manifold, and are the lattice precursors of CFT parafermions.

On the square grid the rotation that enters the discrete CR equation is multiplication by `i` — the same automorphism Q3 proves pins every C3 character to a real line. The literature's local spin on this lattice is C4 (or a parafermion whose spin is a multiple of 1/4 in the Ising/Potts cases), not C3. Q3's no-go is the percolation-site version of that fact. A published escape would have to leave the Gaussian-ideal family; Cardy does not provide one.

---

## 3. Matching neighbourhoods are MC, not algebraic

**Malarz, [2303.10423](https://arxiv.org/abs/2303.10423)** (v4 2025): square site, complex neighbourhoods through the 6th coordination zone. Newman–Ziff + FSS. Power law `p_c(ζ) ∝ ζ^{-γ_2}` with `γ_2 = 0.5454(60)`, `ζ = Σ z_i r_i`. **Ciepłucha et al., [2503.16703](https://arxiv.org/abs/2503.16703):** through the 7th zone, 64 thresholds, range `0.27013` down to `0.11535`.

NN+NNN (Sq8, matching of square) sits in this list. The papers do not claim exact formulae. They are orthogonal to Sykes–Essam / Grimmett–Li (which constrain the *pair*) and to P2 (which excludes low-degree identities for *square* site). Do not enlarge the exact-threshold table with these roots.

---

## 4. Community, and a non-source

**Diskin, Easo, Radhakrishnan, Sudakov, Tassion, [2603.03257](https://arxiv.org/abs/2603.03257)** (3 Mar 2026). For every infinite transitive graph, `p > p_c`:

```text
P_p(n ≤ |C_o| < ∞) ≤ exp(−c Φ(n))
P_p(o ↔ ∂B_n but o ̸↔ ∞) ≤ exp(−c n)
```

`Φ` is the isoperimetric function. On `ℤ^d` this is Aizenman–Barsky / Menshikov. The theorem does not produce a square-site `p_c`, a wrapping formula, or a matching-odd exponent.

Quanta, 31 Aug 2026, popular account. X 4 Sep 2026: Quanta post. A 1 Sep post mentions an Alpoge rumour; ignore until there is a paper.

---

## 5. Opinions for subsequent analysis

1. **Q1 / Astra.** Ask the external model the level-4 question as written. Do not substitute He's `γ = -5/4`. Cite He 2025 as computing energy–hull only.
2. **Q1 fit.** Do not fit matching-odd to `−5/4` or `a_0 = -25/48`. Colour and spin first (q1-noise note).
3. **Q3.** Treat as literature-consistent. Next escape is a lattice that is not a Gaussian ideal quotient, not a new CFT operator.
4. **P2 / Q4 / #566.** Unchanged. Neighbourhood MC papers are not exact thresholds.
5. **#567.** Unchanged (Pinson/Arguin or a non-claim). Diskin–Tassion is not a wrapping competitor.
6. **Do not** open a compute block on supercritical sharpness; `ℤ²` is settled.

## Not established

- a published coupling for `(x=21/4, s=4)`;
- that high-rank Jordan blocks exist in percolation (He: “suggests”);
- an exact NN+NNN site threshold other than matching against square;
- anything in the claim ledger.
