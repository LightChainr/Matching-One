# Literature officer — Q1 / #227 / colour, 2026-09-05

**Role:** theory input. Companion to `notes/literature-officer-brief-20260905.md` (#571) and `notes/literature-officer-20260905-deeper.md` (this PR). **Does not enter** `docs/STATUS.md`.

GitHub code search of this repository on 2026-09-05 for `Garban`, `noise sensitivity`, `Tan Couvreur` returned **zero** hits. Those papers are the lattice language for Q1's non-scalar log pair and for #227.

---

## Ticket map

| Ticket | What this pass changes |
|---|---|
| Q1 / `x=21/4` spin-4 | Tan–Couvreur–Deng–Jacobsen 2019: `Δ_{4s} = 21/4` is the **scalar** (`spin 0`) 4-cluster operator. Matching-odd is orientation-sensitive, so it is not that operator. The N=4 table has spins 0, 1, 2, 3, 6 — **not 4**. |
| Q1 Jordan / `μ = -5/4` | Same paper: energy mixes with the 2-cluster operator; the combination that cancels the power is a pure log with universal slope `2√3/π`. That is a lattice log observable with a number. It is not the level-4 pairing. |
| colour | Nolin–Qian–Sun–Zhuang PRL **134**, 117101 (2025): exact formula for *same-colour* two-path annulus crossing; other roots of the backbone equation are proposed as a CFT bulk spectrum. Monochromatic `k>2` remains open. If matching-odd is same-colour wrapping, `21/4` is the wrong formula. |
| #227 | Garban–Pete–Schramm, Acta Math. **205** (2010) 19: the Fourier spectral sample of a crossing. Pivotals are four-arm (`x=5/4`). Noise sensitivity iff `ε E[|pivotals|] → ∞`. Square bond has exceptional times a.s. This is the theorem #227 is trying to turn into an H4 tomography. |

X this round: no percolation-research posts (rugby wrapping). Community read unchanged.

---

## 1. Q1 — `21/4` in the literature is spin 0

**Tan, Couvreur, Deng, Jacobsen, PRE 99, 050103(R) (2019)** = [arXiv:1809.06650](https://arxiv.org/abs/1809.06650). Critical **bond** percolation, N-cluster operators projected onto irreps of `S_N`. In d=2 the dimensions are Kac, and the conformal spin is the Young-tableau index.

N=4 (four distinct clusters = 8-leg watermelon), Table 1:

| irrep | Δ | spin |
|---|---:|---:|
| `4s` `[4]` | **21/4** | **0** |
| `4m1` `[3,1]` | 339/64 | 1 |
| `4m2` `[2,2]` | 87/16 | 2 |
| `4m3` `[2,1,1]` | 363/64 | 3 |
| `4a` `[1⁴]` | 111/16 | 6 |

N=2: `Δ_{2s} = 5/4` (energy / four-arm, spin 0), `Δ_{2a} = 23/16` (spin 1).

The repository candidate is a matching-odd, orientation-sensitive amplitude with `x=21/4` and spin 4. The literature operator with `x=21/4` is scalar. A spin-4 descendant of the energy (`5/4`) would sit at a different `(h, \bar h)` than `Δ_{4s}`. **Do not identify matching-odd with `P_{4s}`.** Colour-decompose, then check the spin: a `cos 4θ` lattice amplitude is spin 4, which is absent from this N=4 table.

Logarithmic combination, N=2, exact in d=2:

```text
F(r) = [P_0(r) + P_1(r) − (P_≠)²] / P_{2s}(r)  ∼  δ log r
δ = 2√3 / π  ≈ 1.10266
```

MC: `1.15 ± 0.05`. This is Vasseur–Jacobsen–Saleur's lattice log, with a measured slope. Q1's `μ = -5/4` (He 2024) is the *primary-level* indecomposability parameter; `δ` is a two-point slope. They are related in LCFT and are not interchangeable as numbers to fit.

Simulations are square-bond and triangular-bond, up to L=8192, N≤4. Not matching-odd, not Sq8.

---

## 2. Colour — same-colour two-path has a formula; 8-arm mono does not

**Nolin, Qian, Sun, Zhuang, PRL 134, 117101 (2025)** = [arXiv:2410.06419](https://arxiv.org/abs/2410.06419). The first brief cited the 2023 derivation ([2309.05050](https://arxiv.org/abs/2309.05050)) for transcendence of the backbone. The PRL adds the full same-colour two-path annulus probability. Backbone `x_B` is the unique root in `(1/4, 2/3)` of

```text
√(36x + 3)/4  +  sin(2π √(12x+1) / 3)  = 0
```

`x_B = 0.3566668…`, transcendental. The other roots of a companion equation

```text
S = { s ∈ ℂ : sin(4π √(s/3)) + (3/2)√s = 0 }  \  {0, 1/3}
```

are proposed as a CFT bulk spectrum (numerical `{0.440, 2.194±0.601i, …}`). **Monochromatic k>2 remains open.** There is still no closed form for a same-colour 8-arm, which is what a matching-odd wrapping would be if it is monochromatic.

If the N=580 slope is a wrapping of one colour, the continuum prediction is this transcendental family, not `21/4`. That is the colour test the first brief asked for, now with a formula on the mono-2 side.

---

## 3. #227 — the spectral sample of a crossing is four-arm physics

**Garban, Pete, Schramm, Acta Math. 205 (2010) 19–104,** [arXiv:0803.3750](https://arxiv.org/abs/0803.3750). Fourier–Walsh expansion of a planar crossing indicator `f`. The spectral sample `S_f` (a random subset of bits with `P[S=S] = f̂(S)²`) is the object. Sharp lower-tail bounds. Consequences:

- noise sensitivity of `f` iff `ε E[|pivotals|] → ∞`;
- exceptional times of dynamical critical percolation exist a.s. on the **square grid** (bond) as well as the triangle;
- dimension of exceptional times on the triangle: `31/36` (plane), `5/9` (half-plane).

The pivotal density is the four-arm event, `α_4(R) = R^{-5/4+o(1)}` (Smirnov–Werner). So the Fourier spectrum of a wrapping/crossing is controlled by `x=5/4`, not by `x=21/4`.

**For #227.** A noise-operator tomography of the matching H4 sector is a spectral sample of a *matching-odd* Boolean function, not of Cardy crossing. GPS tells you the spectrum of the even wrapping. The odd sector is a different `f`. The cheap control is: apply the same noise semigroup to the even wrapping and recover GPS (four-arm / `5/4`) before claiming an H4 readout on the odd channel. Square bond already has the theorem; square site is believed the same.

Garban–Steif, *Noise Sensitivity of Boolean Functions and Percolation*, CUP 2014, is the book form. Not cited anywhere in this repository.

---

## 4. Opinions for subsequent analysis

1. **Q1.** Treat `21/4` as the scalar 4-cluster dimension until a spin measurement says otherwise. Matching-odd's `cos 4θ` is a spin-4 claim and needs a different `(h, \bar h)` than `P_{4s}`.
2. **Q1 log.** Fit `δ = 2√3/π` only to the Vasseur two-point combination, not to the matching-odd slope. He 2024 still owns `μ = -5/4`.
3. **Colour, before N=580 or any `x=21/4` fit.** Same-colour 2-path has a PRL formula; same-colour 8-path does not. Decompose matching-odd by colour.
4. **#227.** GPS first on the even wrapping (known theorem), then noise on the odd channel. Do not skip the even control.
5. **#567.** Unchanged from the deeper note: Pinson/Arguin or an explicit non-claim, and keep the two `11/4`s apart.
6. **P2 / Q4 / #566.** Unchanged. No new exact-threshold polynomial this round.

## Not established

- identification of matching-odd with any Tan `P_{N∘}`;
- that a spin-4 lattice amplitude must appear in the `S_N` table;
- a monochromatic 8-arm exponent;
- anything in the claim ledger.
