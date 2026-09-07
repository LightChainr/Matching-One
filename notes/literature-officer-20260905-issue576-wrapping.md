# #576 — wrapping on published ground, 2026-09-05

Plumbing and theory input for [issue #576](https://github.com/LightChainr/Matching-One/issues/576). **Does not enter** `docs/STATUS.md`. No claim-ledger row.

The repository already had `π({1,0})(τ=i)` to 50 digits (`notes/pinson-arguin-primitive-baseline.md`, `predictions/p156_pinson_arguin_baselines_20260829.json`). It did **not** have the rectangular values the N=580 ladder needs, and it had never compared a wrapping count to Akhunzhanov’s exact polynomials.

---

## Part 2 — Pinson `π({1,0})(ir)` at `r = 1, 2, 4`

**Formula used** (Pruessner–Moloney [cond-mat/0310361](https://arxiv.org/abs/cond-mat/0310361), eqs. (1)–(2); Pinson 1994; same numbers as Arguin `hep-th/0111193` Table 1). `g = 2/3` in *this* normalization — not Arguin’s Coulomb-gas `g`.

```text
Z_{m,n}(g; r) = √g / (√r η(e^{-2π r})²) * exp[-π g (m²/r + n² r)]

P̂((a,b), ≥1, r)
  = Σ_l Z_{a 3l,     b 3l}
    − ½ Σ_l Z_{a(3l+1), b(3l+1)}
    − ½ Σ_l Z_{a(3l+2), b(3l+2)}
    −   Σ_l Z_{a 2l,     b 2l}
    +   Σ_l Z_{a(2l+1), b(2l+1)}
```

`P̂((1,0), ≥1, r)` is Arguin’s `π_1({1,0})` (rank-1 homology `{1,0}`, **not** Newman–Ziff `R_h`). At `τ = i` it reproduces the frozen baseline to all 12 digits printed here:

```text
0.169415435321  =  predictions/…json  0.169415435321346889…
```

| `r` | `π({1,0})(ir)` | `π({0,1})(ir)` | ratio to `r=1` |
|---:|---:|---:|---:|
| 1 | 0.169415435321 | 0.169415435321 | 1 |
| 2 | 0.503035897695 | 0.024949480334 | **2.969244784222** |
| 4 | 0.855969321054 | 0.000762380303 | **5.052487215408** |

**Cross-check, Pruessner eqs. (15)–(16).** Their truncated large-`r` expansion versus the sum above:

| `r` | relative error of eq. (15) |
|---:|---:|
| 1 | `4.5 × 10^{-4}` |
| 2 | `6.9 × 10^{-9}` |
| 4 | `< 10^{-15}` |

The ticket asked for relative `< 10^{-8}` at `r = 2`. That bound is on the *expansion versus the exact sum*, and it holds. Arguin Table 1 “analytical” entries (`0.16948`, `0.50293`) are five-digit truncations of the same CFT evaluation; they differ from the sum at `~10^{-4}` and are not a second source.

**What this is not.** Newman–Ziff wrapping in a *specified direction*, including simultaneous wrap in the other, is `R_h(i) = 0.521058290…`. That is `π({1,0}) + π(ℤ×ℤ) + diagonal wrapping that covers x`, and it is the limit of the Akhunzhanov polynomials (Part 1). Do not score `R_h` against the table above.

### Written choice (Part 2 of the ticket)

**Matching-odd slope is not claimed to be a Pinson wrapping.** That is a non-claim, not a postponement. `P4[S']` is orientation-sensitive; Pinson `π({1,0})` is a scalar homology class. Nothing in the N=290 fingerprint, and nothing in this formula, identifies them.

**The three numbers do enter the next frozen ladder as a named competitor** `pinson_pi10_ratio`, because that is the trap #572 flagged:

| competitor | `r=2` | `r=4` |
|---|---:|---:|
| area-normalized weight-4 `Ê4(ri)/Ê4(i)` | 2.75 | 10.99 |
| **Pinson `π({1,0})(ri)/π({1,0})(i)`** | **2.969** | **5.052** |
| bare aspect `r` | 2 | 4 |
| plain area `r²` | 4 | 16 |
| no dependence | 1 | 1 |

At `r=2` Pinson sits **8% from `11/4`**. Scoring against `E4` and reading `r` off the leftover is exactly how N=290 acquired a post-hoc `bare_aspect_ratio`. N=290’s measured `1.880 ± 0.177` is itself `6σ` from `2.969`, so Pinson does not *explain* that run; it still has to be named **before** the next one.

The three `11/4`s, kept apart:

| where | what |
|---|---|
| `Ê4(2i)/Ê4(i) = 11/4` | modular weight-4 amplitude ratio |
| Newman–Ziff estimator `L^{-11/4}` | `1/ν + θ = 3/4 + 2`, FSS rate |
| Pinson ratio at `r=2` | `2.969`, a different number that happens to sit nearby |

Mertens–Ziff matching-function root `~ L^{-4}` is a fourth exponent (see the mertens note). None of these predicts that a spin-4 amplitude ratio equals `11/4`.

---

## Part 1 — Akhunzhanov torus polynomials

**Source.** R. K. Akhunzhanov, A. V. Eserkepov, Y. Y. Tarasevich, *J. Phys. A* **55**, 204004 (2022), [arXiv:2204.01517](https://arxiv.org/abs/2204.01517). Ancillary `anc/torus.txt` in the arXiv source tarball. Square **site**, wrapping **along one specified direction** on the `L × L` torus. **Not** Sq8 / NN+NNN.

Coefficients `c_k` are absolute counts: `R_L(p) = Σ_k c_k p^k (1-p)^{L²-k}`. Last coefficient is `1` (full occupancy wraps). First nonzero is `c_L = L` (the `L` full rows in the specified direction).

**Independent exact enumeration, L = 3 and L = 4.** Union-find with covering-space winding, `2^{9}` and `2^{16}` configs, wrapping in `x` including simultaneous wrap in `y`. Coefficients match `torus.txt` **exactly**:

```text
L=3:  [0,0,0, 3, 18, 45, 63, 36, 9, 1]
L=4:  [0,0,0,0, 4, 48, 280, 1008, 2558, 4480, 5088, 3664, 1744, 560, 120, 16, 1]
```

So the published polynomial is Newman–Ziff `R_h`, not Arguin `π({1,0})`.

**`R_L(p_c)` at `p_c = 0.592746`, published coefficients.** Limit should be Pinson/Newman–Ziff `R_h = 0.521058290…`.

| `L` | `R_L(p_c)` | note |
|---:|---:|---|
| 3 | 0.5212736273 | exact enum match |
| 4 | 0.5173034301 | exact enum match |
| 5 | 0.5171953205 | |
| 6 | 0.5178345949 | |
| 7 | 0.5184255041 | |
| 8 | 0.5188927030 | |
| 9 | 0.5192535596 | |
| **10** | **1.28 (impossible)** | **ancillary file is corrupt: 100 coeffs, not 101; `c_N = 100` not `1`. Do not use.** |
| 11 | 0.5197544660 | |
| 12 | 0.5199306392 | approaching `0.521058` |

**Engine.** The committed transfer wrapping channel was **not** scored against these polynomials in this pass (no engine in this session). What is now on published ground is the *observable*: specified-direction wrap, primary square site, `L ≤ 12` except `L = 10`. Comparing the engine at `L = 3` or `L = 4` against the two lines above is a one-shot plumbing check — if it disagrees, the channel is wrong; if it agrees, later Pinson comparisons at production sizes are a checked channel against a continuum formula.

The paper does not treat NN+NNN, so this does not check the matching half of Mertens–Ziff `M_L`.

---

## What to do on the next freeze (#567 / N=580 or N=3380)

1. Name `pinson_pi10_ratio` as a competitor with the two numbers `2.969244784222` (`r=2`) and `5.052487215408` (`r=4`). Prospective, not post-hoc.
2. Write the matching-odd non-claim in the same freeze file.
3. Keep `R_h` (Akhunzhanov / Newman–Ziff, → `0.521`) and `π({1,0})` (→ `0.169`) as different observables.
4. Engine vs `L=3` polynomial, when someone next touches the wrapping binary.

## Not established

- that the transfer engine agrees with Akhunzhanov (not run);
- any identification of matching-odd with Pinson wrapping (explicitly a non-claim);
- a published table for Sq8 wrapping (still none);
- anything in the claim ledger.
