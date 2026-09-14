# Topological pivotal-pair atlas — exact raw joint table for issue #769

**Author:** rev769 (cloud VM `DevEnvC_551oUR`, `/workspace/rev769/`)
**Date:** 2026-09-14
**Scope:** issue #769 Phase A (exact L=5 main product) + Phase B (D4 displacement-orbit
structure). Phase C is **not** completed as a classifier; only raw material is saved.
**Honesty flags** are in §8. Everything below is *computed by me* unless marked
"引用 (quoted, not independently recomputed)". The T1 control values are **not**
引用 — they are independently reproduced.

---

## §T0 What is derivable and what must be recomputed (with reasons)

This list is itself part of the delivery (task spec §0).

### T0-A. Derivable **without** any new enumeration

| # | quantity | why it is derivable |
|---|---|---|
| A1 | `M(p) = E_p[X]`, and **all** its `p`-derivatives `M′(p), M″(p), …` | `X = r_black − 1` depends on `ω` only through the rank `r_black`. Hence the *entire* `p`-dependence of `M` is a functional of the rank marginal `C[L,j,k] = #{r=j, |ω|=k}` alone. Computed **exactly** (rational) from the reused #775 tables `C_L3..5.json`. |
| A2 | the two Russo/translation constraints `N·Σ_{d≠0}J_d = M″(p)` and `N·E[Δ₀X] = M′(p)` | these are identities of the object definition, with the left sides from the pair atlas and the right sides from A1. Used as the **independent** side of the T2 check. |
| A3 | D4 displacement orbits: representative, multiplicity, torus distance | pure lattice combinatorics, independent of any rank computation. |
| A4 | `Σ_orbit mult = N − 1` | arithmetic check on A3 (2/5/5 orbits at L=3/4/5). |
| A5 | `Δ_vX ∈ {0,1,2}` and hence `Δ₀Δ_dX ∈ {−2,…,+2}` | a priori: `r_black` is monotone under adding black sites, so a single-site difference cannot be negative and cannot exceed 2. This is why **no** value outside `{−2..2}` occurs (§T4) and why the sign convention `Δ_vX = X(v=1) − X(v=0)` is *forced*, not chosen. |

### T0-B. **Not** derivable → must be produced by extending the same enumeration

| # | quantity | why the marginal `C[L,j,k]` cannot give it |
|---|---|---|
| B1 | histogram of `Δ₀X` by occupancy of the other `N−1` sites (the counts `{0:161,1:90,2:5}` etc.) | `C[L,j,k]` is a **marginal** histogram; `Δ₀X = X(0=1) − X(0=0)` is a **differential** observable of the pair (site 0, rest). The only functional of it fixed by the marginal is its **mean** (via Russo, = A2). The *distribution* is not fixed. |
| B2 | for every displacement `d`: histogram of `Δ₀Δ_dX` over `{−2..2}` by occupancy of the remaining `N−2` sites | two-site **mixed** discrete difference — no functional of the rank marginal determines it. The sum rule A2 fixes only `Σ_d J_d`, i.e. one number out of the whole orbit table, and only for the *signed* moment (not for `A_d = E|Δ₀Δ_dX|`). |
| B3 | `J_d, A_d, P⁺_d, P⁻_d` | linear functionals of the B2 histograms. |

### T0-C. Reuse, not re-production

The expensive component — the rank / winding kernel — is **copied verbatim** from
`rank775`'s delivery (issue #775):

* path **A** (`pivotal.c`) reuses `rank775-work/ranksec_brute.c::rank_of`
  (union-find with `Z²` offsets, forward edges only);
* path **B** (`pivotal_bis.c`) reuses `rank775-work/rank_bis.c::rank_uf`
  (generic-adjacency union-find) and `rank_bis.c::rank_cov`
  (universal-cover `Z²` BFS, both edge orientations).

Only the **dimension** was added (site-0 state × displacement orbit) plus the
occupancy-count histogram accumulation. No second rank production was written.

### T0-D. Explicitly not done (constraints of #769)

no reconstruction of `A_d` from the one-dimensional `M(p)` polynomial (forbidden);
no Monte Carlo; no L≥6; no STATUS / #275 contract change; no arm-count classifier
(Phase C only saved raw material, §T5); no assertion of an exponent anywhere.

---

## §T1 Control-value reproduction — **all 12 numbers reproduced**

Frozen reference value `p_c = 0.59274605079210` (as given by #769). All mine are
**exact rational** evaluations (`Fraction`, `p_c` treated as the exact decimal
`59274605079210/10^14`); the table prints 20 significant digits.

### T1 control table

| L | quantity | #769 quoted | mine (20 digits) | delta |
|---|---|---|---|---|
| 3 | counts | `{0:161, 1:90, 2:5}` | `{0: 161, 1: 90, 2: 5}` | **0 (exact)** |
| 3 | Ed0 | 0.4458689104616698 | 0.44586891046166977571 | -5.551e-17 |
| 3 | jump2 | 0.1039113209 | 0.10391132091753346805 | +1.753e-11 |
| 3 | A2 | 23.04926674699 | 23.049266746993427819 | +3.428e-12 |
| 3 | Mpp | 2.336146088861 | 2.3361460888611493896 | +1.492e-13 |
| 3 | cancel | 9.86636 | 9.8663636049532082272 | +3.605e-06 |
| 4 | counts | `{0:24639, 1:7840, 2:289}` | `{0: 24639, 1: 7840, 2: 289}` | **0 (exact)** |
| 4 | Ed0 | 0.31159887509994 | 0.31159887509993999399 | +0.000e+00 |
| 4 | jump2 | 0.07434091137 | 0.074340911371403713232 | +1.404e-12 |
| 4 | A2 | 45.48519754007 | 45.485197540074412686 | +4.412e-12 |
| 4 | Mpp | 2.933170180195 | 2.9331701801949563475 | -4.352e-14 |
| 4 | cancel | 15.50718 | 15.507179858568993614 | -1.414e-07 |

Legend: `Ed0 = E[Δ₀X]`; `jump2` = rank-jump-2 share of `M′`; `A2 = A2_abs`;
`Mpp = M″`; `cancel = A2_abs/|M″|`.

Every deviation is **at or below the precision at which #769 quoted the value**
(e.g. the `cancel` deviations are ≤ 3.6e-6 while the quoted values have 6 digits).
The `counts` are equal **exactly**, as integers.

### Convention judgement (T1 asked me to decide this myself)

* the counts sum to `2^{L²−1}` (256 / 32768): they are counts over the configurations
  of the **other `L²−1` sites**, i.e. the distribution of `Δ₀X`;
* `Δ₀X = X(v=1) − X(v=0)` **does not depend** on the state of site 0 itself, so there
  is no "site 0 fixed black vs white" ambiguity for the counts. The state of site 0
  enters only *inside* each term of the difference, which is exactly how I evaluate it;
* the **direction** is forced and not a free choice: `r_black` is monotone under adding
  black sites, hence `Δ₀X ≥ 0` (my enumeration produced **zero** negative values, and
  the L=3 list is `{0:161, 1:90, 2:5}` with no 3-value);
* the decisive evidence for this reading is the simultaneous exact match of
  **3 integer counts + 5 derived quantities at two sizes** with two independent code
  paths. Nothing was substituted.

---

## §T2 The two required identities — **both hold exactly**

| L | `N·E[Δ₀X]` (atlas) | `M′(p)` (exact from reused `C[L,j,k]`) | equal? | `N·Σ_{d≠0}J_d` (atlas) | `M″(p)` (exact from reused `C[L,j,k]`) | equal? |
|---|---|---|---|---|---|---|
| 3 | 4.0128201941550279814 | 4.0128201941550279814 | **yes** | 2.3361460888611493896 | 2.3361460888611493896 | **yes** |
| 4 | 4.9855820015990399038 | 4.9855820015990399038 | **yes** | 2.9331701801949563475 | 2.9331701801949563475 | **yes** |
| 5 | 5.8827716932504600604 | 5.8827716932504600604 | **yes** | 3.4161276168223449860 | 3.4161276168223449860 | **yes** |

The comparison is performed as an **exact `Fraction` equality** (`Mpp_enum − ddEX == 0`
and `Mprime_enum − dEX == 0`), not as a float tolerance. The right column is obtained
by differentiating `M(p) = Σ_k G[k] p^k q^{N−k}` symbolically, with
`G[k] = Σ_j (j−1) C[L,j,k]` from the **reused #775 rank tables** — i.e. it uses only
the rank marginal and **none** of my pair data. Both identities therefore hold to
exact arithmetic, at all three sizes.

`M″` at L=3 equals the #769 control value `2.336146088861` (see §T1): the identities
and the controls are satisfied simultaneously by the same numbers.

> Note for readers: these identities are consequences of Russo's formula plus
> translation invariance, so their holding is a **consistency/validation** result
> (it certifies that the pair atlas is the correct joint object), not an independent
> piece of physics evidence. The non-trivial content is that the *independently
> enumerated* pair table reproduces the second derivative of the rank-marginal
> polynomial exactly.

---

## §T3 Phase A — exact L=5 (main product)

### T3 identity / L=5 table

| L | counts (j=0,1,2) | sum | E[Δ₀X] | jump2 share | J_sum | M′ | M″ | A2_abs | cancel_ratio |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 161,90,5 | 256 (=2^8) | 0.44586891046166977571 | 0.10391132091753346805 | 0.25957178765123882107 | 4.0128201941550279814 | 2.3361460888611493896 | 23.049266746993427819 | 9.8663636049532082272 |
| 4 | 24639,7840,289 | 32768 (=2^15) | 0.31159887509993999399 | 0.074340911371403713232 | 0.18332313626218477172 | 4.9855820015990399038 | 2.9331701801949563475 | 45.485197540074412686 | 15.507179858568993614 |
| 5 | 13800255,2906128,70833 | 16777216 (=2^24) | 0.23531086773001840242 | 0.052085344508862943781 | 0.13664510467289379944 | 5.8827716932504600604 | 3.4161276168223449860 | 70.908531625515400319 | 20.756991418099877127 |

**L=5 facts**

* full enumeration of `2^25 = 33 554 432` configurations; the base enumeration is
  `2^24` configurations of the sites `≠ 0` (33.5 M rank evaluations, 0.88 s wall,
  16 threads). Per-`d` pair totals are each exactly `2^23 = 8 388 608`, and the
  single-site total is exactly `2^24` — the enumeration is complete and exact;
* `Δ₀Δ_dX` support at L=5 is the **full** `{−2,−1,0,1,2}` (at L=3 only `{−1,0,1}`;
  L=4 already reaches the full range);
* L=6+ was **not** attempted (not an acceptance condition for #769).

### T3 cross-validation — **three code paths, byte-identical at L=3, 4 and 5**

| path | rank kernel | algorithm |
|---|---|---|
| A | union-find + `Z²` offsets (`ranksec_brute.c::rank_of`) | build `j0[γ]` once (array), then table-lookup pairs |
| B-uf2 | generic-adjacency union-find (`rank_bis.c::rank_uf`) | recompute **all four** paired-flip ranks `X(0=a,d=b)` directly |
| B-cov | universal-cover `Z²` BFS (`rank_bis.c::rank_cov`) | same direct four-rank recomputation, structurally different kernel |

Result: the **single-site histogram and all per-`d` pair histograms are identical as
integers** for A vs B-uf2 and A vs B-cov, at **L=3, L=4 and L=5** (not only L=3,4).
So the L=5 numbers satisfy the #769 requirement "两个实现路径至少在 L=3,4 完全一致后才信 L=5"
with margin.

---

## §T4 Phase B — D4 displacement-orbit structure

`p_c = 0.59274605079210`. Orbit-level sums; `A_d`, `J_d`, `P±_d` are **per pair**
(i.e. divided by the orbit multiplicity). `hist identical = yes` means every member of
the orbit produced an identical exact integer histogram (internal symmetry check, all
**yes**).

**L3** (N=9, P_pivotal=0.42270349674060299395, P_piv²=0.17867824615673296593)

| rep d | mult | torus dist | J_d·mult | A_d·mult | J_d | A_d | P⁺_d | P⁻_d | A_d/P_piv² | hist identical |
|---|---|---|---|---|---|---|---|---|---|---|
| (0,1) | 4 | 1.0000 | 0.967840 | 1.652645 | 0.241960 | 0.413161 | 0.327561 | 0.085601 | 2.312320 | yes |
| (1,1) | 4 | 1.4142 | -0.708268 | 0.908384 | -0.177067 | 0.227096 | 0.025015 | 0.202082 | 1.270977 | yes |

**L4** (N=16, P_pivotal=0.30001660292132312676, P_piv²=0.090009962028450872520)

| rep d | mult | torus dist | J_d·mult | A_d·mult | J_d | A_d | P⁺_d | P⁻_d | A_d/P_piv² | hist identical |
|---|---|---|---|---|---|---|---|---|---|---|
| (0,1) | 4 | 1.0000 | 0.528192 | 1.130011 | 0.132048 | 0.282503 | 0.206411 | 0.075227 | 3.138573 | yes |
| (0,2) | 2 | 2.0000 | 0.191503 | 0.432418 | 0.095751 | 0.216209 | 0.155789 | 0.060229 | 2.402057 | yes |
| (1,1) | 4 | 1.4142 | -0.452616 | 0.837505 | -0.113154 | 0.209376 | 0.047877 | 0.160971 | 2.326146 | yes |
| (1,2) | 4 | 2.2361 | -0.043599 | 0.350775 | -0.010900 | 0.087694 | 0.038397 | 0.049297 | 0.974267 | yes |
| (2,2) | 1 | 2.8284 | -0.040157 | 0.092115 | -0.040157 | 0.092115 | 0.025979 | 0.061965 | 1.023385 | yes |

**L5** (N=25, P_pivotal=0.22918274392381965818, P_piv²=0.052524730112451096451)

| rep d | mult | torus dist | J_d·mult | A_d·mult | J_d | A_d | P⁺_d | P⁻_d | A_d/P_piv² | hist identical |
|---|---|---|---|---|---|---|---|---|---|---|
| (0,1) | 4 | 1.0000 | 0.384887 | 0.842637 | 0.096222 | 0.210659 | 0.153181 | 0.057219 | 4.010669 | yes |
| (0,2) | 4 | 2.0000 | 0.166644 | 0.514443 | 0.041661 | 0.128611 | 0.084998 | 0.043475 | 2.448575 | yes |
| (1,1) | 4 | 1.4142 | -0.343235 | 0.697267 | -0.085809 | 0.174317 | 0.044125 | 0.129912 | 3.318754 | yes |
| (1,2) | 8 | 2.2361 | -0.016390 | 0.561085 | -0.002049 | 0.070136 | 0.033962 | 0.036061 | 1.335287 | yes |
| (2,2) | 4 | 2.8284 | -0.055261 | 0.220909 | -0.013815 | 0.055227 | 0.020688 | 0.033977 | 1.051453 | yes |

The complete exact integer histograms (all `d`, all `Δ∈{−2..2}`, all occupancy bins
`k = 0..N−2`) are in `results/topological-pivotal-pair/orbit-table-20260914.json`;
they are `p`-independent, so any other reference `p` can be used without re-enumerating.

### T4 readings

1. **Signed interaction changes sign with distance and does not decay monotonically.**
   At L=5 the per-pair signs run `+ (d=1), − (√2), + (2), − (√5), − (2√2)`; the
   diagonal pair is strongly *negative* and the NN pair strongly *positive*. The
   signed magnitudes at `√5` and `2√2` are ~2 orders of magnitude below the NN one.
   ⇒ "fast sign change": `compatible`, but not a clean one-sign decay.
2. **Absolute interaction is concentrated at short distance.**
   per-pair `A_d` at L=5: `0.2107 (1) > 0.1743 (√2) > 0.1286 (2) > 0.0701 (√5) > 0.0552 (2√2)`.
   At L=4 the ordering is *not* monotone in distance (`0.2162` at dist 2 vs `0.2094` at
   dist √2) — a clear finite-size/torus artefact and a warning against reading this as
   a radial profile.
3. **Near/far comparison against "two independent pivotals".**
   `A_d / P_piv²` is **≈ 1.0 at the largest available separation** (L=3: 1.27 at √2,
   which is the L=3 diameter; L=4: 1.02 at 2√2; L=5: 1.05 at 2√2), i.e. the far field is
   at the scale of the product of two single-point pivotal probabilities, while the
   NN pair is **2.31 / 3.14 / 4.01 ×** that product at L=3 / 4 / 5. So there *is* a
   short-distance excess over the independent-pivot scale, and its size grows with L.
4. **Against the named fusion scales — `unresolved`.**
   `A_d` is a bounded lattice diagnostic (`≤ 2`), and the only scale I can compare it to
   honestly is the single-point pivotal product in 3. With 2/5/5 orbits at three sizes
   and a torus diameter of ≤ 2√2, I **cannot** distinguish "single-point pivotal
   probability squared" from "6-arm fusion" from "8-arm fusion" magnitudes:
   ⇒ write `unresolved`. The near/far **pattern** is `compatible`; the **fusion
   identification is unresolved** (it would need the arm witness classifier, §T5).

### Maximal / representative pair comparisons (as requested)

* nearest-neighbour `(0,1)` vs diagonal `(1,1)`: opposite sign of `J`, both with large
  `A`; this sign opposition is the main source of the cancellation in §D1.
* intermediate `(0,2)` (dist 2) and `(1,2)` (dist √5): small signed `J`, moderate `A`
  (the L=5 `(1,2)` orbit has multiplicity 8, the largest).
* maximally separated `(2,2)` (dist `2√2 = 2.828`, the L=5 diameter): smallest `A`,
  negative `J`, and `A_d/P_piv² ≈ 1.05`.
* the `L=3` torus has only **two** orbits (the diameter is √2), so L=3 cannot separate
  "nearest" from "farthest" at all — a hard finite-size limitation.

---

## §T5 Phase C (bonus) — **not** completed as a classifier; raw material saved

The issue explicitly allows, when an exact arm classifier is too hard, saving
representative configurations + lifted homology. That is what was done
(`witness_L3/4/5.json`; 90 witnesses at L=5, 5 per (orbit, `Δ`-value) with
`Δ ≠ 0`). Each witness records, for the **four** paired-flip states
`(s₀,s_d) ∈ {0,1}²`: the rank `r_black`, the number of black components, the explicit
list of black winding vectors in `Z²`, plus the black-site lists of the two extreme
states. Example worth reading (L=5, orbit `(0,1)`, `Δ=+1`): rank is `0` in all three
states except `(s₀=1,s_d=1)`, where a single `+x` winding appears — i.e. the two sites
must be occupied *together* to close a winding cycle.

**No arm count is claimed anywhere.** Local degree counts would be exactly the
"guessing from local degree" that #769 forbids.

---

## §D1–D3 The three required verdicts

### D1 — Is the smallness of `M″` caused by increasingly severe ± pivotal-pair cancellation?

**Verdict: `compatible` (on the three computed sizes).**

| L | A2_abs | abs(M″) | cancel_ratio | A2_abs/N² | rank-jump-2 share |
|---|---|---|---|---|---|---|
| 3 | 23.049266746993427819 | 2.3361460888611493896 | 9.8663636049532082272 | 0.284559 | 0.10391132091753346805 |
| 4 | 45.485197540074412686 | 2.9331701801949563475 | 15.507179858568993614 | 0.177677 | 0.074340911371403713232 |
| 5 | 70.908531625515400319 | 3.4161276168223449860 | 20.756991418099877127 | 0.113454 | 0.052085344508862943781 |

`cancel_ratio` increases monotonically `9.87 → 15.51 → 20.76`: the signed total `M″`
grows only `2.34 → 2.93 → 3.42` while the absolute pair activity grows `23.05 → 45.49
→ 70.91`. So the *surviving fraction* `|M″|/A2_abs` falls `0.1014 → 0.0645 → 0.0482`.
Independently, the share of `M′` carried by rank-jump-2 events falls
`0.1039 → 0.0743 → 0.0521` (as #769's Round-5 note reported for L=3,4 — 引用 for the
trend direction, but my values are independently reproduced).
What I **cannot** say: that `cancel_ratio` follows a power law; three finite sizes
cannot certify an exponent.

### D2 — Does the absolute pair interaction show a clear near "fusion" region and a far "two-pivotal" region?

**Verdict: near/far scale separation `compatible`; fusion identification `unresolved`.**

Supporting facts are in §T4: per-pair `A_d` decreases with torus distance at L=5
(0.211 → 0.055), and `A_d/P_piv²` goes from `4.01` (NN) to `1.05` (maximally
separated). The far field is therefore at the order of magnitude of the product of two
single-point pivotal probabilities — the natural "two (quasi-)independent pivotal"
scale — while the NN field sits a factor 2.3–4.0 above it, with the factor *growing*
with L.

Why `unresolved` for the fusion reading:
* three sizes, and only five orbits (two at L=3) — no exponent may be written;
* the L=4 ordering is non-monotone in distance, showing the "profile" is still
  dominated by finite-size/torus effects;
* no arm classifier was run, so nothing here identifies a 6-arm or 8-arm object;
* `A_d` is order-1-bounded, so "compatibility" between the far-field level and any of
  {pivotal², 6-arm fusion, 8-arm fusion} magnitudes cannot be turned into an argument,
  in either direction. → `unresolved` for all three named scales.

### D3 — Is there a fact that directly opposes #768?

**Verdict: `unresolved` / nothing found.**

No configuration class in this atlas can be certified as "matching-odd and not removed
by signed balance": that classification needs the Phase-C arm witness, which was
deliberately not attempted (§T5). The nearest adjacent observation is that the **total
signed** pair interaction `Σ_{d≠0}J_d` is positive at all three sizes
(0.2596 / 0.1833 / 0.1366) and *decreasing*, while the absolute activity grows — so any
candidate "first non-common correction" that reads `M″` (or `Σ_d J_d`) directly as an
unsuppressed amplitude is, at these sizes, sitting inside a progressively stronger
cancellation. This is a **constraint**, not a counterexample, and I do not claim it
refutes or confirms any arm candidate in #768.

---

## §6 Files

Cloud (`/workspace/rev769/out/`, all downloaded to `rev769-out/`):

* `controls-20260914.json` — T1/T2/T3 exact-rational controls + identity checks
* `orbit-table-20260914.json` — T4 D4 orbit atlas (exact integer histograms + moments)
* `pivotal_L{3,4,5}_pathA.json`, `pivotal_L{3,4,5}_pathB_uf2.json`,
  `pivotal_L{3,4,5}_pathB_cov.json` — raw per-`d` integer histograms, three paths
* `witness_L{3,4,5}.json` — Phase C raw material
* `note-tables.md` — machine-rendered markdown tables used in this note

Scripts (`rev769-work/`, also on the cloud in `/workspace/rev769/scripts/`):

* `pivotal.c` (path A), `pivotal_bis.c` (path B), `witness.c`
* `analyze_pivotal.py` (exact-rational post-processing), `compare_paths.py`,
  `make_note_tables.py`

Reproduction (on the cloud VM; `PYTHONPATH` not needed for these):

```sh
gcc -O2 -march=native -pthread -o pivotal pivotal.c
gcc -O2 -march=native -pthread -o pivotal_bis pivotal_bis.c
gcc -O2 -march=native -o witness witness.c
for L in 3 4 5; do
  ./pivotal     $L --threads 16 --out out/pivotal_L${L}_pathA.json
  ./pivotal_bis $L --kernel uf2 --threads 16 --out out/pivotal_L${L}_pathB_uf2.json
  ./pivotal_bis $L --kernel cov --threads 16 --out out/pivotal_L${L}_pathB_cov.json
done
python3 compare_paths.py        # -> ALL_IDENTICAL
python3 analyze_pivotal.py      # -> controls-20260914.json, orbit-table-20260914.json
```

All computation was performed on `DevEnvC_551oUR`; nothing was computed on the Mac.

---

## §8 What I cannot claim

1. **No exponent is certified.** L=3,4,5 are three finite sizes. Every statement about
   the L-dependence of `cancel_ratio`, `A2_abs`, `A_d` is `compatible / incompatible /
   unresolved` only — never "measured exponent".
2. **The absolute influence is not a CFT operator.** `A_d = E|Δ₀Δ_dX|` and `A2_abs`
   are pure lattice diagnostics; #769 explicitly forbids calling them CFT operators,
   and I attach no field-theoretic identity to them.
3. **The identities of §T2 are consistency checks**, implied by Russo + translation
   invariance; they are not independent evidence about #768.
4. **`M″` is not "small" in absolute terms** — it grows with L; the small quantity is
   the *ratio* `|M″|/A2_abs`.
5. **The near/far structure is not an asymptotic two-region structure.** The L=5 torus
   diameter is `2√2`; the "far" region is 2–3 lattice spacings.
6. **Phase C is not a classifier.** No arm count, 6-arm, 8-arm or four-cluster label is
   claimed anywhere in this delivery.
7. **L=6+ was not computed** (not an acceptance condition). No MC was added; no STATUS,
   `#275` contract, or original data was modified.
8. Numbers quoted from #769/#775 are marked 引用; the T1 control values are **not**
   引用 — they were independently reproduced here.
