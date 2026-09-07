# Probe MID — Exact algebraic controls the new map stands on

**Issue #619 · C1–C10 · 2026-09-07 · exact, no Monte Carlo, no L=4 bond enumeration.**

This probe turns the four load-bearing facts of the new map into re-runnable
exact identities with JSON, tests, and failure counts. All programs ran;
every identity that should hold held; every "must not be ~0" quantity stayed
decisively away from 0.

---

## Standing result

| fact | claim | verdict |
|---|---|---|
| (1) `M_L(1/2) ≠ 0` | NN square-site is **not** self-matching | ✅ `-21/64` (L=3), `-13757/32768` (L=4) |
| (2) `Q(u)+Q(1-u)−1 ≠ 0` | matching involution is a two-model relation, not a single-model even/odd | ✅ max ≈ 0.13–0.18, not ~0 |
| (3) `p_L^H = Q(1/2) ≠ 1/2` | the live statement is `M(1/2)≠0` | ✅ `p_L^H` = 0.5865115… (L=3), 0.5906721… (L=4) |
| (4) `O = 1+(X−X²)/2`, `Cov(O,X)=P(r=2)=Var(X)/2` | #608 stable coupling = algebra + Harris | ✅ 0 config failures; all three equal at 0.287857056 |

Negative answers were the most useful outcome expected; none materialised — the
map's four supports all survived exact reproduction.

---

## C1 — L=3,4 site census reproduction

`exact_torus_enum.py` (imported from #606, not duplicated) and
`verify_independent.py` (two independent winding algorithms) both reproduced:

| L | configs | dual fail | (0,2)/(1,1)/(2,0) | `M(1/2)` | `p_L^H` |
|---|---:|---:|---|---|---|
| 3 | 512 | 0 | 259/162/91 | `−21/64` | 0.586511455113 |
| 4 | 65536 | 0 | 36559/19932/9045 | `−13757/32768` | 0.590672112331 |

Mismatch between the two winding algorithms = 0 at both sizes. Dest paths write
to `results/`, no `parents[2]` fix needed.

## C2 — `M_L` as an exact polynomial

`M_L(p)` in the power basis (exact rational coefficients), degree N=9 and N=16:

* L=3: `M(p) = −1 + 6p³ − 18p⁷ + 18p⁸ − 4p⁹`
* L=4: `M(p) = −1 + 8p⁴ + 32p⁶ − 64p⁷ + 172p⁸ − 704p⁹ + 1104p¹⁰ − 608p¹¹ − 56p¹² + 128p¹³ + 16p¹⁴ − 32p¹⁵ + 6p¹⁶`

Exact values:

| quantity | L=3 | L=4 |
|---|---|---|
| `M(0)` | `−1` | `−1` |
| `M(1)` | `+1` | `+1` |
| `M(1/2)` | `−21/64` | `−13757/32768` |
| `M'(1/2)` | `225/64` | `4209/1024` |

Uniqueness: Sturm sequence gives **exactly one root in (0,1)** at both sizes.

## C3 — Quantile table and the self-symmetry killer

`Q_L(u) = M⁻¹(2u−1)` (exact bisection on the rational polynomial, ≤1e-14):

| u | Q_L(u), L=3 | Q_L(u), L=4 |
|---|---|---|
| 0.1 | 0.324359162 | 0.378527592 |
| 0.2 | 0.413141491 | 0.451787961 |
| 0.3 | 0.479046015 | 0.504784420 |
| 0.4 | 0.535155140 | 0.549584290 |
| 0.5 | **0.586511455** | **0.590672112** |
| 0.6 | 0.636082791 | 0.630677097 |
| 0.7 | 0.686436757 | 0.671941457 |
| 0.8 | 0.741022052 | 0.717774932 |
| 0.9 | 0.807660296 | 0.776105597 |

* `Q(0.5) = p_L^H` exactly (same bisection). ✅
* `F_L(1/2) = [1+M(1/2)]/2`: L=3 `43/128`, L=4 `19011/65536`.
* **Self-symmetry kill**: `max_u |Q(u)+Q(1−u)−1|` = 0.1712 (L=3), 0.1803 (L=4) —
  decisively **not** consistent with 0 within 1e-12. Single-model self-symmetry
  is false; the matching involution is a two-model relation.
* `Q_L(u) − p_c` (p_c = 0.5927460 truncated, a cited convention, not a theorem):
  at u=0.5, L=3 `−0.0062345`, L=4 `−0.0020739`.

## C4 — Joint `(n_black, r_b, r_w)` and Alexander

* full joint histograms `(n_black, r_b)` and `(n_black, r_w)` for L=3,4 in
  `results/probe-exact-controls/joints.json`.
* `r_b + r_w = 2` on **every** config: `dual_fail = 0` (both sizes); no failure
  list to dump.
* Hand configs: the seven in-repo cases plus three NEW declared-before-running
  cases all satisfy `r_b + r_w = 2`. The checkerboard prediction `(r_b,r_w)=(1,1)`
  **lost to intuition**: the four corner sites form horizontal *and* vertical
  wrap cycles, giving `(r_b,r_w)=(2,0)` — a clean "the intuition loses" record.

## C5 — Two algorithms, every rank

Per-rank confusion tables (rows = algorithm A, cols = algorithm B) are diagonal:

* L=3: diag (259, 162, 91); L=4: diag (36559, 19932, 9045).
* `max |r_A − r_B| = 0`; `mismatch = 0` at both sizes.
* Runtime documented: L=4 enumeration 4.6 s, cross-check 3.3 s (well in budget).

## C6 — Square-bond L=3, wrap / `X` configuration-wise

Enumerated all `2^18 = 262144` square-bond configs (`O = wrap_either = int(r>0)`,
`X = r−1`):

* `O == 1 + (X − X²)/2` configuration-wise: **0 failures** (load-bearing identity
  holds exactly).
* Consequences at p=1/2 (uniform, exact counts):
  * `P(r=0,1,2) = 0.287857056, 0.424285889, 0.287857056` (note `P(r=0)=P(r=2)` —
    the duality `X ↦ −X` at self-dual p=1/2).
  * `E[O]=0.712142944`, `Var(O)=0.204995371`, `E[X]=0`, `Var(X)=0.575714111`.
  * **`Cov(O,X) = P(r=2) = Var(X)/2 = 0.287857056`** — all three equal. ✅
  * `Cov(O,X²) = −0.122133687`.

## C7 — Harris / nested-reveal, exact

Using #608's radius-ordered `spatial_filtration_levels` (4 levels at L=3) and the
exact conditional-mean `Γ_j = E[m_j m_jᵀ] − E[m_{j−1} m_{j−1}ᵀ]`:

* **Telescoping** `Σ_j Γ_j = Cov(Y)` to **0.000e+00** (exact). ✅
* `Γ_j` for `(O,X)` per level:

| j | revealed bonds | `Γ_j^{O,X}` |
|---|---:|---:|
| 0 | 4 | +0.039324621 |
| 1 | 12 | +0.103400728 |
| 2 | 14 | +0.040358663 |
| 3 | 18 | +0.104773045 |

* Every cross-increment `Γ_j^{O,X} ≥ 0` (min `+0.039`); Harris association holds
  at every scale — no negative increment to dump. `#608`'s single-signed ladder
  is **algebra + Harris, not a new mechanism**.

## C8 — Matching involution

Coupling used: site black = 4-connect (`G`), white = 8-connect (`Ĝ`, the matching
graph); `r_b + r_w = 2` configuration-wise (Alexander).

* `coeff_w[k] = −coeff_b[N−k]` exactly, for all `k` (L=3 and L=4). ✅
* `M_black(p) + M_white(1−p) = 0` as exact polynomials. ✅
* The three-line lemma `M_Ĝ(p)=−M_G(1−p)`, `F_Ĝ(p)=1−F_G(1−p)`,
  `Q_Ĝ(u)=1−Q_G(1−u)` therefore holds; and it **does not** imply
  `Q_G(u)+Q_G(1−u)=1` (C3 is the numerical support). Alexander and the
  involution were **not** conflated.

## C9 — Protocol card (≤40 lines, no scoring)

```text
declared chart     = span{1, Q_base, g_frozen}  unless HIGH probe overwrites
transport          = #612 identity before any residual
weightings         = spin0 AND equal, always together
forbidden          = second exponent; chart change after residual
Q_N(u) -> p_c      = location only (H1-H3); not a rate
M(1/2) != 0        = square-site is not self-matching   (C1/C2/C3)
Q(u)+Q(1-u) != 1   = single-model self-symmetry is false (C3)
#608 Cov(O,X)      = algebra + Harris, not a mechanism   (C6/C7)
p_L^H = Q(1/2) != 1/2   under F=(1+M)/2 and strict increase (C2/C3)
```

Not applied to N=725.

## C10 — Tests in the tree

`tests/test_probe_exact_controls.py` (7 tests, L=3 fast, L=4 slow-marked) —
all pass (`OK`, 1 slow skip). Covers C1 numbers, C2 fractions, C3 `Q(0.5)=p_L^H`
(1e-12) + self-symmetry kill (>1e-3), C5 mismatch=0, C6 failure count=0, C8
involution.

---

## Scope

No Monte Carlo; no L=4 bond (`2^32`) enumeration; no L=5; no arXiv; no N=725.
`B_even`'s triple joint `(O, X, B_even)` was left out of C6 because `B_even`
requires the `geometric_dual_mask` path (`square_bond_duality_exact`), which is
out of the site enumerator's scope; the `(O,X)` skeleton above is the typed
control #608 needs, and `B_even` can be appended when that path is imported.
Nothing here edits `docs/STATUS.md`; no ticket is closed.
