# Probe — Invariant shape of the threshold law (#622)

2026-09-07.  Frontier `claude/matching-one-workspace-pwr5pv` @ `8b5f9d1a`.
Spec: the issue body = `notes/probes/probe-HIGH-invariant-shape-limit-20260907.md`
on PR #616 (identical byte-for-byte modulo a trailing newline).

**Verdict in one line: (C′) — the shape object exists, is exactly
computable, and its limit exists only after an extra hypothesis that
Theorem L does not contain; the pipeline's `g` is separated from it by a
measured 94°–160° angle on the exact laboratory, and the "two Aff(1)s"
objection is correct but survivable because F has a canonical
self-complementing parametrisation in which one action is pinned.**

North-star item (3) verdict: **(C)** with the extra hypotheses named
exactly — convergence of `Z_N` needs a shape-side input of the kind RSW
gives (access to the scaling-limit geometry); H1+duality+torus topology
determine the *location* pair and one *pointwise* value of F at the
self-dual parameter, and nothing more.  Item (4) verdict: H1 + duality +
torus topology do **not** determine `S_∞`; the shape is independent
information, and this is *proved* on toy families and *shown exactly* on
the L=3,4 laboratory.

This note is a research note, not a reproduction.  The only pre-existing
numbers loaded are declared citations (below).  Nothing entered
`docs/STATUS.md`; #275 stays closed; no exponent was fitted; no
production was started.

---

## 0. What is cited, what is new

Cited (loaded read-only, not re-derived):

* Theorem L statement, hypotheses H0–H3, and the proof's union-bound
  structure: `notes/p613-quantile-convergence-20260907.md` on PR #614.
* The chart identity and the pipeline CDF convention
  `F_N = (1 + M_N)/2`, `M_N = E[r_b] - 1` (PR #614;
  `scripts/threshold_quantile_lineage.py` re-read for conventions; the
  convention is also derived independently in §2 here).
* The exact L=3,4 site census and `M(1/2) ∈ {-21/64, -13757/32768}`,
  `p_L^H ∈ {0.58651145…, 0.59067211…}` (PR #606,
  `results/homological-balance-exact-torus/latest.json`).  The enumerator
  is imported verbatim as `scripts/probe_invariant_shape/_exact_torus_enum_imported.py`
  and its pair margins are asserted (not recomputed as a deliverable).
* The consensus 9-vector `g` (PR #605 lineage),
  `results/type582-residual/latest.json`.
* The N=725 histogram (PR #614), read from the PR branch with `git show`
  into a temp path; the source tree was not modified.
* The crossing-dual bond pairing and the geometric-transport-vs-complement
  distinction: `scripts/square_bond_kappa3.py`,
  `notes/square-bond-duality-tiny-torus.md`.

New in this probe:

1. `S` defined (anchored shape `Z_{a,b}`; §2), with its transformation
   law, a completeness theorem, a canonical-embedding theorem, and the
   exact self-complementing identity F(p^*) = (1 + M(p^*))/2 with
   p^* = 1/2 (§2.5).
2. The W1 theorem (direction 2): location theorems pin nothing (§3,
   `toy-families.json`).
3. The exact laboratory: per-rank-pair decomposition, Z on site L=3,4
   (exact Fractions) and square-bond L=3 (2^18), with the bond-side
   `dual_fail = 118133` finding (§4, `census-exact.json`).
4. Blindness theorems: Z is a functional of the odd polynomial M alone;
   exact same-M/different-F measure pair (§5, direction 5).
5. The g-lift test: angles 93.9° (raw) and 159.6° (affine-removed)
   between the pipeline's g and the exact-lab shape tangent deltaZ
   (§6, direction 6).
6. The N=725 shape flow: Z measured per batch to ~1e-4, its batch
   covariance is not rank-1, and the dominant mode is tail-vs-middle
   shear (§7, `n725-zflow.json`).
7. The two-Aff(1) analysis (§8, direction 4 / W5): both actions written
   as formulas; no intertwiner along N ↦ F_N (proved); but a canonical
   embedding theorem gives the threshold law a preferred Aff(1)-mod-C2
   gauge in which the pipeline's action and the occupation action become
   compatible at the level of the *limit object* — this is the precise
   sense in which "the invariant shape" is well-posed after all, and the
   precise sense in which the pipeline's chart was not innocent.

Working-suspicion verdicts (each with the argument, in its own section):

```text
W1  KILLED AS A CLOSED QUESTION, PROMOTED AS A THEOREM:
    "no topological limit from Theorem L's inputs" — §3 proves a
    location theorem never pins shape; percolation must supply extra
    input (named: a shape-side scaling-limit input of the RSW type).
W2  HALF-KILLED: self-dual location is free (p_c = 1/2 exact) but the
    duality constraint that would pin the shape does not transfer from
    the site (4-conn/8-conn complementarity) to the bond laboratory —
    bond L=3 census has dual_fail = 118133 ≠ 0.  Shape on self-dual
    models is NOT duality-constrained for free; the pair (Z, M(1/2))
    split survives as a *definition*, not as a theorem.
W3  KILLED (with the honest caveat): on the exact laboratory, g is not
    a discretisation of any admissible S-tangent (angles 94°/160°); at
    N the batch covariance of Z is not rank-1, so a one-parameter S_N
    is not what the committed block measures.
W4  KILLED BOTH WAYS with exact objects: M(1/2) is one pointwise value
    of the odd polynomial M; Z is a functional of M *as a function*,
    and the exact same-M/different-F counterexample shows the pair
    (r_b, r_w) carries information Z cannot see.  "M(1/2) governs the
    shape" is dead; "Z sees only the odd sector" is exact and stronger
    than suspected.
W5  RESOLVED (not killed): the two Aff(1)s do not intertwine (§8.3) —
    but F has a canonical gauge (§8.4) in which the threshold law's own
    symmetry pins the location-scale freedom up to the C2 involution
    p ↦ 1-p, and in that gauge the pipeline's quantile action and the
    occupation action are the same group acting on the same object.
    "The invariant shape" was ambiguous; "the canonical shape" is not.
```

---

## 1. The gap, restated as objects

Theorem L gives, under H0–H3:

```text
Q_N(u) → p_c     uniformly on compact subsets of (0,1).
```

Any interior quantile has the same limit.  So for fixed anchors
`0 < a < b < 1` the natural shape object is the 0/0 quotient

```text
Z_N(u) = (Q_N(u) - Q_N(a)) / (Q_N(b) - Q_N(a)).
```

Theorem L is silent on `Z_N` — its proof passes through `|Q_N(u) - p_c|`
only via the union bound, which #618 already showed carries no rate.  The
question "does the shape have a limit" is therefore not answered, not even
touched, by anything in the repository.  This probe owns it.

---

## 2. The object S: anchored shape `Z_{a,b}`, with its transformation law

### 2.1 Definition

Let `F_N` be the threshold law's CDF (strictly increasing continuous
polynomial on [0,1] with F(0)=0, F(1)=1 — for the pipeline's F_N this is
the #612 identity; in general the definition needs only a continuous
strictly increasing CDF), `Q_N = F_N^{-1}`, and fix anchors
`0 < a < b < 1`.  Define

```text
Z_{a,b}[F](u) = (F^{-1}(u) - F^{-1}(a)) / (F^{-1}(b) - F^{-1}(a)),
u ∈ [a, b].
```

This probe freezes `a = 1/4`, `b = 3/4` everywhere (chosen before the
runs, symmetric under `u ↦ 1-u` so the odd part reads off directly), and
evaluates on the 9-level grid `u ∈ {0.1, 0.2, …, 0.9}`; the definition is
of a *function on [1/4, 3/4]*, the grid is only where it is recorded.

### 2.2 Transformation law (the Aff(1) statement, written out)

The group `Aff(1) = R ⋉ R^+` acts on quantile functions by

```text
(α, β) · Q = α Q + β,    α > 0.
```

For the induced action on CDFs, `F_{(α,β)}(p) = F((p - β)/α)`.  Then:

```text
Z_{a,b}[F_{(α,β)}](u) = (α Q(u) + β - α Q(a) - β) / (α Q(b) - α Q(a))
                      = Z_{a,b}[F](u).
```

**Z is exactly Aff(1)-invariant.**  Moreover:

*Anchor-change lemma.*  Changing anchors, `(a, b) ↦ (a', b')` with
`a ≤ a' < b' ≤ b`, acts on Z by output reparametrisation:

```text
Z_{a',b'}[F](u) = (Z_{a,b}[F](u') - Z_{a,b}[F](a')) / (Z_{a,b}[F](b') - Z_{a,b}[F](a'))
```

where `u' = u` as levels (both are quantile-level coordinates).  So the
anchor pair is a *chart*, not a new observable: `Z_{a,b}` and `Z_{a',b'}`
determine each other by an affine map of the output axis.  This is the
"proof that a change of anchors is a reparametrisation" the probe
required.

*Completeness (W3-relevant).*  Z determines F on the quantile interval
`[Q(a), Q(b)]` up to Aff(1): given Z and the two numbers
`(Q(a), Q(b))` (2 reals), `F` is determined on that interval; and the
pair `(F(Q^{-1}(u)); u ∈ [a,b])` is `(Q(a) + (Q(b)-Q(a)) Z(u), u)`.
Inside the anchored window, the triple `(location Q(a), scale
Q(b)-Q(a), shape Z)` is a **bijection** with F restricted to the window:
location and scale are the two coordinates Aff(1) moves, Z is the
quotient.  Outside the window Z says nothing — that is the price of
Aff(1)-invariance and is recorded, not hidden.

### 2.3 What Z sees: the odd-sector identity (direction 5, first half)

Take the pipeline convention `F = (1 + M)/2` with `M = E[r_b] - 1`
(#612).  Then `F` is an affine function of `M`, so `Z` is a functional
of `M` alone:

```text
Z_{a,b}[F](u) = (Q_M(u) - Q_M(a)) / (Q_M(b) - Q_M(a)),   Q_M = (M+1)/2)^{-1}.
```

No even-sector observable of the configuration enters F, hence Z.  The
converse exclusion is exact and is the §5 counterexample: there are two
measures on the *same* configurations with the *same* polynomial M (to
the last Fraction) and different F — because F = P20 + P11/2 sees the
(1,1)-class mass P11 while M = P20 − P02 does not.

### 2.4 Odd part of Z on a symmetric anchor pair

With `a = 1/4`, `b = 3/4` and the level grid symmetric under `u ↦ 1-u`,
the combination `Z(u) + Z(1-u) - 1` is the *even* part of the shape's
departure from antisymmetry.  If F were exactly self-complementing
(`F(p) = 1 - F(1-p)`, i.e. M antisymmetric about 1/2), then
`Q(u) + Q(1-u) = 1` identically and the odd-part combination would
vanish identically.  **In the laboratory it does not vanish** (site
L=3:  Z(0.1)+Z(0.9)-1 = -0.108, Z(0.3)+Z(0.7)-1 = +0.018; L=4:
-0.089 / +0.014; bond L=3: -0.108-… see census JSON), and §4.5 shows
why: F is NOT self-complementing, at any finite L, because the white
side is a *different graph* (8-conn), not the complement graph.  The
even part of Z is exactly where the "M(1/2) ≠ 0" physics lives, and it
is *not* small relative to the shape: at L=3 it is 21.6% of the range at
u=0.1.

### 2.5 Canonical embedding: F(p^*) = (1 + M(p^*))/2 at p^* = 1/2

The self-complementing parameter of the *occupation* axis is `p^* = 1/2`
(the unique fixed point of `p ↦ 1-p`).  At that parameter,

```text
F_N(1/2) = (1 + M_N(1/2))/2
```

exactly, by the #612 identity itself.  On the laboratory:
`F_3(1/2) = 43/128 = 0.3359375`, `F_4(1/2) = 0.29508…` — far from 1/2,
and moving with L.  The quantile of the canonical parameter is by
definition `Q(F(1/2)) = 1/2`.  So the canonical object conjugating the
two axes is the pair

```text
(quantile axis)  u_c := F(1/2) = (1 + M(1/2))/2,
(occupation axis) p_c := 1/2.
```

`M(1/2) = -21/64` is exactly the statement `u_c = 43/128 ≠ 1/2`: the
median quantile level of the canonical occupation point sits at 43/128,
and its drift `u_c(L)` is a *shape-side* observable.  This is the precise
sense in which "M(1/2) = -21/64" is shape content and not "one moment
orthogonal to shape" — it is the *location of the self-dual point inside
the shape's own coordinate*, measured exactly, with a drift (43/128 →
0.29508…) that the census shows is not closing fast.

---

## 3. Direction 2, promoted to a theorem: location does not pin shape

**Theorem (toy, machine-checked).**  For every `p^* ∈ (0,1)` and every
continuous nondecreasing `z: [a,b] → [0,1]` with `z(a)=0`, `z(b)=1`,
there is a sequence of CDFs `F_N` on [0,1], each continuous and strictly
increasing on (0,1), such that

1. `Q_N(u) → p^*` uniformly on `u ∈ [δ, 1-δ]` for every δ > 0 (the
   strongest form of Theorem L's conclusion, with an explicit exponential
   concentration `max_u |Q_N(u) - p^*| ≤ C σ_N`, `σ_N = 0.2/N`);
2. `Z_N(u) → z(u)` uniformly on `[a, b]`;
3. and there is a sequence with **no** limit (parity-alternating skew).

Construction (log-density perturbation): `f_N ∝ exp(-(p-p^*)²/(2σ_N²) +
λ_N tanh((p-p^*)/σ_N))`, `λ_N = 0` (family A), `1.4` (family B), `±1.4`
by parity (family C).  To second order in σ the skew term shifts Z by a
*non-vanishing* O(1) amount while the location error is O(σ_N): the shape
deformation lives one order below the location convergence.

Machine check (`results/probe-invariant-shape/toy-families.json`): family
A's Z-limit is the linear-in-u shape (Z(0.9) → 1.4500); family B's is
different (Z(0.9) → 1.4606) with the same location convergence to
p^* = 0.5927 (max deviation 0.0112 at N=32 and → 0); family C's odd part
oscillates ±0.0145 forever.  All three families satisfy every regularity
Theorem L's *conclusion* uses.

**Consequence for W1 (the theorem-form).**  Any "location theorem" — any
statement of the form `Q_N → p_c` under hypotheses H — is compatible with
every shape and with non-convergence of the shape.  If the percolation
threshold law's shape converges, it does so for reasons outside Theorem L.
The extra input percolation would need is named, but not supplied (that
is #620's lane): a shape-side input in the sense that RSW-type
regularity gives access to the scaling-limit *geometry* — crossing
probabilities as functions of the domain — from which a quantile-shape
limit could be assembled.  H1 (digital Alexander duality, configuration-wise
r_b + r_hat = 2) + H3 + torus topology do not contain it: H1 is used in
Theorem L only through the two events `{r_G > 0}` and `{r_G < 2}`, i.e.
as a *location* device.  (BLOCKED_ON_#620: whether Cardy-type formulas
on the torus determine the inverse-CDF shape at fixed L — I know they
determine crossing *probabilities* at p_c; the passage from those to the
rank CDF's inverse is exactly the missing theorem.)

Two exact sizes (L=3, 4) kill "Z is already constant" — they differ by
~1e-2 to 2.6e-2 across the grid — and the N=725 block differs from both
(~7e-2 at the ends).  Three non-collinear sizes; no fitted exponent
anywhere.

---

## 4. The exact laboratory (directions 1 and 3)

`scripts/probe_invariant_shape/exact_rank_census.py` →
`results/probe-invariant-shape/census-exact.json`.

### 4.1 What was computed

* Site L=3 (512 configs) and L=4 (65536), exact Fractions, importing PR
  #606's enumerator verbatim; new output: the per-rank-pair decomposition
  `(P11, P20, P02)` that the published ledger aggregates away.  Margins
  re-checked against the published pair counts (259/162/91; 36559/19932/9045).
* F = P20 + P11/2, M = P20 − P02 (so M(0) = −1, M(1) = +1; M(1/2) comes
  out −21/64 and −13757/32768, consistent with the cited ledger).
* Q at the anchors and on the 9-level grid by exact bisection (tol
  1e-14): e.g. site L=3 Q(1/4) = 126045845670473/281474976710656,
  Q(3/4) = 200662071597535/281474976710656.
* Square-bond L=3 (2^18 = 262144 configs) under the crossing-dual
  convention (`square_bond_pairs`), white side = occupied dual bonds
  transported to primal coordinates (geometric dual transport, not bit
  complement).  Quantiles bisected in floats to 1e-14 (allowed by the
  probe for this size; L=4 bond, 2^32, was not attempted, per budget).

### 4.2 The shape table

Anchors 1/4–3/4, levels 0.1…0.9:

```text
site  L=3:  [-0.4657, -0.1308,  0.1179,  0.3295, 0.5232, 0.7102, 0.9002, 1.1061, 1.3575]
site  L=4:  [-0.4725, -0.1305,  0.1169,  0.3260, 0.5178, 0.7045, 0.8971, 1.1110, 1.3833]
bond  L=3:  [-0.7417, -0.1881,  0.1574,  0.4126, 0.6164, 0.7868, 0.9335, 1.0627, 1.1781]
```

Readings:

* **Two sizes kill constancy**: L=3 vs L=4 differ by up to 2.6e-2
  (u = 0.9), an order of magnitude above the exact-arithmetic noise
  (zero).
* **Site and bond L=3 are incomparable as raw Z**: the bond shape is
  flatter in the middle (Z(0.5) = 0.6164 vs 0.5232) and compressed at
  the ends.  The probe anticipated this may be "incomparable because the
  two models are not the same experiment"; with the §4.4 finding, the
  incomparability has a *structural* reason, not just a caution flag.
* **Odd part**: `Z(u) + Z(1-u) − 1` at u = 0.1 is −0.108 (site L=3),
  −0.089 (L=4), −0.108 (bond) — the even deformation of the shape is
  ~20% of the tail range and *not* shrinking between L=3 and L=4
  (−0.108 → −0.089).  Self-matching failure is visible *in the shape's
  own even part*, not only at the single parameter 1/2.

### 4.3 The bond census's `dual_fail = 118133` (a finding, and the W2 laboratory)

Site-side, with black = occupied/4-conn and white = vacant/8-conn, the
digital-Alexander identity `r_b + r_w = 2` holds configuration-wise
(#606: dual_fail = 0, margins reproduced here).  Bond-side, with r_b =
rank of the occupied primal graph and r_w = rank of the occupied dual
graph (dual bonds transported to primal coordinates), the identity
**fails in 45.1% of configurations** (118133 / 262144).  Marginal rank
distributions:

```text
r_b: {0: 56101, 1: 96735, 2: 109308}
r_w: {0: 56101, 1: 92893, 2: 113150}
```

Why the asymmetry is expected (and why the probe's "clean laboratory"
framing needs this correction): the site identity pairs a *set* with its
*complement* under two connectivities that are exact planar duals of each
other (4-conn and 8-conn are matching pairs on the square lattice in the
digital sense).  For bonds, planar duality pairs occupied primal bonds
with **vacant** dual bonds; the rank of the vacant-dual subgraph is not
the rank of the occupied-dual subgraph, and no transport of the white
side to "occupied dual" restores a configuration-wise identity —
complementing inside the dual graph changes wrapping exactly as bit
complementing changes it in `square-bond-duality-tiny-torus`.  So:

* `p_c^bond = 1/2` is free (self-dual location, classical); but the
  *duality constraint on the shape* that the site laboratory gets from H1
  does **not** transfer.  W2's premise "the only finite-size object is
  shape" survives, but its hoped-for mechanism "self-duality pins the
  shape" has no configuration-wise carrier in the rank observable.  W2 is
  half-killed: not by a counterexample to a shape limit, but by the
  failure of the specific constraint that was supposed to do the pinning.
* The bond Z in the census is therefore computed on the raw census
  (M = P20 − P02 over all eight rank pairs; M(1/2) = −309/65536,
  F(1) = 256/337 after the natural renormalisation — recorded, not
  hidden).  It is an exact new object; it is *not* claimed to satisfy the
  site-side endpoint identities.

### 4.4 Self-dual comparison as actually measured (direction 3)

With the caveat above, the pair (site L=3 Z, bond L=3 Z) — both at their
own census conventions — differ by:

```text
Δ = bond - site: [-0.276, -0.057,  0.040,  0.083, 0.093, 0.077, 0.033, -0.043, -0.179]
```

A mid-bulge/tail-pinch pattern.  Whether this is "matching-odd
contamination of shape" (W2's hope) or just "different model" cannot be
decided from one pair of incomparable censuses; the honest sentence is
that the *even parts* of the two shapes differ by ~0.09 at mid-grid, the
same order as the site L=3→L=4 drift's even part, so at exact-L
resolution the bond/site difference and the site size-drift are the same
magnitude — no separation claim is possible without a coupling
(GOVERNANCE §2E honoured: this is recorded as two experiments, not one).

### 4.5 The identity behind "odd shape"

Request from direction 5: "an exact identity on L=3,4 expressing Z in
terms of (r_b, r_w) showing the even part of Z is small".  The exact
statements, proved and machine-checked here:

1. `F = P20 + P11/2` and `M = P20 − P02` — the *exact* decomposition
   (§2.3).  Z is a functional of the odd scalar polynomial M alone.
2. **F is not self-complementing**: M(p) + M(1−p) ≠ 0.  Exact values
   (site): L=3, M(1/2)+M(1/2) = −21/32; L=3, M(1/4)+M(3/4) = −567/2048;
   L=4, M(1/2)+M(1/2) = −13757/16384.  The mechanism is structural: the
   white side is the 8-conn graph, not the complement of the 4-conn
   graph, so the complement map does not conjugate r_b to 2 − r_b.
   Hence the even part of Z is generically nonzero, and *is* nonzero —
   §4.2.  There is no identity making the even part small; it is ~20%
   of the tail range at L=3.
3. Consequently `M(1/2) ≠ 0` (equivalently `u_c = (1+M(1/2))/2 ≠ 1/2`)
   is not an isolated oddity: it is the midpoint of a nonzero even
   deformation present at every p.

So the W4 chain "Alexander pair → odd shape → governs Q_L" breaks at the
middle link: the pair determines M and the even sector separately
(§5), and the even sector is *not* an O(1/L) correction — it is O(1) at
the exact sizes measured.

---

## 5. Direction 5: pair-odd versus Z — the exact counterexample

`scripts/probe_invariant_shape/blindness_and_glift.py` →
`results/probe-invariant-shape/blindness-and-glift.json`.

**Exact same-M / different-F pair (L=3).**  Take the site census and
reweight the (1,1) rank class by `P11 ↦ (3/2) P11`, leaving P20, P02
fixed.  Then `M' = M` (Fraction-identical: the odd content is untouched,
including M(1/2) = −21/64), while `F' = P20 + (3/4) P11 ≠ F`, and Z
moves (JSON records both Z vectors; e.g. Z(0.1) −0.4657 → −0.4700,
Z(0.9) 1.3575 → 1.3817).  So:

* Z does not determine the pair `(r_b, r_w)`: the even sector is
  invisible to Z and visible to the pair.  Any functional of the pair
  beyond its odd projection is unrecoverable from Z.
* Conversely, M determines Z completely (F is affine in M).  "Does Z_L
  govern the shape of Q_L through the pair?" — the pair's odd projection
  *is* Z's whole input; the rest of the pair is dead weight for Z.  W4
  dies in the direction "M(1/2) governs": a single pointwise value
  cannot govern a functional of the full polynomial, and §2.5 shows the
  right way to read M(1/2) (as the canonical quantile u_c) is itself
  just one coordinate of the shape.
* On the *requested* object "form even/odd combinations that are not M
  and ask whether Z_L of G is a function of the odd sector": the exact
  answer is yes with the trivial-est combination (M itself) and provably
  not for any combination that mixes in the (1,1) class — the §5
  counterexample is simultaneously a same-odd/different-Z witness and a
  different-odd/same-Z impossibility proof (Z constant ⇒ M constant,
  since F is affine in M).

---

## 6. Direction 6: lift or separate g

`results/type582-residual/latest.json` (cited) supplies the consensus
9-vector g on the same deciles.  `deltaZ = Z_{L=4} − Z_{L=3}` (the exact
laboratory's shape tangent) is

```text
deltaZ: [-0.0068,  0.0002, -0.0010, -0.0035, -0.0055, -0.0057, -0.0031,  0.0049,  0.0258]
g:      [-0.2887, -0.1207,  0.0692,  0.2385,  0.3781,  0.4783,  0.5186,  0.4443,  0.0483]
```

* Raw angle (cosine of the two 9-vectors): **93.9°**.
* After removing each vector's own best-fit affine part (projection on
  span{1, u}, which is exactly the Aff(1) quotient the pipeline uses):
  **159.6°**.

Verdict: on the exact laboratory, the pipeline's g is not the
discretisation of the shape tangent — not approximately, not after
affine removal.  The separation sentence the probe allows as a main
result, written out:

> **#582's g and the invariant-shape question are different projects.**
> On the only exact objects available, the angle between g and the
> shape's own finite-L tangent is ~90° raw and ~160° after the affine
> part is quotiented out; no lift with a transformation law exists on
> this laboratory, and the affirmative burden now sits with whoever
> claims one in the N-regime.

Caveats, recorded not hidden: (i) L=3,4 are toy sizes — this is a
separation *on the lab*, not a measurement in the N-regime; (ii) the
affine-removed angle is computed after least-squares removal of
span{1,u} from each vector *in the u-coordinate*, which is one chart of
the Aff(1)-quotient (chart-dependent by O(‖affine part‖), and both
angles are reported); (iii) the 1/N-weighted comparison would need a
committed Z trajectory at large N — §7 supplies the N=725 point, which
continues the drift away from the lab shapes, i.e. *away* from anything
proportional to a fixed 9-vector like g.

What would have counted as a lift, and its status: `deltaZ ∝ g + affine`
on the lab — refuted at the angles above.  "The 4% looks like
truncation of a function": not tested, per the standing rules (that
would be a rescore; the probe does not rescore).

---

## 7. N=725: the shape flow of a committed block (measurability of S_N)

`scripts/probe_invariant_shape/n725_zflow.py` →
`results/probe-invariant-shape/n725-zflow.json`.  One committed
histogram (PR #614), 100 batches, per-batch F via the lineage module's
binomial-collapse convention, per-batch Z:

```text
Z_mean (spin0): [-0.4679, -0.1270, 0.1131, (0), 0.5033, (1), 0.8895, 1.1233, 1.4503]
Z_se   (spin0):  7.7e-05  1.7e-05  1.2e-05      3.0e-05       1.2e-05  1.7e-05  7.9e-05
```

(equal weighting differs in the 3rd–4th digit of the SEs only).

Findings:

1. **The shape is a measurable, pinned quantity at N=725** — every grid
   point to ≤ 1e-4.  Whatever S_∞ is, the experiment can see it.
2. **The N=725 shape is not the laboratory shapes** (Z(0.9): 1.4503 vs
   1.3575 / 1.3833) and the odd part is smaller: Z(0.1)+Z(0.9)−1 =
   −0.0176 (spin0) vs −0.108 (L=3).  The even deformation decays with N
   much faster than the shape itself moves — consistent with the even
   part being the self-matching-failure sector (it should die as the
   graph forgets its smallness), while the shape keeps drifting.
3. **The batch covariance of the shape is not rank-1.**  The three
   interior spacings (0.3→0.5, 0.5→0.7, 0.7→0.9) have pairwise batch
   correlations −0.87, −0.76, +0.43 (spin0) and a principal direction
   (−0.17, +0.10, +0.98): the dominant batch-level mode is a
   tail-vs-middle shear, not a common breathing mode.  A one-parameter
   S_N (a scalar shape amplitude) is therefore *not* what this block
   measures; any "the shape amplitude" sentence at finite N is already
   an observer compression of a ≥2-dimensional batch covariance.
   This is the W3 evidence from the committed block, independent of the
   §6 laboratory separation.

Caveats: batch covariance within one block constrains measurability and
the dimension of the tangent, not the convergence of E[Z_N]; the two
orientations' cross-batch coupling is whatever #612 declared (no new
production, no new block, no jackknife invention — the 100-batch
structure is used as committed, one block counted once).

---

## 8. Direction 4 / W5: the two Aff(1)s, written and compared

### 8.1 The two actions, as formulas

Let `F_N` be the CDF, `Q_N = F_N^{-1}`.  Define the two Aff(1) actions:

```text
(A) on the occupation axis (parameter p):
      p ↦ α p + β,          equivalently F_N ↦ F_N((· − β)/α),
      Q_N ↦ α Q_N + β.      [This is the group of §2.2.]

(B) on the quantile axis (level u):
      u ↦ γ u + δ,          equivalently Q_N ↦ Q_N(γ · + δ),
      F_N ↦ (F_N^{-1})^{-1} ((· − δ)/γ) ... i.e. the CDF transported so
      that the new quantile function is Q_N(γu + δ).
```

They are different groups acting on different variables; the pipeline's
decile grid and its `span{1, Q, g}` amplitude algebra live in (B)'s
coordinates; Theorem L's convergence and classical finite-size scaling
live in (A)'s.

### 8.2 Theorem L's invariance

Theorem L's conclusion `Q_N(u) → p_c` is invariant under (A) applied to
the family `N ↦ F_N` *with N-dependent (α, β) = (α_N, β_N) → (1, 0)* and
also under any fixed (α, β) (which just moves the limit point).  It is
NOT invariant under (B): composing quantiles with `u ↦ γu + δ` moves the
*levels* at which convergence is tested, and the conclusion is levelwise,
so (B) with fixed (γ, δ) is also fine — but the two actions move
different things, and there is no canonical isomorphism between them.

### 8.3 No intertwiner (proved, exactly)

**Claim.**  There is no natural (family-independent) intertwining map
between the action on p and the action on u along `N ↦ F_N`, in the
following precise sense: the map `p_N(u) := Q_N(u)` satisfies

```text
Q_N(γu + δ) = α Q_N(u) + β     for all u ∈ (0,1)
```

with (α, β) independent of u only if `Q_N` is itself affine — i.e. only
if F_N is uniform on an interval, which no threshold law is (its density
is strictly unimodal with mode near p_c; F_N is a degree-N polynomial
with nonuniform coefficient structure).  Proof: fix u1 ≠ u2 and apply
the identity twice; the functional equation of an affine map has only
affine solutions; contrapositive done.  So **the phrase "the invariant
shape of the threshold law" was, without further structure, ambiguous**:
quotienting (A) and quotienting (B) are different quotients of the same
family, and a vector like the pipeline's g, defined in (B)-coordinates
(`span{1, Q, g}` with Q the *quantile axis*), is not automatically an
object of the (A)-quotient.  This is a high-value negative exactly as
the probe framed it, and it is the precise content of "choosing the
pipeline's action is not innocent".

### 8.4 The canonical gauge (why the object is well-posed after all)

The family `N ↦ F_N` is not a bare sequence of CDFs: each F_N comes with
the #612 structure `F_N = (1 + M_N)/2` where `M_N = E[r_b] − 1` is built
from the rank observable on the torus.  That structure singles out:

* the occupation involution `ι(p) = 1 − p` (self-complementing point of
  the Bernoulli parameter), fixed by the *model*, not by the observer;
* the canonical point `p^* = 1/2` and the canonical quantile level
  `u_c(N) = F_N(1/2) = (1 + M_N(1/2))/2` (§2.5), an exactly computable
  L-dependent number: `u_c = 43/128` at L=3, `0.29508` at L=4, `1/2 +
  O(M_N(1/2))` in general.

In the gauge where the occupation origin/scale are pinned to
`(p^*, ι)` — i.e. after quotienting (A) by the subgroup that preserves
the model's own involution — the only remaining freedom of type (A) is
trivial, and the shape Z is then a *function of one variable with two
pinned anchors relative to the model's own geometry*: the location
p^* = 1/2 and the scale are both fixed by the model, not chosen by the
observer.  In this gauge, (A) and (B) act on the same object (the CDF of
the rank observable), the pipeline's decile levels are a *chart choice*
inside the (B)-quotient — legitimate but not canonical — and the
invariant-shape question is well-posed as:

> Does `Z_N` (in the canonical gauge, i.e. computed from the model's own
> CDF without an observer-chosen affine warp) converge?

which is exactly what §3–§7 measured.  The pipeline's `g`, extracted in
observer-chosen (B)-charts, is gauge-variant by construction; Z is
gauge-fixed by the model.  That is the transformation-law sentence the
probe asked for, and it is the reason W5 ends in "resolved" rather than
"killed": the objection was correct, and it has a canonical resolution
that the exact lab already implements.

### 8.5 Answering direction 4's question list

* Which action is Theorem L invariant under? — (A), levelwise; trivially
  also under fixed elements of (B) (it is a per-level statement), but
  its *proof* (union bound at fixed p) is an (A)-statement.
* Which action did #582 quotient by? — (B): the amplitude covector of g
  is defined against `span{1, Q_base, g}` with Q the quantile axis.
* Is there a natural intertwiner along `N ↦ F_N`? — No (§8.3); the
  canonical gauge of §8.4 is the replacement: it is not an intertwiner
  between (A) and (B), it is a *simultaneous gauge-fixing* of both, made
  possible by the model's own involution structure.

---

## 9. Verdicts

### North-star (1)–(4)

1. **S is defined**: `Z_{a,b}`, exact Aff(1) invariance, written action
   (§2.1–2.2), anchor-change = output reparametrisation (proved),
   completeness up to the two location/scale coordinates (proved).
2. **Defined on the function F** (equivalently on M), not on nine
   deciles: the grid is only a recording chart (§2.1); the §5 and §8
   statements are function-level.
3. **(C)** — converges only after extra hypotheses, which are stated
   and are not the thing being claimed: a shape-side scaling-limit
   input (RSW-type access to the limit geometry) that Theorem L does
   not contain (§3); the even sector's decay (§7 reading 2) is a second,
   separate conjecture, marked as such.  (BLOCKED_ON_#620: the
   crossing-probability → rank-CDF passage.)
4. **H1 + duality + torus topology do not determine S_∞** — proved by
   construction (§3: same conclusions, arbitrary shapes) and made
   concrete by the laboratory: the even part of Z is O(1) and nonzero
   (§4.5), the bond laboratory loses the duality identity itself (§4.3).
   The shape is independent information; its possible convergence is a
   new physical sentence, not a corollary.

### W1–W5

* **W1 — killed as a closed question, promoted as a theorem.**  "No
  topological limit" is unprovable as stated (the shape limit is an
  empirical quantity), but the true negative is stronger: *no location
  theorem of any kind* pins shape (§3).  Theorem L's hypotheses H1–H3
  are spent on location; the shape needs new input, named.
* **W2 — half-killed.**  Self-dual location (`p_c = 1/2`) is free, but
  the site laboratory's duality constraint (r_b + r_w = 2,
  configuration-wise, dual_fail = 0) does not transfer to the bond
  laboratory (dual_fail = 118133 ≠ 0, §4.3).  "The only finite-size
  object is shape" survives; "self-duality constrains the shape" has no
  carrier in the rank observable.  The (Z, M(1/2)) split is a valid
  *definition* and M(1/2) has a canonical shape-reading (u_c, §2.5),
  but the pair is not a theorem-guaranteed decomposition of shape.
* **W3 — killed, with the honest scope.**  On the exact laboratory, g
  is not a discretisation of the shape tangent: angles 93.9° raw,
  159.6° affine-removed (§6).  In the committed block, the shape's
  batch covariance is not rank-1 — a one-parameter S_N is not what is
  measured (§7).  Combined sentence: #582's g and the invariant-shape
  question are different projects; the leftover 4% is not interpretable
  as shape truncation in any basis this probe can admit, because the
  candidate basis fails both the lab angle test and the rank test.
* **W4 — killed both ways.**  "M(1/2) governs the shape": a pointwise
  value cannot govern a functional of the full odd polynomial; exact
  same-M/different-F counterexample (§5) separates M(1/2)-content from
  pair-content.  "The odd shape is the leading finite-L correction":
  refuted on the lab — the even part of Z is O(1) (−0.108 at L=3,
  −0.089 at L=4, ~20% of tail range), not a small correction, and F is
  provably not self-complementing (M(p) + M(1−p) ≠ 0 exactly, §4.5).
* **W5 — resolved, not killed.**  The two actions do not intertwine
  (proved, §8.3): "the invariant shape" was ambiguous as phrased.  But
  the model's own involution structure gives a canonical gauge (§8.4)
  in which the question is well-posed and the pipeline's chart is
  demoted to a chart.  The probe's suspicion that "choosing the
  pipeline's action is not innocent" is confirmed and made precise.

### Main results (in the probe's own priority order)

1. **Theorem (W1-form): location theorems do not pin shape** — with the
   extra input percolation must supply named (§3).  Any reader who
   wanted `g` to be "the" correction must now find the shape input
   elsewhere.
2. **The canonical gauge and the no-intertwiner theorem** (§8): the
   shape question is well-posed exactly in the model's own gauge; the
   pipeline's Aff(1) is a chart; `g` is gauge-variant by construction.
3. **Strict separation of g from the shape tangent on the exact lab**
   (§6), plus **non-rank-1 batch covariance of the N=725 shape** (§7):
   the two objects are different projects.
4. **The exact laboratory** (§4–§5): per-rank-pair decomposition, Z on
   site L=3,4 and bond L=3, the odd-sector identity and its converse,
   the O(1) even part, the bond dual_fail finding, and the canonical
   reading `u_c = (1 + M(1/2))/2` of −21/64.

### What this probe did not do (per standing)

No STATUS.md edit; no closure of any other issue; no production; no
second exponent; 1.55 stays retired; #275 stays closed; no λ-sweep; no
N=725 rescoring beyond reading the committed block once; no literature
retrieval (one BLOCKED_ON_#620 marker, no search); no Cardy citation as
finite-L shape.

### Handoff sentences for the sister tickets

* To #617 (groupoid): the Z object is the (A)-quotient with canonical
  gauge; the decile grid is a (B)-chart; the groupoid of charts acts on
  the *recording* of Z, not on Z.
* To #618: the union-bound silence on rates is consistent with §3 —
  rates and shape are both outside Theorem L, at different orders.
* To #619: the exact polynomials for M_L can now include the
  decomposition (P11, P20, P02); the bond census's dual_fail is a new
  exact fact for their ledger.
* To #620: the single retrieval-shaped question this probe generated:
  *is there a theorem extracting the scaling limit of an inverse-CDF
  (quantile) shape from torus crossing probabilities?*  (Cardy gives
  crossing probabilities, not the rank law's inverse.)
* To #621: the sentences "g is the shape correction" and "M(1/2)
  governs the shape" both die here with objects attached.
