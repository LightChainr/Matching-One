# Probe MID (#624) — location without shape: explicit toy CDF families

2026-09-07. Analysis + numpy. No percolation engine, no arXiv, no N=725,
no #612 rescoring, no fit to #582 (1.55 is retired and irrelevant here).
Does **not** enter `docs/STATUS.md`. Closes nothing.

Brief: `notes/probes/probe-MID-location-without-shape-toys-20260907.md`
(PR #616). Theorem L is cited, not reproved: `notes/p613-quantile-convergence-20260907.md`
(PR #614) — uniform convergence of `Q_N(u)` to `p_c` on compacts of `(0,1)`;
location, not rate, not shape.

Verdicts in one line:

```text
T2  ACHIEVED (the prize). Z_N -> prescribed zeta; in these toys Z_N = zeta
    EXACTLY at every N.
T1  ACHIEVED. Parity switch (two accumulation shapes) and t_N = (1+sin log N)/2
    (dense oscillation, continuum of accumulation shapes), location holding.
T3  NULL. Anchors are gauge: the value-level affine identity between anchorings
    of one profile is exact (residual 0.0); no reparametrisation-free pathology.
D1  recorded as the non-example: fixed profile + shrinking width pins Z.
W5  numerical companion: Z quotients the Q-action; a p-axis kink changes the
    limit shape by O(1) at every N. The two quotients disagree.
```

Consequence for #622 W1: **a location theorem cannot pin the shape.**
Any conclusion of the form `Q_N -> p_*` uniformly on compacts is compatible
with every prescribed `Z_infinity`, with nonconvergent `Z_N`, and with any
anchor gauge — the full analytic laboratory is here.

---

## 1. Setup, and the two Aff(1) actions

Anchors `a = 0.2 < b = 0.8`, `p_* in (0,1)` (toys use `p_* = 1/2`),
`u`-grid `0.05(0.05)0.95`, `N = 2^k`, `k = 4..17` (four decades of `w_N`).

```text
on p:   T(p) = alpha p + beta       acts on laws: F -> F o T^{-1},  Q -> T o Q
on Q:   Q(u) -> alpha Q(u) + beta   (same formula, read on the quantile axis)
```

`Z_N(u) = (Q_N(u) - Q_N(a)) / (Q_N(b) - Q_N(a))` is invariant under the
**Q-action**: `(alpha Q + beta)` cancels exactly in numerator and denominator.
So as written, **`Z` quotients the Q-action**. For affine `T` the two actions
coincide on `(0,1)`-valued quantile functions only after the endpoint
renormalisation `(F - F(0)) / (F(1) - F(0))`, which is NOT affine in `F` — and
a nonlinear p-axis warp is not quotiented by `Z` at all (§7, the W5 companion;
the theorem is #622 W5's).

## 2. Master construction (Route A) and the location proof

Let `psi` be continuous strictly increasing on `(0,1)` (the window profile).
For `w_N > 0` small enough that the window fits in `(0,1)` define the
quantile function

```text
Q_N(u) = p_* + w_N psi(u),    u in (0,1),     Q_N(0) := 0,  Q_N(1) := 1,
```

and let `F_N := Q_N^{-1}` (the generalised inverse is a genuine continuous
strictly increasing bijection `[0,1] -> [0,1]` once `Q_N` is monotone with the
endpoint pins; cadlag is never needed — every `F_N` here is continuous and
strictly increasing on a neighbourhood of `p_*`).

**Location (uniform, on every compact — no rate claimed).**
`|Q_N(u) - p_*| = w_N |psi(u)|` for all interior `u`. On `[eps, 1-eps]` the
profile is bounded, `m_eps = sup_{u in [eps,1-eps]} |psi(u)| < oo`, so

```text
sup_{u in [eps,1-eps]} |Q_N(u) - p_*|  <=  w_N m_eps  ->  0.
```

The whole statement is the window: all interior quantiles sit in
`[p_* - w_N|min psi|, p_* + w_N|max psi|`, a window of width `O(w_N) -> 0`
about `p_*`. The endpoint convention moves only `u in {0,1}`; it is exactly
Theorem L's shape — location with a free profile.

**Constraint compliance** (the issue's fairness list): each `F_N` is a CDF on
`[0,1]` with `F_N(0)=0, F_N(1)=1`, nondecreasing (here strictly increasing and
continuous); strictly increasing near `p_*`; closed form throughout; not the
"mass parked at 0 and 1" cheat — the window `w_N -> 0` carries all interior
quantiles (mass of `[p_*-w, p_*+w]` is `1 - 2 u_tr -> 1` under the truncated
profiles, `= 1` for bounded `zeta` on the committed grid).

**Shape.** If `psi` is normalised at the anchors, `psi(a) = 0`, `psi(b) = 1`,
then

```text
Z_N(u) = (Q_N(u) - Q_N(a)) / (Q_N(b) - Q_N(a))
       = (psi(u) - psi(a)) / (psi(b) - psi(a))  =  psi(u)     for EVERY N.
```

The `w_N` cancels. The shape of the limit is the profile, not the location,
and not the width. That single line is the whole probe; everything below is
furniture.

## 3. T2 — the prize, stated and proved

**Proposition (prescribed shapes).** Let `zeta` be continuous strictly
increasing with `zeta(a) = 0`, `zeta(b) = 1`.

(a) If `zeta` is bounded on `(0,1)`, the family of §2 with `psi = zeta` is a
CDF family with `w_N -> 0` (any exponent), `sup_u |Q_N - p_*| <= w_N ||zeta||_oo
-> 0`, and

```text
Z_N(u) = zeta(u)   for every N and every u in (0,1)      (exact).
```

(b) If `zeta` is unbounded at one or both endpoints (logistic, Gumbel/GEV
tails), take the anchor-normalised truncations `zeta_tr` of §2 with
`psi = zeta_tr`: `zeta_tr = zeta` on `[tr, 1-tr]`, continuous strictly
increasing, flat extension outside, and `Z_N = zeta` exactly on the grid
whenever `tr <= u_min` (we use `tr = 0.025 < 0.05`); `Z_N -> zeta` locally
uniformly as `tr -> 0`.

*Proof.* §2 with `psi` replaced by the normalised/truncated profile; the
identity `Z_N = psi` is the anchor-normalisation cancelling `w_N`. Boundedness
of `zeta` on compacts gives uniform location. ∎

Surjectivity is by construction: the map `psi -> Z[psi]` from strictly
increasing profiles (normalised at the anchors) onto shapes is the identity.
**There is no "the" shape of a sharp threshold; there is a profile, and the
profile is free data.** If one insists on a one-parameter location-scale
window, see §5: the parameter must be promoted to a function, which is the
two-parameter ("width + profile") licence the issue grants.

Three committed instances (`results/probe-location-without-shape/latest.json`):

| family | profile | window | `max_N,u |Z_N - zeta|` |
|---|---|---|---|
| `t2_logistic` | `(logit u - logit a)/(logit b - logit a)` (symmetric) | `w_N = 1/N` | 4.6e-12 |
| `t2_gumbel` | `(-log -log u - ...)/(...)` (strongly skew, left tail) | `w_N = N^{-3/4}` | 5.5e-13 |
| `t2_kink` | piecewise-linear, kink at `u=0.6`, slope ratio 3:1 | `w_N = N^{-3/4}` | 6.5e-13 |

The `shape_err` is float cancellation, not convergence: the identity is exact
algebra. Location envelope `|Q_N - p_*| <= 3 w_N` holds at every committed `N`
(`verification.*.location_envelope`).

## 4. T1 — oscillation while location holds

**T1a (parity).** `psi_N = zeta_logistic` for even `k`, `zeta_kink` for odd
`k` (`N = 2^k`; the switch is on parity of the log, since `2^k` itself is
even — same mechanism for arbitrary integer `N`). Both profiles are
anchor-normalised, so `Z_N = psi_N` exactly and

```text
lim_k Z_{2^{2k}} = zeta_logistic,     lim_k Z_{2^{2k+1}} = zeta_kink,
```

two distinct accumulation shapes; the gap at the last committed pair is
`0.412` in sup norm (`t1_parity_switch.shape_gap_persists`). Location:
`|Q_N - p_*| <= w_N max(||zeta_l||, ||zeta_k||)` → the same envelope check
passes (`envelope_3w: true`). "Sharp ⇒ shape converges" is dead.

**T1b (dense oscillation).** One family, blend
`psi_N = (1 - t_N) zeta_logistic + t_N zeta_kink`, `t_N = (1 + sin(log N))/2`.
Convex combinations of anchor-normalised strictly increasing profiles are
again anchor-normalised and strictly increasing, so `Z_N = psi_N` exactly.
`t_N` is quasi-periodic in `k` (`log 2 / pi` irrational) and dense in `[0,1]`:
the accumulation set of `Z_N` is the **continuum** `{(1-t)zeta_l + t zeta_k :
t in [0,1]}`, and the pairwise spread over the committed grid is `0.203`
(`t1_sin_log.nonconvergence`). The issue's alternative `gamma_N = sin(log N)`
skew parameter is the same mechanism through a skewness coordinate.

## 5. D1 — the non-example (fixed profile pins the shape)

`F_N(p) = sigma((p - p_*)/w_N)` renormalised on `[0,1]`, `sigma` FIXED
(logistic; probit behaves identically), `w_N = 1/N`. Then `Q_N = p_* +
w_N sigma^{-1}(renorm(u))`, so `Z_N = Z[sigma^{-1}]` up to the renormalisation
drift `|sigma(±c/w_N)| <= e^{-c/w_N}` — constant in `N` once `w_N` is small
(`d1_window.shape_constant`: tail drift over `N >= 256` is `3.4e-12`; the
first-decade spread `2.12` is the `u=0.05` row sitting outside the `N=16`
window, a convention effect, not shape motion).

**Reading.** This is the family people silently imagine when they say "sharp
threshold ⇒ universal shape". It is the **non-example**: location PLUS a fixed
window profile pins `Z`. Theorem L gives no profile. Within one-parameter
location-scale windows the limit shape is determined by `sigma` and the window
cannot oscillate without oscillating `sigma` — T2 needs the profile promoted
to a free function (our Route A), which is precisely the freedom Theorem L
permits and the D1 imagination denies.

## 6. T3 / D4 — the anchors are gauge (nullish pathology)

One profile (`zeta_gumbel`), three anchor pairs `(0.1,0.9)`, `(0.2,0.8)`,
`(0.3,0.7)`. Write `p_a = psi(a)` etc. Then the three `Z_infinity` obey the
EXACT value-level affine identity

```text
Z_{a',b'}(u) = ( Z_{a,b}(u) * (p_b - p_a) + (p_a - p_{a'}) ) / (p_{b'} - p_{a'}),
```

residual `0.0` in float at every grid point (`d4_anchor_orbit.max_affine_
residual = 0.0`). Equivalently, as functions of `u`:
`Z_{a',b'} = Z_{a,b} o m` with `m = psi^{-1} o Aff_psi o psi` an increasing
map, verified on the grid wherever `m` lands inside `(0,1)` (residual `0.0`
for the widening `(0.1,0.9)`). Two honest caveats, recorded in the JSON:

1. `m` **fixes `{a,b}` pointwise only when the affine renormalisation is the
   identity** (`m_maps_default_anchors` false for genuine anchor changes,
   `u_map_fixes_anchors: false`). Under the issue's strict reading ("one orbit
   under increasing maps that fix `{a,b}`"), the curves are trivially NOT on
   one anchor-fixing orbit — at `u = a` one curve is `0` and the other is
   `Z_{a',b'}(a) != 0`. That failure is pure relabelling: the curves differ by
   the exact affine map above, i.e. by the anchor convention itself.
2. For anchor intervals strictly inside `(a,b)` (here `(0.3,0.7)`), the orbit
   reparam `m` is defined only on the band of `u` whose stretched profile
   stays in the profile range; off the band the committed curve is the affine
   extrapolation of a different restriction of `psi`. This is an extrapolation
   domain, not a shape disagreement: value-level the affine identity is exact
   everywhere (residual `0.0`).

**Verdict.** "Anchors are gauge" is the correct sentence for this family:
the invariant content of `Z` is the Aff-class (in profile values) of the
window shape, independent of anchor choice; the numerical curve is anchor-
affine by an exact identity. `Z_N` converging to anchor-dependent functions
that are NOT related by increasing reparametrisation does not occur — T3 not
exhibited, and we see no Route-A mechanism to exhibit it: the anchors enter
only through the normalisation of one and the same `psi`.

## 7. The two quotients disagree (numerical companion of #622 W5)

- **Q-action**: `Q -> alpha Q + beta`, `alpha = 0.9`, `beta = 0.12`:
  `Z` unchanged, sup deviation `1.6e-14` (`warp.q_action_invariance`) —
  exact invariance, as §1 says it must be.
- **p-axis warp with a kink inside the window** (the nonlinear part is the
  point; affine warps act identically on both axes): `T` piecewise-linear with
  slopes `1` / `1.5` about the window point `p_* + w zeta(1/2)`. The warped
  quantiles `T o Q_N` generate `Z_warped` with sup deviation `0.123` from the
  unwarped `Z` — at **every** `N`, exactly `w`-independent (§7 mechanism:
  piecewise scaling of the profile does not cancel). The two quotients are
  different equivalence relations on the same data.
- A smooth warp (`T' = 1` at `p_*`) would agree in the limit and differ by
  `O(w_N)` at finite `N`; the kink makes even the limit shape move. Theorem
  status belongs to #622 W5; this is the witness.

## 8. D5 — what extra input would pin `zeta` in percolation

Named inputs already in the #613/#606 notes. "pin `sigma`" = pin the window
profile (shape); "width" = the exponent/size of the window. Unknown is allowed;
nothing here is quoted from new papers.

| input | pins profile `sigma`? | pins width only? | note |
|---|---|---|---|
| RSW (two-sided power-law arm bounds) | **N** | **N** (bounds, not rate) | exponent statements about probabilities at fixed `p`; the union bound of #613 already gives no width rate (#618). A fortiori no profile. |
| four-arm / F1 (#321 ledger; open on square-site) | **unknown** | only if it comes with a rate | if F1 produced `|Q_N - p_c| <= C N^{-alpha}` it would pin the WIDTH; a profile needs more: convergence of the rescaled process in the window. |
| self-duality (`p_c = 1/2`, exact) | **partial** | N | forces the profile equivariance `psi(u) + psi(1-u) = 1` (equivalently `Z` symmetric about `(1/2, 1/2)`) — one constraint, infinitely many degrees of freedom left. Location is free here, shape is exactly what survives unconstrained. |
| full scaling-limit / conformal invariance (LSW-type, assumed for the laboratory) | **Y** (by assumption) | Y (gives `nu`) | this is the D1 non-example imported as physics: assuming the window process converges is assuming a fixed profile. Available triangular, open square-site (F1). |
| exact finite-L enumeration (L=3,4) | Y at those L | N | two points are not a sequence; kills constancy, not the profile-freedom. |

The honest summary: the repository's existing named inputs either pin
nothing (RSW, F1 as used), pin a single symmetry constraint (self-duality),
or *assume* what would have to be proved (scaling-limit ⇒ fixed profile).
A profile-pinning input for square-site does not currently exist in-tree.

## 9. Interface and claims

- **#622 W1: promoted (for the analytic laboratory).** "Location theorems
  do not pin shapes" now has a proof for the toys: Prop §3 + §2. The transfer
  to percolation is a quantifier move — Theorem L's conclusion
  (`Q_N -> p_c` uniform on compacts) is satisfied by the §2 family verbatim
  with any `zeta` (and by T1's with no `Z_infinity` at all), so no theorem
  with Theorem L's conclusion strength can pin the percolation shape. The
  counterweight W1 asked for is delivered; the definition of percolation `S`
  remains #622's own job.
- **#618:** may cite Prop §3 as "location does not imply a common `omega` for
  the 9-vector" — the T2 machinery is proved here (the exact identity
  `Z_N = psi`, with the truncation clause for unbounded profiles).
- **#620:** not involved (no retrieval performed). **#613:** cited, not
  reproved. No STATUS entry; issue left open.

## 10. Files

```text
notes/probe-location-without-shape-20260907.md      (this note)
scripts/probe/toy_cdf_families.py                   (numpy; reruns all tables)
results/probe-location-without-shape/latest.json    (Q_N, Z_N per family; 16/16 checks pass)
```

Run: `/Users/lc/.workbuddy/binaries/python/envs/default/bin/python
scripts/probe/toy_cdf_families.py` (numpy only). Grid `u = 0.05(0.05)0.95`,
`N = 2^4 .. 2^17`.
