# Probe re-run — Invariant shape of the threshold law (#622), 2026-09-08

Frontier `claude/matching-one-workspace-pwr5pv` @ `8b5f9d1a`; base `main`
(at `07350444`). This is the second run of the #622 probe. The first run
(merged as PR #628, with the W4 kill in #629 and the direction-2 toy
theorem in #634) defined `S`, adjudicated W1–W5, and separated `g`. This
run **re-adjudicates those verdicts on the repaired inputs** — PR #653
(repaired bond laboratory, `dual_fail = 0`, exact `M(1/2) = 0`) and PR
#655 (corrected spin-0 N725 Q/Z, pooled histogram, correct `cos 4θ`
weights) — and reports what actually changes. Nothing is re-derived that
is cited; nothing entered `docs/STATUS.md`; no issue is closed; no
production, no MC, no exponent fit; N=725 is read as the committed block,
not scored; 1.55 stays retired; no literature search.

**Verdict in one line: the repaired bond laboratory converts W2 from
"half-killed" to "split precisely" — duality pins `Z` to the antisymmetric
class (proved, T1) but pins no member of that class (proved, T2) — and the
four-size decomposition (T3) locates the surviving shape-limit question in
the antisymmetric part `A_N` of `Z`, under a now-exactly-stated extra
hypothesis C3. The north-star verdict stays (C), but the hypothesis is
sharper than the one #628 named.**

---

## 0. Inputs (cited, imported, not redone)

| Object | Source | Status |
|---|---|---|
| Exact site L=3,4 census, rank-pair decomposition | #606 + merged #628 census | cited |
| Repaired bond L=3 census (`dual_fail=0`, pairs 75460/111224/75460, `M(1/2)=0`, `F(1/2)=1/2`) | **PR #653**, `results/probe-invariant-shape/census-exact.json` | imported; replaces the #628 bond numbers, whose `dual_fail=118133` was withdrawn as a bug compound (#632) |
| Corrected spin-0 N725 pooled Q/Z, 11 levels, delete-one SEs | **PR #655**, `results/probe-invariant-shape/n725-zflow-corrected.json` | imported; `n725-zflow.json` (#628) stays withdrawn on disk |
| Toy location/shape separation theorem | #628 `toy-families.json` | cited |
| Consensus 9-vector `g` | #605 lineage, `type582-residual/latest.json` | cited |
| Theorem L statement + union-bound silence on rates | #613/#614, #618 | cited |

The two #633 retractions that constrain this run are honoured exactly:
(1) batch-covariance rank does **not** bound the mean shape's dimension;
(2) the g-vs-ΔZ affine-removed angle is reported **unoriented (20.4°)**,
not "159.6° ≈ orthogonal". No sentence below revives either.

## 1. New machine-checked statements

Script: `scripts/probe_invariant_shape/rerun_shape_limit_20260908.py` →
`results/probe-invariant-shape/rerun-20260908.json`.
Tests: `tests/test_rerun_shape_limit_20260908.py` (all pass; the pre-existing
`tests/test_probe_invariant_shape.py` is updated to the repaired census
numbers, as #653 already did on its branch).

### T1 (theorem). Duality symmetry of the bond laboratory

If every configuration satisfies `r_b + r_w = 2` under the geometric dual
transport (i.e. `dual_fail = 0`), then, writing `P11, P20, P02` for the
per-rank-pair count polynomials in the number of occupied bonds `k`:

```text
P02[k] = P20[NB−k],   P11[k] = P11[NB−k]
⇒  M(p) = −M(1−p)  (M antisymmetric about 1/2)
⇒  F(p) = 1 − F(1−p)  (self-complementing; F(1/2) = 1/2)
⇒  with the symmetric anchors a = 1/4, b = 3/4:
   Z(u) + Z(1−u) = 1 and Z(1/2) = 1/2 exactly.
```

Proof: complementing the occupied set swaps `r_b ↔ r_w` (vacant primal =
occupied dual under transport), giving the pair-count symmetry; the two
function identities are the reindexed binomial sums (mass
`P20+P11+P02 = C(NB,·)` is exactly the `dual_fail = 0` condition); the `Z`
identities are the affine chart applied to `Q(u) = 1 − Q(1−u)`.

Machine-checked on the repaired census: `Z_antisymmetry_max_err ≈ 1e-15`
(float side of an exact object), `M(1/2) = 0`, `Z(1/2) = 1/2`.

**This reverses #628's §4.3.** The "bond laboratory loses the duality
identity" finding was the broken sign, not physics. The site `r_b + r_w = 2`
identity *does* have a bond carrier, exactly.

### T2 (theorem, by construction). Duality pins the class, not the shape

`F` self-complementing ⟺ density symmetric about 1/2. So the T1 class is
exactly the symmetric-density class — and inside it `Z` is still free:
two members, `f ∝ e^{−(p−1/2)²/2s²}` and the same with a `(1 + 3·tanh²((p−1/2)/s))`
tail-weighting, both satisfy every T1 constraint to integration error and
have shapes differing by **0.1503** at `u = 0.9` (both exactly
antisymmetric). Recorded in the JSON with both `Z` vectors.

**Consequence for W2 (final form).** W2 split in two, both halves now
exact: (i) "self-duality constrains the shape" is *true but only to the
antisymmetric class* — one linear constraint per reflected pair, no
member pinned; (ii) "the only finite-size object is shape" survives
intact. W2's hoped-for mechanism ("self-duality pins the shape") is dead
as a pinning mechanism; its `p_c = 1/2` freeness was always trivial.

### T3 (observation, four sizes, no fit). The A/S decomposition

```text
A(u) = (Z(u) − Z(1−u) + 1)/2   (antisymmetric part: A(u) + A(1−u) = 1)
S(u) = (Z(u) + Z(1−u) − 1)/2   (symmetric part: zero iff T1 holds)
Z = A + S.
```

Across site L=3, site L=4, bond L=3 (self-dual), and the corrected spin-0
N725 block, on the decile grid, `u = 0.1 … 0.5`:

```text
S(0.1):  site3 −0.0541 → site4 −0.0446 → bond3 0 (exact) → N725 −0.0090
A(0.1):  site3 −0.4116 → site4 −0.4279 → bond3 −0.4352 → N725 −0.4591
```

Pairwise max spread over all four objects: **A: 0.0475**, S: 0.0541. The
ordering at `u = 0.1` is monotone in the sequence (site L3, site L4, bond
L3, site N725) — recorded as a trend, explicitly **not** claimed as
convergence, and no rate is fitted (two exact sizes + one self-dual size +
one committed block do not make a limit).

Reading: the *symmetric* part of the shape — the sector #628 measured as
the "even deformation" — collapses toward 0 as the objects gain size or
duality; the *antisymmetric* part is the slowly-drifting survivor. The
shape-limit question, in this chart, therefore lives in `A_N` alone.

### T4 (conjecture, stated to be killable). C3 — self-matching restoration

```text
C3:  sup_{p ∈ K} |M_N(p) + M_N(1−p)| → 0  for every compact K ⊂ (0,1).
```

Equivalent at the level of the grid to `S_N → 0` (both sides exact at
finite N via `F = (1+M)/2`). Under C3, any shape limit is antisymmetric
and the open object is `A_∞ = lim A_N`. Status: **open** — not implied by
Theorem L (whose hypotheses are spent on location), not proved here, not
refuted here. Its finite-N witnesses are exact: `M(1/2) = −21/64` (site
L=3), `−13757/32768` (site L=4), `0` (bond L=3), and the corrected N725
block's `S(0.1) = −0.0090 ± 2.4e-5`.

## 2. Re-adjudicated verdicts (deltas from #628 only)

| Verdict | #628 | This run | Why |
|---|---|---|---|
| North-star (3) | (C), extra input "RSW-type shape-side access" | **(C), sharpened**: the concrete, rank-moment-level hypothesis is C3 (T4); the shape-side input #628 named is sufficient-not-necessary language for the same gap | T1–T3 |
| W1 | killed as question, promoted as theorem | stands unchanged | #628 §3 + #634; T2 is the duality-class refinement: even *within* the antisymmetric class, location + duality leave the shape free |
| W2 | half-killed ("duality constraint does not transfer to bonds") | **split precisely**: the transfer *does* happen (T1); what fails is pinning (T2) | PR #653's repair removed the false negative |
| W3 | killed (94°/160° + non-rank-1 covariance) | **downgraded to "undetermined at lab sizes"** per #633: report the unoriented 20.4°, and the covariance argument is retracted; the separation sentence is withdrawn to "no lift exhibited; the lab angle is 20.4° unoriented, which does not separate" | #633 retractions 1–2 |
| W4 | killed both ways | stands, with the §5 counterexample *downgraded* per #632: `P11 ↦ 3/2 P11` leaves the simplex (mass 593/512), so it is not a probability counterexample; the surviving exact sentence is the weaker "equal M does not determine the full joint `(P11, P20, P02)`" — the W4 kill on "M(1/2) governs" (pointwise value vs functional of the full polynomial) is untouched | #632 table |
| W5 | resolved (canonical gauge) | stands unchanged | no new input |

## 3. The g-lift question after the corrections

The corrected N725 JSON records the only admissible comparisons:

* tiny-torus `ΔZ` vs `g`: raw 93.9°, affine-removed **unoriented 20.4°** —
  not a separation;
* `DZ_Q[g]` vs `Z_725` in the 0.2/0.8 chart: unoriented 81.4° — not the
  g-vs-ΔZ test (no second committed size), not used as a claim;
* N725 g-vs-ΔZ: **undetermined** without another committed size. Stopping
  there, per #633.

So the merged "g is separated from the shape tangent" headline of #628 is
**withdrawn to**: *no lift with a transformation law has been exhibited;
the one lab-size angle compatible with the repaired inputs is 20.4°
unoriented, which neither supports nor refutes a lift.* The main result of
this run does not depend on that question either way.

## 4. What this run contributes, in priority order

1. **T1**: the repaired bond laboratory is an exactly self-complementing
   census — `M` antisymmetric, `F` self-complementing, `Z` antisymmetric,
   all proved from `dual_fail = 0` and machine-checked. This is the clean
   laboratory the issue asked for, restored.
2. **T2**: duality pins `Z` to the antisymmetric class and nothing more
   (spread 0.15 inside the class). W2's mechanism dies a *precise* death.
3. **T3 + T4**: the shape-limit question is localised to `A_N` under the
   explicit conjecture C3 (`M_N + M_N(1−·) → 0`), with exact finite-N
   witnesses at four sizes and no fitted rate. This is the concrete form
   in which the next probe (or a production ticket, if one is ever
   funded) should attack the limit: measure `S_N`'s decay and `A_N`'s
   drift on the committed blocks that already exist.

## 5. What was *not* done (per standing)

No `docs/STATUS.md` edit; no issue closed; no merge; no production or MC;
no second scale (1.55 retired); no N=725 scoring — the block is read once,
through PR #655's corrected object; no λ-sweep, no jackknife invention; no
literature retrieval (no `BLOCKED_ON_#620` marker was needed: every
theorem used is elementary combinatorics on the census itself); no Cardy
citation as finite-L shape; GOVERNANCE §2E honoured (site and bond remain
two experiments; the T3 comparison is cross-model and labelled as such).

## Handoff

* To #619: the bond pair-count symmetry `P02[k] = P20[NB−k]`,
  `P11[k] = P11[NB−k]` is an exact ledger fact with `dual_fail = 0`.
* To #620: unchanged single question — is there a theorem from torus
  crossing probabilities to the inverse-CDF shape? Now with a sharper
  target: the *antisymmetric part* `A_∞`, and the testable conjecture C3.
* To any future shape probe: `S_N`'s decay is measurable on existing
  committed histograms (`S(0.1) = −0.0090 ± 2.4e-5` at N=725); `A_N`'s
  drift across L=3,4,N=725 spans only 0.048 — both are cheap re-reads, and
  the next size is the only thing that can turn the T3 trend into a
  verdict.
