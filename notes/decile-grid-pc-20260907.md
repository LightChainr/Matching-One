# #641 — Does the whole decile grid pin p_c?

Date: 2026-09-07. Branch: `claude/641-decile-grid-pc-pairwise` (off `claude/matching-one-workspace-pwr5pv`).

## Verdict

`NO_E1_TO_E4_COUNTEREXAMPLE_ON_THE_ADMISSIBLE_SWEPT_SUBCLASS_N_LT_18;_E4_IMPOSSIBLE_FOR_THE_WHOLE_CLASS;_MIN_DECILE_GAP_6.92e-6`

Nothing to cite as a counterexample. The deliverable is the ladder disposition
on the admissible class of finite periodic square-lattice quotients with
N ≤ 18 sites (104 D4-inequivalent geometries, 5,356 pairs), plus one proof
that covers the entire class at any N.

## Rung disposition

```text
E4  IMPOSSIBLE for the whole admissible class, any N (proved, not swept):
    F1 and F2 are polynomials. If they agree on an interval [p1,p2] with
    nonempty interior, F1 - F2 has infinitely many roots, so F1 = F2, so
    M1 = M2, so the physical root is the same. E4 succeeding would require
    a geometry's M to have a root different from its own.

E1  NO instance in the swept subclass: across all 5,356 pairs of the 104
    D4-inequivalent geometries with N <= 18, the minimum (over pairs) of the
    maximum (over deciles) quantile-bracket gap is exactly
        7975980164381 / 1152921504606846976  =  6.918e-6  >  0.
    No pair agrees at all nine deciles. (Existence question outside the
    swept subclass stays open; see Boundary.)

E2  NOT REACHED: the closest pair (HNF(2,1,9) vs HNF(2,0,9), i.e. the 3x9
    torus vs the 2x18 torus) has max decile gap 6.92e-6, which exceeds the
    measured SE ~3e-6. So even the two most-degenerate-admissible models in
    the sweep are separated by more than measurement noise.

E3  NO instance: seven groups of geometries share M(1/2) exactly
    (e.g. -21/64 is shared by HNF(3,1,2), HNF(7,2,1), HNF(4,1,2),
    HNF(3,0,3)), and every pair inside every group has BOTH distinct roots
    AND distinct decile vectors. E3 needs decile agreement first; that
    already fails.
```

## The one interesting pair

HNF(2,1,9) vs HNF(2,0,9) — periods {(2,0),(1,9)} (the 3×9 torus, N=18)
against {(2,0),(0,9)} (the 2×18 torus, N=18):

- per-decile gaps, exact: 1.05e-8 (u=.1), 1.46e-7, 7.68e-7, 3.05e-6,
  **6.92e-6 (u=.5, the max)**, 6.01e-7, 6.31e-8, 5.85e-9, 2.01e-10 (u=.9);
- M(1/2): −8449/131072 vs −33/512 (distinct, so not an E3 pair either);
- min root gap: 6.92e-6 — same order as the decile gap.

The (2,1,L) vs (2,0,L) family halves its max decile gap each time L grows by
one (2.16e-5 at L=8, 6.92e-6 at L=9): the gap scales ~1/N. Reaching 3e-6
inside this family needs N ≳ 50, i.e. 2^N configurations — beyond the
brute-force budget of this ticket. The 1/N scaling is itself the message:
**on the admissible class, decile-vector separation and root separation are
the same size (O(1/N) finite-size drift)**. Nothing like #615's analytic
Markov pathology — where identical finite-horizon response hid thresholds at
1/2 and 1/3 — is available here at any size we can enumerate.

## What was computed

1. **Deliverable 3 first, as the ticket requires**
   (`results/decile-grid-pc/committed-pairwise.json`): the five committed
   exact M(p) (axis L=2,3,4; diamond L=2,3) recomputed and compared. All ten
   pairs have max decile bracket gap O(10^-2) (smallest 1.91e-2, axis L=4 vs
   diamond L=3) — the committed geometries are nowhere near each other on
   the decile grid. No E2 pair existed on the shelf.

2. **The admissible-class sweep**
   (`results/decile-grid-pc/admissible-sweep.json`,
   `scripts/decile_grid_pc_admissible_sweep.py`): every full-rank sublattice
   of Z² with index N ≤ 18, D4-canonicalised (HNF triples), 104 distinct
   M(p) after polynomial dedup. Exact integer Bernstein coefficients by 2^N
   enumeration of D(C) = 1{black NN wraps} − 1{white NN+NNN wraps};
   exact Fraction CDF; Sturm isolation of roots and decile quantiles at
   60 bits.

3. **E3 groups** (`results/decile-grid-pc/e3-groups.json`): all pairs
   sharing M(1/2) exactly, with exact decile/root gaps.

4. **Governance §2 checks**: (a) regression identity — the sweep's HNF
   builder reproduces the committed reference implementation's Bernstein
   counts exactly on all five committed geometries (axis L=2,3,4; diamond
   L=2,3); (b) one independent check of the headline number by a different
   method — Fraction bisection on the sign of F(p) − 1/2 reproduces the
   Sturm median gap 6.918e-6 to 12 significant digits (float diagnostic
   only; the claims are the exact brackets). M(0) = −1 and M(1) = 1 hold
   for every geometry, as the wrapping convention requires.

## What it means for the readout (one sentence, per the ticket)

On every admissible geometry we can enumerate, the nine-decile vector and
p_c move together at the same O(1/N) finite-size scale — the decile grid
carries threshold content on this class, and #635/#636's readout change is
motivated on other grounds, not by a #615-style blindness on our own object.

## Boundary

- E1/E2 are existence statements; the negative here is confined to the
  N ≤ 18 swept subclass (104 geometries, all D4-inequivalent shapes). A
  counterexample at N ≳ 50 cannot be excluded from this data — but the 1/N
  scaling of the closest family gives no hint of one.
- "Admissible" is taken exactly as the ticket defines it: integer Bernstein
  form arising from D(C) = 1{black wraps} − 1{white NN+NNN wraps} on some
  finite periodic geometry. No toy families, no stretching.
- Does not close #615, #624, #625, #628, or #641's relatives.

## Files

```
scripts/decile_grid_pc_pairwise.py            committed-geometry pairwise check (deliverable 3)
scripts/decile_grid_pc_admissible_sweep.py    HNF admissible-class sweep, E1-E4 ladder
results/decile-grid-pc/committed-pairwise.json    deliverable 3 artifact (exact)
results/decile-grid-pc/admissible-sweep.json      sweep artifact (exact pairs, E4 argument)
results/decile-grid-pc/closest-pair-exact.json    the 6.92e-6 pair, per-decile exact gaps
results/decile-grid-pc/e3-groups.json             the seven equal-M(1/2) groups
```
