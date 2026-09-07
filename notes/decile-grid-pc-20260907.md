# #641 — Does the whole decile grid pin p_c?

Date: 2026-09-07. Branch: `claude/641-decile-grid-pc-pairwise` (off `claude/matching-one-workspace-pwr5pv`).

## Verdict

`E2_REACHED_AT_N_20_HNF_2_1_10_VS_HNF_2_0_10_DECILE_GAP_2.248e-6_LT_3e-6_ROOTS_DISJOINT;_E4_IMPOSSIBLE_FOR_THE_WHOLE_CLASS;_KNIFE_EDGE_ON_ROOT_VS_QUOTED_PRECISION`

Two admissible finite models share the nine-decile quantile vector to better
than the published SE while having different `p_c`: HNF(2,1,10) vs
HNF(2,0,10), N = 20, max decile bracket gap exactly
648033958315/288230376151711744 = 2.2483e-6 < 3e-6, with disjoint p_c
brackets. The deliverable is the ladder disposition on the admissible class
of finite periodic square-lattice quotients with N ≤ 20 (104 D4-inequivalent
geometries swept at N ≤ 18, plus the N = 20 pair from the closest family),
plus one proof that covers the entire class at any N.

**Self-check correction (2026-09-07, second revision).** The first commit of
this note claimed `E2 NOT REACHED` and "gap ~ 1/N, 3e-6 needs N ≳ 50". Both
claims were wrong and were caught by the post-commit self-check: extending
the (2,1,L)/(2,0,L) family one rung past the sweep cutoff crosses the 3e-6
SE scale at N = 20, and the family's gap decay is ~N^−9 to N^−10 locally
(roughly one order of magnitude per rung), not 1/N.
`scripts/decile_grid_pc_selfcheck.py` is the checking suite that caught
this; `results/decile-grid-pc/e2-n20-pair.json` is the exact artifact for
the E2 pair.

**Self-check correction (2026-09-07, third revision).** A second pass caught
three more defects: (a) the first two commits called HNF(2,1,L) "the 3xL
torus" — wrong, its cell is 2×L with a helical vertical wrap (winding
(1,L)); (b) the second commit projected a "decisive" E2 instance at
N ≈ 25–30 in this family — backwards, since decile gap = root gap exactly
in the family, larger N shrinks both and moves AWAY from the
root-beyond-precision criterion; no family member satisfies both E2
criteria at once (see Family scaling); (c) "seven equal-M(1/2) groups,
all with distinct roots" — the seventh group (eighteen (1,0,L) chains,
M(1/2) = 0) has p_c = 1/2 for every member, so E3's distinct-roots clause
fails there; the artifact's six nonzero groups are the ones where the
clause holds.

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

E2  REACHED at N = 20: HNF(2,1,10) (the 2x10 torus with a helical vertical
    wrap, winding vector (1,10)) vs HNF(2,0,10) (the plain 2x10 torus)
    have max decile gap
        648033958315 / 288230376151711744  =  2.2483e-6  <  3e-6,
    the SE scale of the published nine-decile vector, while their p_c
    brackets are disjoint (root gap also 2.2483e-6, exact). Caveat,
    stated honestly: on the stricter reading "different p_c by more than
    our quoted precision", the case is knife-edge -- the root gap
    2.248e-6 is of the same order as the "a few parts in 10^6" precision
    quoted for Q(0.5) in notes/wasserstein-shape-flow-20260906.md. E2 is
    reached on the decile-within-SE criterion; the root-vs-precision
    criterion is not cleanly separated, it is borderline. Sharper: within
    this family the decile gap and the root gap are EQUAL (exact, L=7..10;
    same order at L=5,6), so the two E2 criteria trade off against each
    other as L grows -- decile-within-SE needs N >= 20, root-beyond-
    quoted-precision needs N < 20. No member of this family satisfies
    both criteria at once. A clean E2 instance needs a pair whose decile
    and root separations differ in order of magnitude, not more rows of
    this one.

E3  NO instance: six groups of geometries with >1 member share M(1/2)
    exactly at nonzero values (e.g. -21/64 is shared by HNF(3,1,2),
    HNF(7,2,1), HNF(4,1,2), HNF(3,0,3); -5/16 by HNF(5,1,1), HNF(5,1,2),
    HNF(5,2,1)), and every pair inside every one of those six groups has
    BOTH distinct roots AND distinct decile vectors. A seventh group
    shares M(1/2) = 0 -- the eighteen (1,0,L) chain geometries, L=1..18 --
    but every member has p_c = 1/2 exactly, so E3's "different p_c" fails
    there regardless of deciles; results/decile-grid-pc/e3-groups.json
    records only the six groups where roots can differ. E3 needs decile
    agreement first; that already fails at N <= 18.
```

## The headline pair

HNF(2,1,10) vs HNF(2,0,10) — periods {(2,0),(1,10)} (the 2×10 torus with a
helical vertical wrap, shortest non-horizontal winding (±1,10), Manhattan
length 11; N=20) against {(2,0),(0,10)} (the plain 2×10 torus, vertical
winding 10; N=20), from the closest pair family:

- max decile gap, exact: 648033958315/288230376151711744 = 2.2483e-6 (< 3e-6 SE scale);
- min root gap: 2.2483e-6, exact — same order as the decile gap;
- independent check (governance §2): the headline gap was recomputed by a
  second route (raw Bernstein counts → separate CDF → Sturm, no shared
  decile-bracket state) and agrees exactly; Sturm and Fraction bisection
  agree to 12 digits on every decile of the pair
  (`scripts/decile_grid_pc_selfcheck.py`, C4/C5).

## Family scaling (measured, not fitted)

The (2,1,L) vs (2,0,L) torus family, exact max decile bracket gaps:

```text
L=5  (N=10)  8.18e-4
L=7  (N=14)  6.82e-5
L=8  (N=16)  2.16e-5
L=9  (N=18)  6.92e-6   <- prior sweep cutoff; read as "E2 not reached"
L=10 (N=20)  2.25e-6   <- E2 threshold 3e-6 crossed
```

The decay is roughly one order of magnitude per rung, ~N^−9 to N^−10 locally
(N=14→18 exponent ≈ 9.1, N=16→20 ≈ 10.1, N=18→20 ≈ 10.7) — NOT 1/N as the
first commit claimed. The first commit's "halves per unit L" and "3e-6
needs N ≳ 50" statements were arithmetic errors on these same measured
numbers (the L=9→L=10 step drops the gap by 3.1x, not 2x, and 2.25e-6 <
3e-6 already at N = 20). Consequence: **on the admissible class,
decile-vector separation and root separation are the same size and collapse
super-polynomially fast in N for this family** — the decile grid is a far
weaker probe of p_c than the first commit stated. And the trade-off cuts
sharper than "go to larger N": in this family the decile gap and the root
gap are exactly equal (L=7..10), so growing N shrinks BOTH together. The
decile-within-SE criterion is met at N ≥ 20, but the root-beyond-quoted-
precision criterion (root gap ≫ a few parts in 10^6) is met only at
N < 20 — no member of this family satisfies both at once. A decisive E2
instance requires a geometry pair whose decile and root separations have
different orders of magnitude; this family cannot supply one at any
enumerable N.

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

3. **The E2 pair at N = 20** (`results/decile-grid-pc/e2-n20-pair.json`,
   `scripts/decile_grid_pc_e2_n20_pair.py`): the (2,1,L)/(2,0,L) family at
   L = 5..10 computed exactly, the L = 10 pair crossing the 3e-6 decile-SE
   threshold, with an exact independent recomputation of the headline gap.

4. **E3 groups** (`results/decile-grid-pc/e3-groups.json`): all pairs
   sharing M(1/2) exactly with roots able to differ (six nonzero groups),
   with exact decile/root gaps. The seventh equal-M(1/2) group — the
   eighteen (1,0,L) chains, all with M(1/2) = 0 and p_c = 1/2 — is excluded
   there because E3's "different p_c" is impossible inside it; the
   exclusion is verified, not an oversight.

5. **Governance §2 checks**: (a) regression identity — the sweep's HNF
   builder reproduces the committed reference implementation's Bernstein
   counts exactly on all five committed geometries (axis L=2,3,4; diamond
   L=2,3); (b) independent checks of the headline numbers by different
   methods — Sturm vs Fraction bisection agree to 12 digits on all nine
   deciles and the root of the headline pair; the N = 20 gap is recomputed
   from raw Bernstein counts via a separate route and agrees exactly
   (`scripts/decile_grid_pc_selfcheck.py`, checks C1–C6 + 7–16, all green).
   M(0) = −1 and M(1) = 1 hold for every geometry, as the wrapping
   convention requires.

## What it means for the readout (one sentence, per the ticket)

Two admissible finite models at N = 20 share the nine-decile vector within
measurement noise while having different p_c — the decile grid does NOT
deterministically pin p_c on the admissible class (E2), though the sharpest
available instance is knife-edge relative to our quoted precision, so the
practical-motivation case for #635/#636's readout change rests on the
#615-style analytic blindness (identical response, thresholds 1/2 and 1/3),
not on a clean finite-geometry counterexample at enumerable sizes.

## Boundary

- E1 is a negative confined to the N ≤ 18 swept subclass (104 geometries);
  E2 is a positive at N = 20, one rung past the cutoff. The knife-edge
  character of the N = 20 pair is structural for this family: decile gap =
  root gap exactly (L=7..10), so the two E2 criteria cross without a
  simultaneous solution anywhere in the family. A *decisive* E2 instance
  needs a pair with decile and root separations of different orders, which
  this ticket does not supply and this family cannot.
- "Admissible" is taken exactly as the ticket defines it: integer Bernstein
  form arising from D(C) = 1{black wraps} − 1{white NN+NNN wraps} on some
  finite periodic geometry. No toy families, no stretching.
- Does not close #615, #624, #625, #628, or #641's relatives.

## Files

```
scripts/decile_grid_pc_pairwise.py            committed-geometry pairwise check (deliverable 3)
scripts/decile_grid_pc_admissible_sweep.py    HNF admissible-class sweep, E1-E4 ladder
scripts/decile_grid_pc_e2_n20_pair.py         N=20 E2 pair + family L=5..10 (exact)
scripts/decile_grid_pc_selfcheck.py           post-commit self-check suite C1-C6 + 7-16
results/decile-grid-pc/committed-pairwise.json    deliverable 3 artifact (exact)
results/decile-grid-pc/admissible-sweep.json      sweep artifact (exact pairs, E4 argument)
results/decile-grid-pc/closest-pair-exact.json    the 6.92e-6 pair, per-decile exact gaps
results/decile-grid-pc/e2-n20-pair.json           the E2 pair at N=20 + family scaling (exact)
results/decile-grid-pc/e3-groups.json             the six nonzero equal-M(1/2) groups (the
                                                  M(1/2)=0 chain group is excluded: same p_c)
```
