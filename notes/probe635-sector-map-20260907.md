# probe635 — is M(p) the difference of two topological sector amplitudes?

Ticket: #635 (P1 long-horizon). Date: 2026-09-07. Branch
`research/probe635-sector-map-20260907`, base
`claude/matching-one-workspace-pwr5pv` @ 39097254.

Scripts: `scripts/probe635_sector_decomposition.py` (site),
`scripts/probe635_bond_sector_map.py` (bond wrap-gap),
`scripts/probe635_bond_rank_audit.py` (rank audit, three implementations).
Results: `results/probe635-sector-map/latest.json`,
`results/probe635-bond-rank-audit/latest.json`.

Everything below is exact integer/Fraction arithmetic; no floats anywhere in
either artifact.

## 1. Candidate decomposition (deliverable 1)

For a finite torus configuration C (site axis L=3/4 or diamond L=2, and the
bond analogue), define the winding label of a cluster by the homology class
of its displacement lattice, computed greedily on a spanning forest:

- `x` / `y`: the cluster winds exactly one axis (up to sign);
- `both-same`: a *single* cluster whose displacement lattice already spans
  H_1(T²; Q) with rank 2 (it wraps both axes simultaneously);
- `both-two`: two distinct clusters, one winding each axis;
- `none`: no winding.

`D(C) = 1{black NN wraps} − 1{white NN+NNN wraps}`. The candidate map is

    M(p) = Σ_k [ #(both-same | black, |C|=k) − #(both-same | white, |C|=k) ]
           · p^k (1−p)^{N−k}.

That is: the two sector amplitudes would be the black-winding and
white-winding indicators restricted to the `both-same` homology class.

## 2. Verdict (deliverable 2): **A, at finite L, but the sectors collapse**

The map above reproduces, bit-for-bit, the integer Bernstein coefficients of
the committed baseline on every geometry enumeration reaches:

- axis L=3 (2^9): `[-1,-9,-36,-78,-90,-36,36,36,9,1]` — exact;
- diamond L=2 (2^8): `[-1,-8,-28,-56,-42,8,20,8,1]` — exact;
- axis L=4 (2^16): `[-1,-16,-120,-560,-1812,-4272,-7448,-9424,-7874,-2896,1720,2832,1660,560,120,16,1]` — exact.

On all three geometries the per-label diagnostic tables show that **x, y and
both-two configurations never carry a nonzero D**: their contribution to M
is identically 0 at every k. Every unit of M comes from `both-same`
configurations. So the identity

    M(p) = A_black^{both-same}(p) − A_white^{both-same}(p)

holds **exactly at finite L** (verdict A in the issue's taxonomy), not merely
asymptotically.

## 3. The structural caveat (why this is not yet Jacobsen's two sectors)

At L=3/4 and diamond L=2, `both-two` is *empty as an event class*: on these
sizes, whenever a configuration wraps, a single cluster already does it
twice (or the wrap event never engages D). The "two topological sectors"
therefore degenerate, at reachable L, into one label — single-cluster
double-winding — measured twice (black on NN, white on NN+NNN). Whether the
both-two sector acquires nonzero weight at L ≥ 5 (beyond 2^16 exact
enumeration) is open and would decide whether the finite-L identity is a
degenerate special case or a structural law. This caveat is a statement
about evidence reach, not about the exactness of A at the sizes tested.

Independent check (GOVERNANCE §2): the sector tables were recomputed through
the separate wrap-label diagnostic path (union-find merge of wrap flags,
independent of the displacement-lattice rank path), and the sum
`p-both-same + m-both-same` equals the committed baseline coefficients
position-by-position; likewise `verify()` in the script asserts checksum
equality against `scripts/exact_matching_polynomial.py` output before any
sector table is printed.

## 4. Bond lab and the #628 rank defect (near-miss ruled out)

The bond analogue was probed with the wrap-gap observable
`G = 1{primal wraps} − 1{dual wraps}` (L=2, L=3; antisymmetric around 1/2,
G(1/2)=0 exactly). While reconciling with #628's bond census we found that
the published per-rank-pair bond numbers rest on an order-dependent rank
implementation.

**The bug.** `bond_ambient_rank` in `scripts/probe_invariant_shape/
exact_rank_census.py` (PR #628) closes each fundamental cycle as
`(ax, ay) + (dx, dy)` — it *adds* the stored edge displacement to the
tree-lift displacement. The cycle is tree path u→v followed by the stored
edge back v→u, so the correct winding is `(ax, ay) − (dx, dy)`. With the
wrong sign the computed "winding" depends on the chosen spanning tree, so
the function is not a function of the edge set: the L=3 mask-63 dual graph
scores 2 in ascending edge order and 1 in descending order (the latter,
and the independent linear-algebra computation below, agree with the
hand-checked topology: that graph winds only along x).

**Causal chain closed.** Re-running the verbatim #628 function over all
2^18 bond configurations reproduces the published numbers exactly:
`dual_fail = 118133`, pair counts (0,2):51704, (1,2):49263, (0,1):4397,
(1,1):41839, (1,0):5633, (2,2):12183, (2,1):46657, (2,0):50468. So those
published numbers are the direct output of the buggy function, and the
function changes rank on **25,090 of 262,144** configurations when the edge
iteration order is reversed — i.e. 25,090 of the published per-config rank
assignments are tree-order artifacts.

**Independent recomputation.** With the sign corrected
(`bond_ambient_rank_fixed`), certified against an independent
spanning-tree-free implementation (exact rational elimination on the
incidence matrix: rank[winding image] = rank[B; W] − rank[B]; 42,048
comparisons, 0 disagreements), the full L=3 bond census is:

    dual_fail = 115608   (still nonzero — the qualitative #628 finding survives)
    (0,0):4356  (0,1):26724  (0,2):44380
    (1,0):26724 (1,1):57776  (1,2):26724
    (2,0):44380 (2,1):26724 (2,2):4356

The corrected census satisfies the sanity laws the buggy one could not:
pair counts symmetric under (a,b) ↔ (b,a); M_bond = P20 − P02 antisymmetric
around k=9 with M(0)=−1, M(1)=+1, M(1/2)=0 exactly, coefficients

    [-1,-18,-153,-804,-2880,-7254,-12552,-13356,-6498, 0,
      6498,13356,12552,7254,2880,804,153,18,1].

(The earlier probe's wrap-gap polynomial `gap` differs from this rank-gap
M_bond from k=3 up — e.g. −810 vs −804 at k=3 — which is expected: wrap
*events* and rank-*pairs* are different observables; the discrepancy is
recorded here so nobody reads one as the other.)

**Impact statement for #628.** The qualitative claim "`r_b + r_w = 2` breaks
on the bond torus" survives the fix (115,608 failures), but every published
per-(r_b, r_w) count is wrong, wrong by tree order, and the correct table
has a different shape (nonzero (0,0)/(2,2) mass, all mixed pairs exactly
26,724). Any downstream use of #628's bond rank-pair table should be
recomputed from `results/probe635-bond-rank-audit/latest.json`.

## 5. Claim boundary

Verdict A is a finite-size exact statement on L ≤ 4 site (2^16) plus
diamond L=2; it is not a proof for all L, and the both-two emptiness at
reachable L is an observation, not a theorem. The bond rank defect is a
statement about PR #628's artifact `results/probe-invariant-shape/
census-exact.json` (bond section) and about the #622-inherited
`bond_ambient_rank` implementation; the site side of #628/#606 is *not*
implicated — the site `ambient_rank` closes its cycles with the correct
sign (`lift_step(w, v)`), which is why its census satisfies symmetry and
the bond one could not.
