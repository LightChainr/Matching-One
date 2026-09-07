# Results — probe-invariant-shape (issue #622)

Directory of the probe's committed artifacts.  Nothing here re-verifies a
pre-existing JSON; every file is a new object computed for #622.

## census-exact.json

`scripts/probe_invariant_shape/exact_rank_census.py`

Exact (Fraction) anchored shape Z on the site L=3 and L=4 honest tori —
importing PR #606's enumerator verbatim as a library — and the square-bond
L=3 census (2^18 configs) under the crossing-dual convention of
`scripts/square_bond_kappa3.py`, with geometric dual transport for the white
side.

Contents:

* per-rank-pair decomposition (P11, P20, P02) of the site censuses; the
  margins reproduce the published #606 pair counts (259/162/91 at L=3,
  36559/19932/9045 at L=4) — margins checked, numbers not re-derived;
* F = P20 + P11/2, M = P20 - P02; M(1/2) = -21/64 (L=3) and -13757/32768
  (L=4) come out of the same decomposition (cited facts, here consistent);
* Q at u = 1/4 and 3/4 by exact bisection (tol 1e-14), and Z on the frozen
  9-level grid, for site L=3, L=4 and bond L=3;
* Z odd part (Z(u) + Z(1-u) - 1) per grid point.

Key numbers (Z, anchors 1/4–3/4):

```text
site  L=3:  [-0.4657, -0.1308,  0.1179,  0.3295, 0.5232, 0.7102, 0.9002, 1.1061, 1.3575]
site  L=4:  [-0.4725, -0.1305,  0.1169,  0.3260, 0.5178, 0.7045, 0.8971, 1.1110, 1.3833]
bond  L=3:  [-0.7417, -0.1881,  0.1574,  0.4126, 0.6164, 0.7868, 0.9335, 1.0627, 1.1781]
```

The bond run reports `dual_fail = 118133` out of 262144.  That number is a
finding, not a bug, and the note analyses it: the identity r_b + r_w = 2
site-side comes from Alexander duality for complementary *vertex* sets in
complementary *connectivities* (4-conn vs 8-conn).  The bond side has no
such pairing with r_w defined as "rank of the occupied dual-bond graph":
duality for bond percolation swaps occupied primal bonds with *vacant* dual
bonds, and the rank of a vacant set is not the rank of the complementary
occupied set (complementarity fails for bonds the same way the
`square-bond-duality-tiny-torus` note shows bit-complement fails).  The
site-side dual_fail = 0 (reproduced from #606's ledger margins) versus
bond-side dual_fail ≠ 0 is itself the W2 laboratory statement: self-dual
location does not buy the duality constraint that site models get from
H1, so "p_c = 1/2 is free" is true while "the shape is duality-constrained"
is not automatic.  The bond Z is computed on the *raw* (r_b, r_w)
census (M = P20 - P02 over the eight rank pairs); it is exact as a census
and is NOT asserted to satisfy M(0) = -1, M(1) = +1 — indeed its M(1/2) is
-309/65536 ≈ -0.0047, and F(1) = 256/337 after renormalisation is reported
in the JSON (the F(1) ≠ 1 caveat is recorded rather than hidden).

## toy-families.json

`scripts/probe_invariant_shape/toy_location_theorem.py`

Machine check of the W1 theorem family: CDFs F_N with Q_N → p_*
uniformly on compacts (max deviation 0.008 at the largest N shown and → 0)
whose Z_N converges to a prescribed non-linear shape (family B: skew
deformation, Z limit ≠ linear in u, e.g. Z(0.9) → 1.4606 while the
linear-shape family A gives 1.45), and a parity-alternating family C whose
Z oscillates forever (Z(0.5) odd part ±0.0145 alternating) — no limit.
Verdict: a location theorem, even in its strongest uniform form, contains
zero information about shape.  W1 is a theorem at the level of "location
theorems in general".

## blindness-and-glift.json

`scripts/probe_invariant_shape/blindness_and_glift.py`

Direction 5: exact same-M / different-F counterexample on L=3
(P11 → (3/2)·P11 with P20, P02 fixed: M unchanged to the last Fraction, F
and Z both move).  Consequence: Z is a functional of the odd scalar
polynomial M alone; the even sector (the (1,1) rank-class mass) is
invisible to Z and visible to the pair (r_b, r_w).  An exact
two-measures-same-odd-content-different-Z pair — the probe's requested
counterexample object.

Direction 6: consensus g (cited from results/type582-residual/latest.json)
versus deltaZ = Z_{L=4} - Z_{L=3} on the same deciles: angle 93.9° raw,
159.6° after affine removal (i.e. after subtracting each vector's own
best-fit affine part).  Not a lift — the objects separate on the exact
laboratory.  Caveat (recorded, not hidden): L=3,4 are toy sizes; this
separates the pipeline's g from the exact-lab shape tangent, it does not
measure the angle in the N-regime.  The honest sentence: on the only exact
laboratory available, the naive lift fails; no affirmative evidence for a
lift exists anywhere in this probe.

## n725-zflow.json

`scripts/probe_invariant_shape/n725_zflow.py`

Per-batch anchored shape Z at N=725 from PR #614's committed histogram
(read-only, via git show), 100 batches, both orientation weightings:

```text
Z_mean (spin0): [-0.4679, -0.1270, 0.1131, (0), 0.5033, (1), 0.8895, 1.1233, 1.4503]
Z_se   (spin0):  7.7e-05  1.7e-05  1.2e-05       3.0e-05        1.2e-05  1.7e-05  7.9e-05
```

Readings:

* The shape at N=725 is measured to ~1e-4 or better per grid point.  It is
  NOT the L=3/L=4 shape (e.g. Z(0.9) = 1.4503 vs 1.3575/1.3833): two exact
  small sizes kill "S is already constant", and N=725 continues the drift.
* The three interior spacings have a strongly anisotropic principal
  direction (0.98 weight on the (0.7, 0.9) spacing) and adjacent spacings
  are anti-correlated (-0.87).  A single scalar common mode — one shape
  parameter breathing — is not what the batch structure shows; the
  dominant batch-level mode is a tail-vs-middle shear.  This is the W3
  evidence *against* a one-parameter S_N at finite N: the shape's own
  batch covariance is not rank-1.

Caveat: these are batch-covariance statements about one committed block;
they constrain measurability of S_N, they do not decide convergence of
E[Z_N].
