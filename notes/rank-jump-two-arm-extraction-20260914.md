# From rank-jump-two attachment spines to macroscopic black arms

Date: 2026-09-14. Follow-up proof draft to `rank-jump-two-attachment-spines-20260914.md` on draft PR #773. The black-arm extraction below follows directly from covering-space geometry. The final black/white alternating-arm statement is isolated behind a standard planar site/matching separation lemma which should be checked in the repository's exact convention before promotion.

## 1. Covering-space setup

Let `T=Z^2/Lambda` be an honest NN square torus with systole `ell`, and let `v` be vacant in a configuration whose black ambient homology rank is zero. Let `C` be one black component touching neighbours of `v`.

Since every closed walk in C has zero deck displacement, the inverse image of C in the universal cover is a disjoint union of connected copies, each mapped isomorphically to C. Fix a lift `v_tilde` of v. Every attachment germ `a` of C at v enters one of these lifted copies; equivalently it has a deck-offset label `tau(a) in Lambda`, defined up to a common translation of all labels for C.

For attachments a,b, a path in C between them, closed by the two incident edges at v after opening v, has deck displacement

    tau(a)-tau(b).                                        (1)

This is the attachment-offset lemma from the companion note.

## 2. Type T forces three disjoint black arms to the systolic scale

Assume one outside component C has three attachment germs `a1,a2,a3` whose offset differences generate rank two. Then the three offsets are pairwise distinct. Therefore, when viewed from the fixed local lift `v_tilde`, the three occupied neighbour germs lie in **three distinct lifted copies**

    C_1, C_2, C_3                                          (2)

of the same quotient component C. Distinct lift copies are vertex-disjoint.

Consider C_i. Because the quotient component touches v through all three attachment types, the lifted copy C_i contains, besides the local attachment adjacent to `v_tilde`, a lift of another attachment adjacent to

    v_tilde + tau(a_j)-tau(a_i)                            (3)

for any chosen `j != i`. The vector in (3) is a nonzero period, hence has Euclidean length at least `ell`.

Since C_i is connected, there is a black path inside C_i from the local neighbour of `v_tilde` to a neighbour of the translated v in (3). Stop that path at its first exit from the ball `B(v_tilde,ell/3)`. Doing this for i=1,2,3 gives three black paths from the microscopic neighbourhood of v to radius `ell/3`.

They are mutually vertex-disjoint because they lie in the three distinct lift copies (2). Thus:

> **Proposition T-black.** A minimal three-contact Type-T rank-jump-two configuration forces three disjoint black NN arms from v to distance `ell/3` in the universal cover.

The constant `1/3` is arbitrary; any fixed number below `1/2` works after harmless lattice-distance bookkeeping.

The same argument works for a one-component rank-two offset set with more than three contacts: choose any three whose differences span rank two.

## 3. Type R forces four disjoint black arms

For a minimal split-loop Type-R event, there are two different quotient components C,D, each with two attachment offsets whose difference is a nonzero period. The two local attachment germs of C enter two different lifted copies of C; likewise for D. Quotient components are different and their lift copies are disjoint.

Each of the four copies contains a path from its local germ to a neighbour of a nonzero translate of v. Truncate at radius `ell/3` as above.

> **Proposition R-black.** A minimal Type-R split rank-jump-two configuration forces four mutually disjoint black NN arms from v to distance `ell/3`.

Nonminimal variants contain one of these mechanisms after choosing a generating subset of attachment offsets.

## 4. Alternating white separators: the remaining convention check

Inside the injective disk `B(v_tilde,ell/3)`, the problem is planar square-site percolation. In the closed-v configuration the black arms constructed above belong to distinct lifted black components. Choose simple representatives and order them cyclically around a small inner boundary.

The standard site/matching separation statement should imply that between consecutive disjoint black crossings of the annulus there is a white matching crossing connecting the inner and outer boundaries. Therefore one expects:

    Type T  => 3 black + 3 white alternating arms,
    Type R  => 4 black + 4 white alternating arms,          (4)

to a fixed fraction of the systole.

Before calling (4) a theorem in this repository, check explicitly:

1. the inner boundary convention around the toggled site;
2. white NN+NNN diagonal replacement in each square face;
3. simultaneous choice of all separating white crossings without accidental identification through a corner;
4. the chosen torus metric versus the injectivity radius.

These are local digital-planarity details, not an unresolved large-scale probability estimate.

## 5. What this does and does not imply for exponents

If the separator lemma is completed, Type T is a subset of a macroscopic polychromatic six-arm event and Type R a subset of an eight-arm event. This proves a **necessary arm count** for each attachment-spine class.

It does **not** prove that the Type-T probability has the bare six-arm exponent. The event also carries prescribed deck/homology information and may sit in a smaller map/connectivity sector with a higher scaling dimension or zero matching-odd amplitude.

That distinction is exactly why #768 remains necessary. The new finite result kills only the naive argument

    rank jump two => four black arms => 8-arm first.

It replaces it with the sharper question

    why does the abundant Type-T / six-arm geometry fail to contribute to
    the matching-odd scalar one-point sector, if indeed it fails?

## 6. Computational interface

The exact L=4 control already finds 176 Type-T three-contact configurations among 289 single-site rank-jump-two configurations (60.07% of the jump-two probability at the diagnostic p used in #769). L=5 in #769 should add the explicit deck-offset vectors and verify directly that the three selected lift-copy paths reach a fixed fraction of the systole.

No finite count is promoted to an arm exponent in this note.
