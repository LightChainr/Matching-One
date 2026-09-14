# A conditional KdV/Ward route to the common `x=4` cancellation

Date: 2026-09-14

Status: conditional continuum mechanism / proof target.  It does not assume a matching OPE parity.  The Virasoro KdV eigenvalue formula is standard; the unresolved Matching-One inputs are the sector dictionary and the identification/coupling of the lattice `x≈4` correction.

## 1. The lattice fact to explain

The two critical safe/topological sector energies share a large leading irrelevant correction compatible with a per-length `L^-3` term, i.e. total scaling dimension

```text
x≈4.
```

Yet their difference, which moves the Matching-One root, has no visible root term of exponent

```text
x-x_t = 4-5/4 = 11/4.
```

The first visible difference instead occurs much later, at root exponent four.

In the difference-tangent language the theorem target is

```text
D_4^(4)=0.
```

## 2. Identity-family spin-four correction as a KdV charge

The chiral level-four identity-family quasiprimary / quantum KdV charge has a zero mode whose eigenvalue on a Virasoro primary `|h>` is, up to the conventional cylinder scale,

```text
q3(h,c)
 = h^2
   -(c+2)h/12
   + c(5c+22)/2880.                                (2.1)
```

The anti-chiral charge has the same formula in `hbar`.

Thus a real square-lattice H4 perturbation built from the chiral/anti-chiral identity-family quasiprimaries has first-order diagonal energy shift

```text
delta E_a
 = g4 [q3(h_a,c)+q3(hbar_a,c)] L^-3             (2.2)
```

for a primary state `a`, modulo normalization of `g4`.

The important point is structural: this shift is fixed entirely by the Virasoro highest weights and the common microscopic coupling.

## 3. Conditional cancellation theorem

Assume:

1. the rank-0 and rank-2 critical long-cylinder sectors flow to two copies of the **same** Virasoro highest-weight state

```text
(h_0,hbar_0)=(h_2,hbar_2);
```

2. the leading `x=4`, H4 lattice perturbation is the same identity-family KdV/quasiprimary coupling `g4` in the physical square-site action for both topological sectors;
3. no additional map-resolved operator of the same `(x=4,H4)` block contributes a different matrix in the topological multiplicity space.

Then (2.2) gives exactly

```text
boxed:
delta E_0^(x=4,H4)=delta E_2^(x=4,H4),
```

so the sector difference has no first-order `x=4` contribution:

```text
boxed:
D_4^(4)=0.                                             (3.1)
```

This explains the absence of an `L^-11/4` root shift without assigning a matching-even scalar sign to the field.

## 4. Why the “same physical model” formulation matters

The matching-safe transfer is a computational representation of the complementary/rank-2 sector.  Physically, rank0 and rank2 are sectors of the **same square-site Bernoulli model** related by the finite digital-Alexander/matching dictionary.

Therefore a local lattice anisotropy of the original action should be treated as one microscopic coupling whose matrix elements are evaluated in two topological states, rather than as unrelated couplings in two different theories.

This is the correct setup for the KdV cancellation argument.

A transfer implementation must nevertheless transport the physical perturbation consistently through the matching representation; otherwise an apparent amplitude difference can be a source-coordinate artifact.

## 5. Map/multiplicity is the real remaining adversary

Equal conformal weights alone do not prove assumption 3.

If the continuum state space contains several copies/map sectors with the same Virasoro weights, a generic local perturbation may act as a nontrivial matrix in that multiplicity space.  A **pure Virasoro KdV charge** acts identically on equal highest-weight copies, but the lattice `x≈4` correction could contain additional map-resolved operators with the same radial/angular quantum numbers.

Thus the useful falsification question is:

```text
is the measured common x≈4,H4 block exhausted by the identity-family KdV direction?
```

This interfaces directly with #585's warning that modular covariance/spin alone do not define a one-dimensional physical solution space.

## 6. A concrete Ward/rank-projection test

Let `F_r` be the actual finite/continuum rank-sector functional.  For a fixed quasiprimary convention, the current Ward programme has a schematic combination

```text
A_r
 = F_r(U4+Ubar4)
   - C_E4(tau) F_r(epsilon),
```

where the second term removes coordinate/thermal contamination.

The targeted theorem is simply

```text
A_rank0 - A_rank2 = 0                                      (6.1)
```

for the identity-family `x=4,H4` block, including seam/contact terms appropriate to the nonlocal rank projection.

Equation (6.1), not an abstract OPE parity assignment, is exactly what the root mechanism needs.

## 7. What a failure would mean

If an independently normalized calculation finds a nonzero `x=4,H4` difference block, then at least one assumption fails:

```text
rank0/rank2 are not the same continuum state,
lattice x≈4 contains another map-resolved operator,
physical coupling is transported incorrectly,
or seam/contact terms distinguish the sectors.
```

A genuine nonzero block should produce an `L^-11/4` root contribution unless another cancellation removes it.  Its absence in data would then demand a second mechanism rather than be silently ignored.

## 8. Relation to the leading `x=21/4` block

The KdV argument only explains why the ordinary identity-family `x=4` anisotropy is common.  It does **not** identify the first noncommon block.

The observed H4 root exponent four still says that the first visible difference block has effective total dimension `21/4`.  Whether that block is the thermal-family level-four quasiprimary, an eight-arm/map-resolved sector, or a degenerate mixture remains a separate question.

## 9. Literature anchor

The quantum KdV charge `I3` can be written in Virasoro modes as

```text
I3 ~ 2 sum_(k>0) L_-k L_k
     + L0^2
     -(c+2)L0/12
     + c(5c+22)/2880,
```

so on a primary the positive-mode terms vanish and (2.1) follows.  This standard formula is enough for the conditional argument; no large-c or thermal assumption is used.

## 10. Claim boundary

Standard/exact CFT algebra:

- KdV `I3` primary eigenvalue depends only on `(h,c)`.

Conditional Matching-One consequence:

- `D_4^(4)=0` if both extreme sectors are the same highest-weight representation and the lattice x=4 block is the common identity-family KdV perturbation.

Open:

- rigorous/transfer proof of the sector dictionary at the required map resolution;
- exclusion of additional x=4 H4 operators in the multiplicity space;
- correct seam/contact term in the actual rank projection.

This is a deliberately narrow route: prove the amplitude equality actually needed by the root, rather than a stronger unconstructed matching automorphism.