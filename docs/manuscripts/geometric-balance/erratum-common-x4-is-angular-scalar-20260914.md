# Erratum: the observed common `x≈4` correction is angular-scalar, not established H4/KdV

Date: 2026-09-14

Status: explicit correction to `kdv-identity-family-cancellation-target-20260914.md`.  That earlier note over-assigned the angular character of the common `x≈4` correction and should not be cited as the current mechanism.  The radial-dimension observation remains; the H4/KdV identification is withdrawn.

## 1. What was wrong

The earlier note took the large common `x≈4` correction in the average magnetic/safe gap and treated it as an H4/spin-four identity-family KdV correction.

The already committed oblique physical-gap data do not support that assignment.

For a direction `u=(a,b)`, use

```text
ell=n|u|,
E_phys=|u| I0,
2 pi x_m = 0.6544984694978736...
```

and form the dimension-four radial diagnostic

```text
C_common(ell,theta)
 = [E_phys - (2 pi x_m)/ell] ell^3.                    (1.1)
```

Using existing committed values gives approximately

```text
diagonal (1,1), n=4 : 0.2900
diagonal (1,1), n=5 : 0.3010
(2,1), n=3          : 0.3089
(2,1), n=4          : 0.3269
(3,1), n=3          : 0.3405
(3,2), n=2          : 0.3065
```

These values stay **positive and near one common radial amplitude** across orientations with both positive and negative `cos(4theta)`.

Dividing by `cos(4theta)` destroys the collapse and changes signs.  Thus the visible common correction is much more naturally angular H0/scalar than H4.

## 2. Corrected leading interpretation

The current hierarchy should be written as

```text
common average gap:
  angular H0,
  x≈4,
  candidate identity-family scalar such as T Tbar / scalar descendant mixture;

sector difference / root numerator:
  angular H4,
  x≈21/4,
  first visible noncommon block.
```

The ordinary common `x=4` correction and the leading noncommon H4 correction are therefore different angular sectors.

## 3. What survives from the “same representation” idea

A weaker conditional cancellation route remains valid in spirit.

If rank0 and rank2 critical sectors are two copies of the same Virasoro highest-weight representation, then a **scalar identity-family perturbation** whose diagonal matrix element is fixed by that representation has the same first-order shift in both copies, provided the microscopic coupling is common and no map/multiplicity operator acts differently.

For a `T\bar T`-type scalar perturbation the diagonal primary matrix element factorizes through the cylinder stress-tensor zero modes and depends only on `(h,\bar h,c)` under the standard assumptions.  Equal sector weights then imply equal scalar shift.

This is the corrected theorem template:

```text
same continuum primary copies
+ same scalar identity-family coupling
+ no nontrivial map-space matrix
=> common x=4 H0 shift cancels from the sector difference.
```

It is **not** a KdV/H4 theorem.

## 4. The actual null to prove is `D_0^(4)=0`, not `D_4^(4)=0`

In the difference-tangent notation, the visible common block is now typed as

```text
angular index a=0,
radial dimension x≈4.
```

Therefore the structural cancellation needed by the root is

```text
boxed:
D_0^(4)=0.                                             (4.1)
```

not the previously written `D_4^(4)=0`.

The first visible nonzero angular H4 difference coefficient remains

```text
D_4^(21/4) != 0
```

as the empirical/scaling target.

This distinction is important because it separates ordinary isotropic finite-size corrections from the anisotropic Matching-One root mechanism.

## 5. Consequence for the radial-dressing story

An H0 `x=4` scalar can still dress the H4 `x=21/4` coupling multiplicatively and generate an `ell^-2` radial correction **inside the H4 coefficient**:

```text
P4(ell)
 = a4 ell^-4 [1+c2 ell^-2+...].
```

That earlier dressing idea remains possible.

But, as already corrected elsewhere, an exact H4 angular projector removes the entire `P4(ell)H4(theta)` term including such dressing.  It cannot explain a post-H4 angular-orthogonal residual.

Thus two separate statements are now retained:

```text
common H0 x=4 field may radially dress P4;
post-H4 residual requires H0/H8/... or nonlinear H4^2 response.
```

## 6. Relation to Ward theory

The present high-value Ward question is no longer “is the common x4 H4 block a KdV charge?”

It is instead:

1. identify the actual scalar `x≈4` correction in the common magnetic/topological gap;
2. show why its matrix element is equal in rank0/rank2 extreme sectors;
3. separately derive the H4 `x=21/4` difference response.

A Virasoro identity-family argument may still solve item 2, but the relevant scalar quasiprimary/operator must be typed correctly before applying a charge formula.

## 7. Supersession statement

`kdv-identity-family-cancellation-target-20260914.md` is superseded as a Matching-One mechanism note because it assigned H4 character to the wrong observed block.

The standard KdV eigenvalue formula quoted there is correct CFT algebra; its application to the common `x≈4` lattice correction was unsupported by the angular data.

Do not use that note as evidence for `D_4^(4)=0`.

## 8. Claim boundary

Data-level correction:

- oblique common-gap `ell^-3` residual is approximately orientation-independent rather than proportional to `cos4theta` over existing controls.

Working interpretation:

- common block: H0, `x≈4`;
- leading noncommon block: H4, `x≈21/4`.

Open:

- exact continuum identity of the common H0 block;
- proof of equal rank0/rank2 matrix element;
- possible map-resolved scalar multiplicity.

The purpose of this erratum is to keep angular typing ahead of field naming, as required by the reverse audit.