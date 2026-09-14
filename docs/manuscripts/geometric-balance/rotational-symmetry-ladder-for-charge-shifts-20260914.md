# A rotational/Kac-module selection rule for semi-infinite charge-root shifts

2026-09-14.  Corrected conjectural synthesis for the sector-odd charge correction.  This note replaces the earlier naive rule `s_*=lcm(2,k)`: point-group symmetry alone is not enough.  One must also remove descendants that are null or redundant in the percolation thermal conformal family.

The correction was forced by the level-two null vector of the thermal field `phi_{2,1}` and by the new same-model oblique-cylinder spin-four controls.

## 1. Charge-root shift as a sector-odd thermal-family anisotropy

Assume the leading primal/matching sector difference at criticality comes from a nonredundant chiral descendant of the thermal primary.  Percolation has

\[
h_t=\bar h_t=5/8,
\qquad x_t=5/4.                                               \tag{1.1}
\]

A chiral descendant at level `s` and its reflected anti-chiral partner have total dimension

\[
x_{odd}=x_t+s                                                \tag{1.2}
\]

and spins `+s,-s`.  A real lattice perturbation uses the reflection-even combination and hence produces an angular harmonic such as `cos(s theta)`.

The critical sector mismatch scales as

\[
\Theta_w(p_c)\asymp w^{1-(x_t+s)},                           \tag{1.3}
\]

while the thermal slope scales as

\[
\Theta'_w(p_c)\asymp w^{1-x_t}.                              \tag{1.4}
\]

Therefore a surviving level-`s` anisotropy gives

\[
\boxed{p_w-p_c\asymp w^{-s}.}                               \tag{1.5}
\]

The exponent is a descendant **level/spin only after null and redundant states have been removed**.

## 2. The thermal level-two state is not an independent spin-two field

The thermal primary is the degenerate Kac field

\[
\phi_t=\phi_{2,1}.                                            \tag{2.1}
\]

Its level-two singular-vector relation is

\[
\left(
L_{-2}-\frac{3}{2(2h_t+1)}L_{-1}^2
\right)|t\rangle=0.                                         \tag{2.2}
\]

For `h_t=5/8`,

\[
\boxed{L_{-2}|t\rangle=\frac23L_{-1}^2|t\rangle.}           \tag{2.3}
\]

Thus after quotienting by the null state there is no independent level-two quasiprimary.  The would-be spin-two thermal descendant is only a derivative/coordinate redundancy.

This removes the previous proposed `C2 -> Delta=2` mechanism.

It also aligns with the geometric warning already present in the first version of this note: ordinary spin-two anisotropy mostly changes the continuum metric.  The Virasoro null relation shows more sharply why there is no independent thermal-family spin-two correction to promote after metric calibration.

Primary null-vector references are standard Virasoro/Kac theory; the generic `phi_{2,1}` level-two relation is also used in BPZ differential-equation derivations.

## 3. First nonredundant low levels of the `phi_{2,1}` quotient

At a purely counting level, the Verma dimensions at levels `n=0,1,2,3,4` are

\[
1,1,2,3,5.                                                    \tag{3.1}
\]

Removing the level-two null module subtracts the partition numbers at level `n-2`, leaving

\[
1,1,1,2,3.                                                    \tag{3.2}
\]

Modulo total `L_{-1}` derivatives, the first new quasiprimary content is therefore

- no independent level-two quasiprimary;
- one level-three quasiprimary;
- one new level-four quasiprimary.

This is the representation-theoretic reason the point-group rule must start from the **actual thermal Kac module**, not from all integer spins.

A full logarithmic-module treatment may refine higher levels; only the low-level null structure needed for the selection statements below is used here.

## 4. Corrected point-group ladder

Let `C_k` be the microscopic rotational subgroup.  The candidate level must both

1. be invariant under the point group, `s=0 mod k`;
2. exist as a nonnull/nonredundant thermal-family quasiprimary.

Reflection does **not** by itself exclude odd spin: it exchanges the `+s` and `-s` chiral descendants, and their real sum is reflection even.

For the low symmetries relevant here this gives

| microscopic rotational symmetry | first thermal-family candidate | predicted charge-root shift |
|---|---:|---:|
| no nontrivial rotation (`C1`, reflection allowed) | spin 3 | `w^-3` |
| `C2` / `D2` | spin 4 | `w^-4` |
| `C3` / `D3` without 60-degree rotation | spin 3 | `w^-3` |
| `C4` / `D4` | spin 4 | `w^-4` |
| `C6` / `D6` | spin 6 | `w^-6` |

Exact self-matching/self-duality can annihilate the complete sector-odd amplitude and override the table.

The old claim `C2 -> w^-2` is withdrawn.

## 5. Square site: exponent, sign and angular amplitude now all resolve spin four

The square lattice has `C4`.  The axial safe transfer gives

\[
p_w^{axis}-p_c\sim-\frac{A}{w^4},
\qquad A\approx0.30.                                         \tag{5.1}
\]

The new oblique safe transfers keep the microscopic model fixed and rotate only the homology direction.  If the correction has spin four, the leading law is

\[
\boxed{
p_{n,u}^{ch}-p_c
\sim
-\frac{A\cos(4\theta_u)}{(n|u|)^4}.}                         \tag{5.2}
\]

This prediction is supported quantitatively by several independent directions.

### Axis `(1,0)`

`cos(4theta)=+1`.  Widths 8--9 give

\[
A_{est}=0.30020,\ 0.29804.                                   \tag{5.3}
\]

### Diagonal `(1,1)`

`cos(4theta)=-1`; the root moves to the **opposite side** of `p_c`.  At `n=5`, after using physical circumference `ell=5sqrt2`,

\[
A_{est}=0.29916.                                              \tag{5.4}
\]

### Direction `(2,1)`

\[
\cos4\theta=-7/25=-0.28.                                    \tag{5.5}
\]

At `n=4`,

\[
A_{est}=0.30077.                                              \tag{5.6}
\]

### Direction `(3,2)`

\[
\cos4\theta=-0.7041420118\ldots                              \tag{5.7}
\]

and already at `n=2`,

\[
A_{est}=0.29956.                                              \tag{5.8}
\]

These are much more discriminating than an exponent-four fit: a scalar field of dimension `x_t+4` would not produce the observed orientation sign/magnitude law.

Data and scripts:

- `diagonal-spin4-charge-root-20260914.md`;
- `scripts/diagonal_charge_transfer.py`;
- `scripts/oblique_charge_transfer.py`;
- `results/geometric-consistency/oblique-spin4-angular-controls-20260914.json`.

## 6. Kagome/hexagonal symmetry remains a spin-six control

Jacobsen's kagome bond eigenvalue sequence has the square-like exponent-four amplitude absent and a leading shift compatible with exponent six.  The kagome/hexagonal bulk symmetry includes 60-degree rotation, which forbids spin four and allows spin six.

This remains consistent with the corrected Kac-module rule:

\[
C_6:\quad s_*=6.                                              \tag{6.1}
\]

The triangular-site self-matching model is stronger still: complement symmetry fixes the charge root to `1/2` at every width, so every sector-odd amplitude vanishes.

## 7. A new distinction: `C3` is not `C6`

The old table bundled three-fold symmetry plus reflection with six-fold symmetry.  That is not generally justified.

A spin-three pair transforms trivially under a `120 degree` rotation:

\[
e^{\pm i3(2\pi/3)}=1,                                       \tag{7.1}
\]

and the real combination is reflection even.  Since a nontrivial level-three thermal quasiprimary survives the level-two null quotient, a genuine `D3` lattice without 60-degree rotation can in principle have

\[
\boxed{p_w-p_c\asymp w^{-3}.}                               \tag{7.2}
\]

A `D6` lattice forbids it and first admits spin six.

This gives a sharper future lattice-symmetry test than the earlier `C2` proposal.

## 8. Pell directions create a same-model spin-four null experiment

The square harmonic factorizes arithmetically:

\[
a^4-6a^2b^2+b^4
=(a^2-2ab-b^2)(a^2+2ab-b^2).                                 \tag{8.1}
\]

The exact spin-four zero direction is

\[
\frac ab=1+\sqrt2,
\qquad \theta=\pi/8.                                         \tag{8.2}
\]

Choose primitive Pell approximants satisfying

\[
a^2-2ab-b^2=\pm1.                                           \tag{8.3}
\]

Then

\[
\cos4\theta=O(|u|^{-2}).                                     \tag{8.4}
\]

Consequently the nominal `ell^-4` spin-four root shift is geometrically suppressed to order `ell^-6` along this sequence.  The two Pell signs flip the residual spin-four contribution.

This suggests a powerful same-model experiment:

- average the `+1` and `-1` Pell subsequences to expose a genuine spin-six term;
- difference them to isolate the residual spin-four anisotropy.

The existing `(5,2)` control already lies near this null direction (`cos4theta=41/841`) and its charge root is within about `10^-6` of the reference `p_c` at only `n=2`, despite physical circumference about `10.77`.

A wider Pell computation is a targeted follow-up, not a generic angle scan.

## 9. Revised general rule

The correct conjectural selection principle is

\[
\boxed{
\Delta_{charge}=s_*,
\quad
s_*=\min\{s>0:\ s\text{ is point-group allowed and is a nonredundant thermal-family quasiprimary level}\}.} \tag{9.1}
\]

This is stronger and safer than the old arithmetic rule `lcm(2,k)`.

It separates three mechanisms that must not be conflated:

1. continuum metric/stress-tensor anisotropy;
2. Virasoro-null or total-derivative descendants;
3. genuine sector-odd thermal-family anisotropy.

Only the third sets the charge-root shift after the first two are removed.

## 10. Claim boundary

The `phi_{2,1}` level-two null vector is standard CFT structure.  The oblique charge roots are deterministic lattice calculations.  Point-group invariance of a spin harmonic is exact.

The identification of the leading charge correction with a thermal-family quasiprimary, the higher-level logarithmic-module content, the `D3 -> 3` and Pell spin-six separation predictions remain conjectural.  The previous `C2 -> 2` prediction is explicitly superseded by this corrected note.
