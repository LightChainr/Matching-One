# Same-model diagonal sign flip of the `w^-4` charge-root correction

2026-09-14.  Deterministic safe-transfer control on the **same square-site NN / complementary matching model**, with the short cylinder period changed from an axis vector to a diagonal vector.

This is currently the sharpest direct test of the sector-odd spin-four anisotropy conjecture.  It avoids comparing different lattices and does not infer the spin from the exponent alone.

## 1. Spin-four prediction

Suppose the leading sector-odd lattice correction at criticality transforms as spin four.  If the short physical period makes angle `theta` with the square-lattice x-axis, the leading amplitude must carry the square harmonic

\[
\cos(4\theta).                                                \tag{1.1}
\]

The axial semi-infinite charge root has

\[
p_w^{axis}-p_c
\sim-\frac{A}{w^4},
\qquad A>0.                                                   \tag{1.2}
\]

For the diagonal direction

\[
\theta=\pi/4,
\qquad
\cos(4\theta)=-1.                                            \tag{1.3}
\]

Therefore the spin-four conjecture predicts **both**

\[
\boxed{p_n^{diag}-p_c>0}                                     \tag{1.4}
\]

and, after using the physical circumference

\[
\ell=n\sqrt2,
\]

\[
\boxed{
(p_n^{diag}-p_c)\ell^4
\longrightarrow A,
}                                                             \tag{1.5}
\]

the same positive amplitude that appears as `(pc-p_axis) w^4` on the axis.

A scalar sector-odd correction of the same scaling dimension would not predict this orientation sign flip.

## 2. Exact diagonal coordinate transfer

Use the integer basis

\[
u=(1,1),\qquad v=(0,1),\qquad \det(u,v)=1,                 \tag{2.1}
\]

so

\[
(x,y)=s u+t v=(s,s+t).                                       \tag{2.2}
\]

Quotient by

\[
s\sim s+n,                                                   \tag{2.3}
\]

which is exactly the physical period `n(1,1)`.  No rotation of the microscopic interaction is made.

In `(s,t)` coordinates the NN edges with positive `t` displacement are

\[
(\Delta s,\Delta t)=(0,1),\ (-1,1).                         \tag{2.4}
\]

The matching diagonals add

\[
(1,0),\qquad(-1,2).                                          \tag{2.5}
\]

Hence two frontier rows are sufficient for the full matching graph.  The transfer keeps integer lifted-`s` gains and rejects a transition immediately when a component acquires a cycle with nonzero deck gain.

The implementation is

`scripts/diagonal_charge_transfer.py`.

This is the same safe/void semantics as the axial transfer, only sliced in an integer oblique basis.

## 3. Charge roots

Using the safe Perron equality

\[
\lambda^0_{4,n}(p)
=\lambda^0_{8,n}(1-p),                                       \tag{3.1}
\]

the diagonal roots are

| `n` | period | `p_n^diag` | `p_n^diag-pc_ref` |
|---:|---|---:|---:|
| 2 | `(2,2)` | `0.5997254073143432` | `+6.9793565e-3` |
| 3 | `(3,3)` | `0.5937572147651211` | `+1.0111640e-3` |
| 4 | `(4,4)` | `0.5930452429807460` | `+2.9919219e-4` |
| 5 | `(5,5)` | `0.5928657129692908` | `+1.1966218e-4` |

The diagnostic reference is

\[
p_c^{ref}=0.59274605079.                                    \tag{3.2}
\]

It is not used to locate the roots.

The sign is already the opposite of the axial sequence.

## 4. Physical-circumference amplitude

Since

\[
\ell=n\sqrt2,
\qquad
\ell^4=4n^4,                                                  \tag{4.1}
\]

the scaled diagonal shifts are

| `n` | `(p_n^diag-pc_ref) ell^4` |
|---:|---:|
| 2 | `0.4466788` |
| 3 | `0.3276171` |
| 4 | `0.3063728` |
| 5 | `0.2991554` |

The independent axial transfer gives

\[
(p_c^{ref}-p_8^{axis})8^4=0.3001967,                          \tag{4.2}
\]

\[
(p_c^{ref}-p_9^{axis})9^4=0.2980447.                          \tag{4.3}
\]

Thus already at the available widths,

\[
\boxed{
(p_n^{diag}-p_c)\ell^4
\approx
(p_c-p_w^{axis})w^4
\approx0.30.}                                                 \tag{4.4}
\]

The diagonal sequence approaches this value from above while the axial sequence approaches it in the opposite charge-root direction.

Machine-readable values are in

`results/geometric-consistency/diagonal-charge-transfer-n2-n5-20260914.json`.

## 5. Interpretation

The same-model result is exactly the leading angular behavior expected from

\[
\boxed{
p_{\ell,\theta}^{ch}-p_c
\sim-\frac{A\cos(4\theta)}{\ell^4}.}                         \tag{5.1}
\]

At `theta=0`, (5.1) is negative.  At `theta=pi/4`, it is positive with equal magnitude.

This supplies three pieces of evidence simultaneously:

1. exponent four;
2. orientation sign flip;
3. physical-length amplitude agreement.

The last two are not explained by merely assigning an arbitrary scalar operator dimension `x_t+4`.

## 6. Relation to the square/kagome evidence

Jacobsen's square/kagome comparison suggested that rotational symmetry may remove the exponent-four correction on a three-/six-fold lattice.  The present experiment is conceptually cleaner in one respect: no lattice universality comparison is needed.

The microscopic model is held fixed and only the physical homology direction is changed.  A genuine spin-four correction must rotate with the cylinder orientation; a scalar one cannot.

Therefore the diagonal sign flip materially strengthens the spin-four interpretation in

`sector-odd-spin4-anisotropy-20260914.md`.

## 7. New directional conjecture

The natural all-angle refinement is

\[
\boxed{
A_4(\theta)=A\cos(4\theta)}                                  \tag{7.1}
\]

for the leading square-lattice dual-odd correction after expressing circumference in the isotropic physical correlation-length metric.

For a primitive direction `(a,b)`,

\[
\cos(4\theta)
=\frac{a^4-6a^2b^2+b^4}{(a^2+b^2)^2}.                       \tag{7.2}
\]

This gives strong future controls without an angle scan:

- `(1,0)`: `+1`;
- `(1,1)`: `-1`;
- `(2,1)`: `-7/25`;
- `(5,2)`: `41/841`, close to a leading-amplitude zero.

A `(5,2)` or nearby rational direction should strongly suppress the `ell^-4` term and expose the next correction, but its transfer has a larger row memory and should only be attempted if needed.

## 8. Important geometric normalization

The transfer step in an oblique integer basis need not have unit physical distance perpendicular to the circumference.  That factor affects individual excitation energies, but it cancels from the **root location**, which is a ratio of the sector mismatch and its thermal derivative.

For comparing root-shift amplitudes across orientations, the essential geometric normalization is the physical circumference `ell=|n u|`, as used above.

A future direct comparison of `Theta(pc)` amplitudes themselves must additionally normalize the longitudinal row spacing.

## 9. Claim boundary

The diagonal roots are deterministic finite-state transfer outputs.  The coordinate edge list is an exact rewrite of the square NN/matching graph.  The observed sign flip and physical-length amplitude agreement are numerical finite-width facts.

The asymptotic formula (5.1), its interpretation as a spin-four thermal-family correction, and the all-angle law (7.1) remain scaling conjectures.  Wider diagonal widths or one additional rational direction would strengthen the asymptotic case, but the existing same-model sign/amplitude test is already substantially more discriminating than an exponent-only fit.
