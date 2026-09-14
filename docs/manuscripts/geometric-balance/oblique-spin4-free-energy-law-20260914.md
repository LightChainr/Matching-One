# Oblique-cylinder evidence for a spin-four sector-odd free-energy operator

2026-09-14.  Deterministic safe-transfer evidence that the square-lattice primal/matching critical sector difference transforms as a genuine spin-four correction, not merely that the resulting pseudo-critical root happens to have exponent four.

## 1. Physical normalization in an integer oblique basis

Let `u=(a,b)` be a primitive short-period direction and choose a Bezout complement `v` with `det(u,v)=1`.  Quotient by `n u` and transfer by one integer `v` step.

The physical circumference is

\[
\ell=n|u|.                                                     \tag{1.1}
\]

One transfer step advances by perpendicular physical height

\[
h_\perp=\frac1{|u|}.                                         \tag{1.2}
\]

Therefore a per-row transfer excitation `I_row` corresponds to physical energy

\[
E_{phys}=|u|I_{row}.                                          \tag{1.3}
\]

This normalization is independently checked by the magnetic primary:

\[
\boxed{n|u|^2 I^0\to2\pi\frac5{48}}                          \tag{1.4}
\]

for axis, diagonal and several oblique directions; see `oblique-magnetic-metric-controls-20260914.json`.

## 2. Spin-four critical mismatch prediction

Let

\[
\Theta_{u,n}(p)
=I^0_{4,u,n}(p)-I^0_{8,u,n}(1-p)                              \tag{2.1}
\]

be the per-transfer-step charge free-energy difference.

If the leading dual-odd critical correction is a spin-four thermal-family descendant of total dimension

\[
x_{odd}=x_t+4=21/4,                                          \tag{2.2}
\]

then in physical units

\[
E_{odd}^{phys}(p_c)
\sim B\cos(4\theta_u)\,\ell^{1-x_{odd}}
=B\cos(4\theta_u)\,\ell^{-17/4}.                             \tag{2.3}
\]

Using (1.3), this predicts the per-row quantity

\[
\boxed{
\Theta_{u,n}(p_c)
\sim
\frac{B}{|u|}\cos(4\theta_u)\,\ell^{-17/4}.}                 \tag{2.4}
\]

Thus the orientation-independent coupling estimate is

\[
\boxed{
B_{est}
=\frac{|u|\Theta_{u,n}(p_c)\ell^{17/4}}
       {\cos(4\theta_u)}.}                                    \tag{2.5}
\]

A scalar correction of dimension `21/4` would have no `cos(4theta)` factor and fails this test.

## 3. Deterministic transfer values

Using reference

\[
p_c^{ref}=0.59274605079,                                     \tag{3.1}
\]

which is not used in locating any charge root, the oblique safe transfers give:

| direction `u` | `n` | `Theta_row(pc_ref)` | `B_est` |
|---|---:|---:|---:|
| `(1,1)` | 4 | `-4.7096255e-4` | `1.05183` |
| `(1,1)` | 5 | `-1.7703612e-4` | `1.02068` |
| `(2,1)` | 3 | `-4.0366939e-5` | `1.05058` |
| `(2,1)` | 4 | `-1.1559828e-5` | `1.02175` |
| `(3,1)` | 3 | `+6.1566873e-6` | `0.98845` |
| `(3,2)` | 2 | `-4.5029268e-5` | `1.02168` |

The axis controls give approximately

\[
B_{axis}(w=8)=1.0216,
\qquad
B_{axis}(w=9)=1.0124.                                        \tag{3.2}
\]

The sign of `Theta` flips exactly with `cos(4theta)`:

- axis `(1,0)`: positive;
- diagonal `(1,1)`: negative;
- `(2,1)`: negative;
- `(3,1)`: positive;
- `(3,2)`: negative.

More importantly, division by the exact angular harmonic collapses the amplitudes to one number near `B~1.02` across unrelated row memories and physical circumferences.

Machine-readable values:

`results/geometric-consistency/oblique-spin4-free-energy-amplitude-20260914.json`.

## 4. The charge-root angular law follows as a ratio

The thermal derivative is a scalar to leading order after physical metric normalization:

\[
\partial_p E_{charge}^{phys}(p_c)
\sim C\ell^{-1/4}.                                            \tag{4.1}
\]

Therefore solving `Theta(p_root)=0` gives

\[
\boxed{
p_{root}-p_c
\sim
-\frac BC\frac{\cos(4\theta)}{\ell^4}.}                     \tag{4.2}
\]

The observed root amplitude is

\[
A=B/C\approx0.30,                                             \tag{4.3}
\]

consistent with the independent axial pivotal/thermal slope amplitude.

Thus the root-level `cos(4theta)` law is not the primitive observation.  It is the quotient of

1. a spin-four critical sector mismatch;
2. a leading scalar thermal response.

The free-energy data directly resolve item 1.

## 5. Strong exclusion of the scalar `x=21/4` alternative

Before the orientation test, the exponents alone admitted a numerical coincidence:

\[
x_{odd}=21/4                                                  \tag{5.1}
\]

could have been assigned to a scalar object, for example by over-reading a Kac-table dimension.

The oblique data disfavor that explanation in three independent ways.

### Sign

A scalar mismatch cannot reverse sign under a 45-degree change of the cylinder direction while the microscopic lattice and thermal parameter are unchanged.

### Magnitude

The `(2,1)` and `(3,2)` amplitudes follow the nontrivial rational values of `cos4theta`, not merely `+/-1` between axis and diagonal.

### Near-node suppression

The `(5,2)` direction has

\[
\cos4\theta=41/841\approx0.04875.                            \tag{5.2}
\]

Its charge root is already within about `1e-6` of `p_c_ref` at physical circumference only about `10.77`, consistent with suppression of the leading spin-four term.

These are operator-transformation signatures, not exponent fitting.

## 6. Connection to the thermal Kac module

The thermal primary is `phi_{2,1}`, with `h=5/8` and a level-two null relation.  There is no independent level-two thermal quasiprimary after quotienting by the null state.  The first relevant even-spin nonredundant chiral descendant occurs at level four.

This makes the spin-four interpretation representation-theoretically natural:

\[
\boxed{
\mathcal O_{4}^{odd}
\sim
Q_4\phi_t\otimes\bar\phi_t
+\phi_t\otimes\bar Q_4\bar\phi_t,}                           \tag{6.1}
\]

where `Q_4` denotes the level-four quasiprimary descendant, schematically.  The real square-lattice combination produces the `cos(4theta)` harmonic.

Equation (6.1) is a structural candidate, not a normalized operator identification; logarithmic mixing and precise quasiprimary normalization remain to be worked out.

## 7. Pell null experiment

The angular numerator factorizes:

\[
a^4-6a^2b^2+b^4
=(a^2-2ab-b^2)(a^2+2ab-b^2).                                 \tag{7.1}
\]

Pell directions satisfying

\[
a^2-2ab-b^2=\pm1                                             \tag{7.2}
\]

approach the exact spin-four node `theta=pi/8` with

\[
\cos4\theta=O(|u|^{-2}).                                     \tag{7.3}
\]

Along such a sequence the spin-four contribution is demoted from `ell^-4` to `ell^-6` in the root shift.  Alternating Pell signs reverse that residual contribution, providing a clean future separation of spin-four and true spin-six terms inside one microscopic model.

## 8. Claim boundary

The oblique coordinate transforms, safe Perron values, physical normal-height factor and angular trigonometric factors are deterministic.  The collapse of `B_est` is finite-width numerical evidence.

The asymptotic operator statement (2.3), the normalized level-four thermal quasiprimary identification, and the Pell separation of the next spin-six coefficient remain scaling/LCFT conjectures.  They are now substantially more constrained than an exponent-only hypothesis.
