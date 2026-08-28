# Hexagonal modular spin filter and Pell superconvergence

Status: C0/C1 theoretical prediction, frozen before any dedicated square-site Pell-hexagonal production.

## 1. Symmetry filter

At the square torus `tau=i`, the 90-degree automorphism allows spin multiples of 4. At

\[
\rho=e^{i\pi/3}=\frac12+i\frac{\sqrt3}{2},
\]

the torus has a 60-degree automorphism. A torus one-point amplitude of spin `s` can survive only if

\[
e^{is\pi/3}=1.
\]

Therefore spin 4 is forbidden at the hexagonal point while spin 12 is allowed. This is an operator discriminator independent of the norm-5 Gaussian multiplier.

## 2. Potts thermal one-point and the E4 conjecture

Roux, Ribault & Jacobsen, arXiv:2604.24491, show that the torus one-point function of the Potts energy operator `V^d_<1,2>` is a single conformal block,

\[
\langle V^d_{\langle1,2\rangle}\rangle=|\mathcal F|^2,
\qquad
\mathcal F=\eta(\tau)^{2\Delta_{(1,2)}}.
\]

At `Q=1`, `\Delta_(1,2)=5/8`, the percolation thermal weight. Standard torus Ward/Zhu recursion writes descendant one-point functions as modular/quasimodular differential operators acting on the primary one-point. `L_-4` contributes the Weierstrass invariant `g2`, i.e. weight-4/E4 structure.

For the repository's non-null thermal level-4 spin-4 quasiprimary a sharp conjecture is

\[
\frac{\langle Q_4 V_t\rangle}{\langle V_t\rangle}=C E_4(\tau).
\]

The exact proportionality is a theory target; the zero at `rho` follows already from the 60-degree automorphism.

## 3. Integer Pell sequence

Use positive Pell solutions

\[
x^2-3m^2=1
\]

and period vectors

\[
v_1=(2m,0),\qquad v_2=(m,x).
\]

Then

\[
\tau_m=\frac12+i\frac{x}{2m}\to\rho
\]

with error `O(m^-2)`. The first useful quotients are

| x | m | `N=|det P|=2mx` | Im tau |
|---:|---:|---:|---:|
| 7 | 4 | 56 | 0.875 |
| 26 | 15 | 780 | 0.8666666667 |
| 97 | 56 | 10864 | 0.8660714286 |
| 362 | 209 | 151316 | 0.8660287081 |

The committed generator gives

\[
E_4(\tau_m)/E_4(i)=0.0362565,\ 0.00265354,\ 0.000190777,\ldots
\]

and `m^2 E4(tau_m)/E4(i)` tends to about `0.59837`, as expected from a simple zero combined with Pell `m^-2` error.

## 4. Frozen H4 versus H12 consequence

Generic matching-odd H4 behavior is

\[
M_{odd}\sim L^{-13/4}.
\]

If its shape amplitude has the spin-4 zero at `rho`, the Pell family adds two powers:

\[
M_{odd}^{H4}(\tau_m)\sim L^{-21/4},
\qquad
\boxed{\delta p^{H4}_{Pell}\sim L^{-6}}.
\]

A same-radial-exponent spin-12 alias is allowed by the 60-degree symmetry and remains

\[
M_{odd}^{H12}\sim L^{-13/4},
\qquad
\delta p^{H12}\sim L^{-4}.
\]

Thus the Pell-hexagonal torus makes H4 versus H12 an exponent-level challenge rather than a small amplitude comparison.

## 5. Link to the historical L^-7 scalar sector

If H4 is correct, the finite Pell mismatch leaves an `L^-6` root term. A rotational-scalar matching-odd sector at `L^-7` would be next. An E4-weighted interpolation of shapes on opposite sides of the hexagonal point could cancel the linear shape mismatch and expose this scalar term. That interpolation must be frozen before its target is observed.

## 6. Execution gate

Norm-5 #57 remains the first compute priority. In parallel:

1. derive the exact modular operator for the level-4 quasiprimary;
2. validate the integer period matrices on exact controls;
3. after #57, use `N=56` as a small control and `N=780` as the first serious modular-spin target;
4. score fixed H4 `L^-6` and H12 `L^-4` root laws before any free exponent;
5. proceed to larger Pell members only if measured information-per-cost justifies it.

## References

- P. Roux, S. Ribault, J. L. Jacobsen, *Torus one-point functions in critical loop models*, arXiv:2604.24491 (2026).
- E. M. Brehm, I. Runkel, *Lattice models from CFT on surfaces with holes I*, J. Phys. A 55 (2022) 235001, torus descendant recursion.
