# Hexagonal modular spin filter and Pell superconvergence

Status: C0/C1 theoretical prediction, frozen before any dedicated square-site Pell-hexagonal production.

## 1. Why this is stronger than a generic shape scan

The current square-torus program lives at `tau=i`, where a 90-degree torus automorphism allows spin multiples of 4. A second elliptic point of the modular group is

\[
\rho=e^{i\pi/3}=\frac12+i\frac{\sqrt3}{2},
\]

whose torus has a 60-degree automorphism.

A torus one-point amplitude of spin `s` can be nonzero at `rho` only if

\[
e^{is\pi/3}=1.
\]

Therefore spin 4 is symmetry-forbidden at the hexagonal point, while spin 12 is allowed. This gives an operator discriminator orthogonal to the norm-5 Gaussian multiplier.

## 2. Connection to the 2026 Potts torus one-point result

Roux, Ribault & Jacobsen, arXiv:2604.24491, show that the torus one-point function of the Potts energy operator `V^d_<1,2>` is a single conformal block,

\[
\langle V^d_{\langle1,2\rangle}\rangle=|\mathcal F|^2,
\qquad
\mathcal F=\eta(\tau)^{2\Delta_{(1,2)}}.
\]

For `Q=1`, the loop weight is `n=1`; their parametrization gives `beta^2=2/3`, `c=0`, and `Delta_(1,2)=5/8`, the percolation thermal weight.

Standard torus Ward/Zhu recursion writes descendant one-point functions as modular/quasimodular differential operators acting on the primary one-point. In particular, `L_-4` contributes the Weierstrass invariant `g2`, i.e. weight-4/E4 structure. Since the repository's non-null thermal level-4 quasiprimary has spin 4, a sharp conjecture is

\[
\frac{\langle Q_4 V_t\rangle}{\langle V_t\rangle}
=C\,E_4(\tau)
\]

for an appropriate modular-covariant descendant normalization. The exact proportionality still needs derivation; the zero at `rho` follows already from the 60-degree automorphism.

## 3. Integer square-lattice Pell sequence

An integer square-lattice period matrix cannot realize `rho` exactly. Use positive Pell solutions

\[
x^2-3m^2=1
\]

and periods

\[
v_1=(2m,0),\qquad v_2=(m,x).
\]

Then

\[
\tau_m=\frac12+i\frac{x}{2m}\to\rho,
\]

with error `O(m^-2)`. The first useful quotients are

| x | m | site count `|det P|=2mx` | Im tau |
|---:|---:|---:|---:|
| 7 | 4 | 56 | 0.875 |
| 26 | 15 | 780 | 0.8666666667 |
| 97 | 56 | 10864 | 0.8660714286 |
| 362 | 209 | 151316 | 0.8660287081 |

The committed script computes the modular diagnostic. Numerically,

\[
E_4(\tau_m)/E_4(i)
=0.0362565,\ 0.00265354,\ 0.000190777,\ldots
\]

and

\[
m^2E_4(\tau_m)/E_4(i)\to0.59837\ldots
\]

as expected from the simple zero and Pell `m^-2` error.

## 4. Frozen H4 versus H12 prediction

At generic shape the empirical matching-odd H4 law is

\[
M_{\rm odd}\sim L^{-13/4}.
\]

If its shape amplitude has the expected spin-4 zero at `rho`, then along the Pell sequence

\[
M_{\rm odd}^{H4}(\tau_m)\sim L^{-13/4}L^{-2}=L^{-21/4}.
\]

After division by `M'~L^{3/4}`, the root bias becomes

\[
\boxed{\delta p^{H4}_{\rm Pell}\sim L^{-6}}.
\]

A same-radial-exponent spin-12 alias is allowed by the 60-degree torus symmetry and therefore remains

\[
M_{\rm odd}^{H12}\sim L^{-13/4},
\qquad
\delta p^{H12}\sim L^{-4}.
\]

Thus Pell-hexagonal tori turn H4 versus H12 into an exponent-level challenge, not a small amplitude comparison.

## 5. Relation to the historical L^-7 sector

If H4 is correct, the imperfect Pell shape leaves an `L^-6` root contribution. The conditional rotational-scalar matching-odd sector discussed in #47/#74 would enter at `L^-7`. Therefore:

1. first Pell sequence: expect `L^-6` if H4 survives and no larger scalar term intervenes;
2. an E4-weighted interpolation of shapes on opposite sides of `rho` could cancel the linear shape error and expose the scalar `L^-7` sector;
3. this would be a modular-shape analogue of Symanzik improvement.

The interpolation is a later target; do not fit it after looking at Pell outcomes.

## 6. Execution gate

This route should not delay norm-5. Recommended order:

1. derive the exact torus descendant operator and check whether the quasiprimary ratio is indeed `const*E4`;
2. validate the period-matrix geometry on exact controls;
3. after #57, run `N=56` as an engineering/control point and `N=780` as the first serious spin-filter point;
4. compare H4 `L^-6` and same-exponent H12 `L^-4` using frozen rules, not a free exponent first;
5. proceed to larger Pell members only if the information-per-cost is competitive.

## References

- P. Roux, S. Ribault, J. L. Jacobsen, *Torus one-point functions in critical loop models*, arXiv:2604.24491 (2026), especially the Potts energy one-point formula (5.11).
- E. M. Brehm, I. Runkel, *Lattice models from CFT on surfaces with holes I*, J. Phys. A 55 (2022) 235001, torus descendant recursion and `L_-4` one-point formula.
