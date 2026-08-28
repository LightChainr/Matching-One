# Exact E4 torus one-point of the c=0 thermal level-4 quasiprimary

Status: exact Virasoro/torus one-point identity for the `c=0, h=5/8` degenerate module, with the convention for `g2` stated below. This does **not** by itself prove that the lattice matching residual couples only to this field.

## 1. Input states

At `c=0, h=5/8`, the level-2 null state is

\[
\chi_2=\left(L_{-2}-\frac23L_{-1}^2\right)|h\rangle=0.
\]

Modulo its null submodule, the non-null level-4 quasiprimary used in the repository is

\[
Q_4|h\rangle=
\left(40L_{-2}^2-60L_{-3}L_{-1}-9L_{-4}\right)|h\rangle.
\]

Its formal Shapovalov norm is nonzero (`4930` in the repository normalization).

## 2. Translation identities on a torus one-point function

For any descendant state `psi`, translation invariance of a one-point function gives

\[
\langle L_{-1}\psi\rangle_{\tau}=0.
\]

Using `[L_m,L_n]=(m-n)L_{m+n}`,

\[
L_{-3}L_{-1}=L_{-1}L_{-3}-2L_{-4},
\]

hence

\[
\langle L_{-3}L_{-1}\phi\rangle=-2\langle L_{-4}\phi\rangle.
\]

Likewise `[L_{-2},L_{-1}]=-L_{-3}` gives

\[
L_{-2}L_{-1}^2
=L_{-1}^2L_{-2}-L_{-1}L_{-3}-L_{-3}L_{-1},
\]

so all terms beginning with `L_-1` vanish in the torus one-point function and

\[
\langle L_{-2}L_{-1}^2\phi\rangle
=2\langle L_{-4}\phi\rangle.
\]

Now apply `L_-2` to the null relation. Since `L_-2 chi_2=0`,

\[
\langle L_{-2}^2\phi\rangle
=\frac23\langle L_{-2}L_{-1}^2\phi\rangle
=\frac43\langle L_{-4}\phi\rangle.
\]

## 3. Collapse of Q4 to L_-4 in a torus one-point function

Therefore

\[
\begin{aligned}
\langle Q_4\phi\rangle
&=40\frac43\langle L_{-4}\phi\rangle
-60(-2)\langle L_{-4}\phi\rangle
-9\langle L_{-4}\phi\rangle\\
&=\frac{493}{3}\langle L_{-4}\phi\rangle.
\end{aligned}
\]

Brehm--Runkel's torus descendant recursion gives, for a primary of weight `h`,

\[
\langle L_{-4}\phi\rangle_{\tau}
=\frac{h g_2(\tau)}{20}\langle\phi\rangle_{\tau}.
\]

At `h=5/8`,

\[
\boxed{
\frac{\langle Q_4\phi\rangle_{\tau}}
{\langle\phi\rangle_{\tau}}
=\frac{493}{96}g_2(\tau).
}
\]

For periods `(1,tau)`, with the standard convention

\[
g_2(\tau)=\frac{4\pi^4}{3}E_4(\tau),
\]

this is

\[
\boxed{
\frac{\langle Q_4\phi\rangle_{\tau}}
{\langle\phi\rangle_{\tau}}
=\frac{493\pi^4}{72}E_4(\tau).
}
\]

Thus the weight-4 modular factor is not merely dimensional reasoning: it follows exactly from the null relation plus torus translation/Ward recursion.

## 4. Independent consistency with the Potts energy block

Roux--Ribault--Jacobsen (arXiv:2604.24491, eq. 5.11) find that the Potts energy one-point conformal block is

\[
\mathcal F(\tau)=\eta(\tau)^{2h}.
\]

At percolation `h=5/8`. Its Ramanujan--Serre derivative vanishes,

\[
\mathcal D_h\eta^{2h}=0,
\]

which is the modular form of the same level-2 null relation. A modular-covariant level-4 descendant divided by the primary block must therefore be a weight-4 modular form, and the explicit calculation above fixes it to `E4` with nonzero coefficient.

## 5. Consequences

1. At the exact hexagonal modulus `rho=exp(i*pi/3)`, `E4(rho)=0`, so this quasiprimary's torus one-point vanishes exactly.
2. On positive Pell approximants to `rho`, the simple `E4` zero gives an additional `O(L^-2)` suppression.
3. If the lattice matching-odd residual is dominated by this field, its generic `L^-13/4` residual becomes `L^-21/4` on the Pell sequence, and the associated root bias becomes `L^-6` after division by the thermal slope.
4. A spin-12 same-radial alias is not killed by the 60-degree torus automorphism and does not inherit this `E4` zero.

## 6. Evidence boundary

The identity proved here concerns the continuum Virasoro descendant one-point function. The remaining nontrivial bridge is **lattice coupling**: showing that the matching-odd lattice observable has a nonzero overlap with this quasiprimary and that lower/competing nonlocal sectors do not dominate. Norm-5 and Pell-modulus spectroscopy are direct tests of that bridge.

## References

- P. Roux, S. Ribault, J. L. Jacobsen, *Torus one-point functions in critical loop models*, arXiv:2604.24491 (2026), eq. (5.11).
- E. M. Brehm, I. Runkel, *Lattice models from CFT on surfaces with holes I*, J. Phys. A 55 (2022) 235001, eq. (4.52) and Appendix A.
