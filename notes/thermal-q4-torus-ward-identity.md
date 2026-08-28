# Exact torus Ward fingerprint of the thermal level-4 spin-4 candidate

**Status:** exact CFT descendant identity conditional on the ordinary `c=0, h=5/8` thermal module; application to the square-site matching residual remains a C0 lattice-to-CFT identification hypothesis.

This note upgrades one part of the current `x=21/4` story from a modular-form guess to an explicit Virasoro calculation. It also sharpens how the norm-5 H4/H12 experiment should be interpreted.

The central conclusion is:

\[
\boxed{
\frac{\langle Q_4\,\bar\phi\rangle_{\tau}}
     {\langle \phi\bar\phi\rangle_{\tau}}
=\frac{493}{96}\,g_2(\tau)
}
\]

for the repository normalization

\[
Q_4=40L_{-2}^2-60L_{-3}L_{-1}-9L_{-4}
\]

acting on an ordinary primary `phi` with `c=0`, `h=5/8` and null vector

\[
(L_{-2}-\tfrac23L_{-1}^2)\phi=0.
\]

Here `g2` is the Weierstrass invariant of the torus period lattice. For periods `(1,tau)`, `g2` is proportional to the weight-4 Eisenstein series `E4(tau)`. Therefore this particular spin-4 torus one-point fingerprint vanishes exactly at the equianharmonic/hexagonal elliptic point.

The identity is about the candidate continuum operator. It does **not** by itself prove that the matching residual couples only to this operator.

## 1. Inputs already present in this repository

The existing note `thermal-level4-spin4-candidate.md` establishes that the percolation thermal primary has

```text
c = 0
h = hbar = 5/8
x_t = 5/4
```

and has the level-2 null relation

\[
\chi_2=(L_{-2}-\tfrac23L_{-1}^2)|h\rangle=0.
\]

After quotienting the null submodule, the level-4 chiral module contains the non-null quasiprimary

\[
Q_4=(40L_{-2}^2-60L_{-3}L_{-1}-9L_{-4})|h\rangle,
\]

which gives the bulk weights `(37/8,5/8)` and hence

```text
x = 21/4
spin = +4
```

plus the antiholomorphic conjugate.

## 2. External ingredient: torus Ward recursion

Brehm and Runkel, arXiv:2112.01563, give a recursive torus Ward identity for Virasoro descendants. Two consequences are enough here.

For a primary `phi` of holomorphic weight `h`,

\[
\langle L_{-4}\phi\rangle
=\frac{g_2 h}{20}\langle\phi\rangle.
\]

Also, one-point translation invariance implies that an outer `L_-1` insertion vanishes:

\[
\langle L_{-1}\Psi\rangle=0
\]

for every descendant state `Psi` inserted at the single marked point.

These identities are local to the Virasoro module and do not require assuming that the level-4 answer spans a one-dimensional modular-form space.

## 3. Exact reduction of the repository Q4

Write

\[
C=\langle L_{-4}\phi\rangle.
\]

### 3.1 The `L_-3 L_-1` term

Using

\[
[L_{-1},L_{-3}]=2L_{-4}
\]

and translation invariance,

\[
0=\langle L_{-1}L_{-3}\phi\rangle
=2C+\langle L_{-3}L_{-1}\phi\rangle.
\]

Hence

\[
\boxed{\langle L_{-3}L_{-1}\phi\rangle=-2C.}
\]

### 3.2 The `L_-2^2` term

The null relation gives

\[
\langle L_{-2}^2\phi\rangle
=\frac23\langle L_{-2}L_{-1}^2\phi\rangle.
\]

Since

\[
L_{-2}L_{-1}=L_{-1}L_{-2}-L_{-3},
\]

translation invariance removes the first term after taking the one-point function, so

\[
\langle L_{-2}L_{-1}^2\phi\rangle
=-\langle L_{-3}L_{-1}\phi\rangle
=2C.
\]

Therefore

\[
\boxed{\langle L_{-2}^2\phi\rangle=\frac43 C.}
\]

### 3.3 Assemble Q4

Thus

\[
\begin{aligned}
\langle Q_4\phi\rangle
&=40\left(\frac43C\right)-60(-2C)-9C\\
&=\frac{493}{3}C.
\end{aligned}
\]

For `h=5/8`, the torus Ward result gives

\[
C=\frac{g_2}{32}\langle\phi\rangle.
\]

Hence

\[
\boxed{
\langle Q_4\phi\rangle
=\frac{493}{96}g_2(\tau)\langle\phi\rangle.
}
\]

Tensoring the untouched antiholomorphic thermal primary gives the bulk spin `+4` result; complex conjugation gives spin `-4`.

## 4. Why the 2026 Potts torus result matters

Roux, Ribault and Jacobsen, arXiv:2604.24491, show that the Potts energy field `V^d_<1,2>` has a torus one-point function built from a single conformal block. Their chiral block is

\[
F_{\rm energy}(\tau)=\eta(\tau)^{2\Delta_{1,2}}.
\]

At percolation `Delta_(1,2)=5/8`, so the ordinary energy block is proportional to `eta^(5/4)`.

This is useful for two reasons:

1. the primary one-point function that appears in the Ward identity is a concrete Potts torus object rather than an abstract nonzero function;
2. the level-4 descendant shape is not an arbitrary weight-4 function multiplying that block: for the explicit repository `Q4`, the Ward/null-state algebra fixes it to `g2` exactly.

So the previous minimal hypothesis

```text
F_Q4 / F_energy proportional to E4
```

is stronger now: within the ordinary thermal module it follows directly from Ward recursion.

## 5. Exact elliptic-point annihilation of this Q4 fingerprint

Let

\[
\omega=e^{i\pi/3}=\frac12+i\frac{\sqrt3}{2}.
\]

The equianharmonic torus has

\[
g_2(\omega)=0
\]

or equivalently `E4(omega)=0`.

Therefore

\[
\boxed{\langle Q_4\bar\phi\rangle_{\omega}=0.}
\]

This is exactly what the order-3 elliptic stabilizer suggests: a scalar torus background cannot support a first-order spin-4 one-point response at the hexagonal fixed point, while spin 12 is not removed by that stabilizer.

This makes the hexagonal point an operator filter orthogonal to the current norm-5 experiment:

```text
norm 5 at tau=i:     changes the microscopic angle, fixed continuum modulus
hexagonal filter:    changes continuum modulus so the spin-4 torus response vanishes
```

## 6. Integer square-lattice route: Pell approach to the elliptic zero

An exact hexagonal period basis cannot be generated by two integer square-lattice vectors, but it can be approached unusually cleanly.

Take periods

```text
u = (2b,0)
v = (b,a)
```

so that

\[
\tau_{a,b}=\frac12+i\frac{a}{2b},
\qquad
N=\det(u,v)=2ab.
\]

Choose a fixed-defect Pell family

\[
a^2-3b^2=D.
\]

Then

\[
\tau_{a,b}-\omega
=i\frac{a-\sqrt3b}{2b}
=i\frac{D}{2b(a+\sqrt3b)}
=O(N^{-1}),
\]

and more sharply

\[
\boxed{N(\tau_{a,b}-\omega)\to iD/2.}
\]

Using the Ramanujan derivative identity

\[
\frac{1}{2\pi i}E_4'
=\frac13(E_2E_4-E_6),
\]

at `omega`, where `E4(omega)=0`, gives

\[
E_4'(\omega)=-\frac{2\pi i}{3}E_6(\omega).
\]

Therefore

\[
\boxed{
N E_4(\tau_{a,b})
\longrightarrow
\frac{\pi D}{3}E_6(\omega).
}
\]

With the standard q-series normalization,

```text
E6(omega) = 2.881541100790945623...
(pi/3) E6(omega) = 3.017542784420626905...
```

so the two useful Pell families have parameter-free opposite asymptotic signs:

```text
D=+1:  N E4 -> +3.0175427844206269...
D=-2:  N E4 -> -6.0350855688412538...
```

A direct q-series evaluation gives:

| `(a,b)` | `D` | `N` | `E4(tau)` | `N E4(tau)` |
|---|---:|---:|---:|---:|
| `(5,3)` | -2 | 30 | -0.217005339961346 | -6.51016019884 |
| `(7,4)` | +1 | 56 | +0.0527808750147651 | +2.95572900083 |
| `(19,11)` | -2 | 418 | -0.0145179091920255 | -6.06848604227 |
| `(26,15)` | +1 | 780 | +0.00386292016315337 | +3.01307772726 |
| `(71,41)` | -2 | 5822 | -0.00103701130037872 | -6.03747979080 |
| `(97,56)` | +1 | 10864 | +0.000277726626463678 | +3.01722206990 |

The opposite-side ratio is itself a useful no-amplitude target. If the matching H4 shape is proportional to this ordinary Q4 torus fingerprint, then

\[
\boxed{
\frac{N_-^{21/8}\Delta M_-}
     {N_+^{21/8}\Delta M_+}
\to -2
}
\]

for matched generations of the `D=-2` and `D=+1` families, up to the common microscopic coupling and declared orientation convention.

At finite N the stronger conditional target uses the actual `E4(tau)` values rather than the asymptotic `-2` ratio.

## 7. Consequence for the H4-induced matching-root bias

The surviving square-torus empirical law is approximately

\[
\Delta M_{H4}\sim N^{-13/8}.
\]

The thermal slope scales as

\[
\bar M'\sim N^{3/8}.
\]

If the H4 shape factor has the simple `g2/E4` zero above, then on a fixed-D Pell family

\[
\Delta M_{H4}
\sim N^{-13/8}N^{-1}
=N^{-21/8}.
\]

The **H4-induced component** of the root shift then scales as

\[
\boxed{
\Delta p_{H4}\sim N^{-3}=L^{-6}
}
\]

instead of the square-torus `N^-2=L^-4` behavior.

The word “component” is essential. After H4 is suppressed, another scalar, higher-spin, composite or logarithmic sector can become the leading total root bias. A total `L^-6` law is a further prediction only after those surviving sectors are bounded.

## 8. Norm-5 H12 should be interpreted as a severe adversary, not a peer local operator

The current norm-5 preregistration wisely treats H12 with the same radial power as an empirical angular alias. The CFT spectrum shows that this is a deliberately harsh adversary.

### 8.1 Same thermal family

Starting from the thermal primary with `x_t=5/4`, a descendant of spin 12 requires chiral level difference 12. Since total descendant level is at least 12,

\[
\boxed{x_{\rm thermal,spin12}\ge 5/4+12=53/4.}
\]

So a thermal-family spin-12 field cannot have `x=21/4`.

### 8.2 Standard non-diagonal Potts/loop primaries

Critical loop models have non-diagonal fields with

\[
\Delta(r,s)=\frac14(\beta r-\beta^{-1}s)^2
-\frac14(\beta-\beta^{-1})^2,
\]

and spin `rs`. At percolation `Q=1`, one may take `beta^2=2/3`. The total scaling dimension is then

\[
 x(r,s)
 =-\frac1{12}+\frac{r^2}{3}+\frac{3s^2}{4}.
\]

For fixed absolute spin `m=|rs|`, AM-GM gives the universal lower bound

\[
\boxed{x\ge |m|-1/12.}
\]

Thus a standard local non-diagonal spin-12 primary obeys

\[
 x\ge 143/12\approx11.9167,
\]

with allowed half-integer-r states reaching `x=12` near the continuous optimum. Again this is nowhere near `21/4=5.25`.

Therefore:

> if norm-5 selects the **same-radial-exponent H12 sign law**, the result should not be narrated as “H12 replaces H4” inside the same ordinary local-CFT picture. It would be evidence that the simple local thermal-Q4 identification is wrong or incomplete.

Possible explanations would include a nonlocal/topological/defect sector, a finite-size mixture that masquerades as H12 on the tested lineages, or a more complicated lattice-to-CFT map.

This interpretation does **not** change the frozen Issue #57 scoring order or targets.

## 9. Research assessment and recommended order

My assessment after the N=185/265 result and the cross/either erratum is:

1. **The matching-odd orientation sector is no longer a discovery-stage signal.** It has enough independent support that the main question is operator identification and representation structure.
2. **The strongest missing theory bridge is now observable identification.** Exponent arithmetic alone is no longer the bottleneck. We need to show how the matching observable/topological projector couples to the thermal `Q4` torus sector.
3. **Norm 5 remains the right expensive next experiment.** But its H12 branch should be treated as a deep falsification branch, not as a cosmetically different harmonic fit.
4. **After norm 5 and the frozen 145->290 full curve, a moderate Pell-to-hexagonal pilot has unusually high information value.** `N=418` (`D=-2`) and `N=780` (`D=+1`) already have `N E4` values close to the asymptotic `-2:1` ratio and are dramatically cheaper than the gated N=1105 four-angle production.
5. **Do not spend the Pell experiment merely estimating an exponent.** Preserve the full curve and score the signed cross-family `E4` shape ratio, H4 suppression, root response, and any surviving H12/scalar sideband jointly.
6. **The q=2-versus-Jordan problem remains orthogonal.** The Gaussian `Q=2,5,10` plaquette/cocycle test from the isogeny program is still the cleanest way to decide whether the S-prime drift is an analytic correction or a non-semisimple logarithmic action.

## 10. Concrete validation tasks

### A. Algebraic verification — zero simulation

Implement an independent symbolic/`Fraction` checker of the coefficient `493/96` from:

- the Virasoro commutator;
- the level-2 null relation;
- the Brehm-Runkel `L_-4` torus Ward identity.

This is not a numerical fit.

### B. Observable bridge — theory first

For the actual matching observable, derive whether the first-order lattice anisotropy correction is:

- literally a `Q4` one-point insertion;
- an integrated correlator of the matching/topological projector with `Q4`;
- a defect-sector matrix element;
- or a logarithmic/Jordan extension.

The hexagonal zero is strongest if the matching observable is scalar under the elliptic-point stabilizer. This symmetry statement should be derived directly for the relevant torus homology combination.

### C. Pell geometry engine — tiny exact gate

Before production, for `(a,b)=(5,3),(7,4),(19,11),(26,15)` verify:

- period determinant and Smith invariants;
- homology/wrapping conventions;
- exact channel maps;
- orientation conventions;
- tiny enumeration where feasible.

### D. Frozen moderate-N pilot

If A-C succeed, the first informative pair is

```text
D=-2: (a,b)=(19,11), N=418
D=+1: (a,b)=(26,15), N=780
```

Score the no-fit geometry prediction before any free shape fit.

## Claim boundary

What is exact here:

- the reduction of the repository `Q4` one-point function to `(493/96) g2` in the ordinary `c=0,h=5/8` Virasoro module;
- the zero `g2(omega)=0`;
- the Pell `O(1/N)` approach and the asymptotic `N E4` constant;
- the lower dimension bound showing that a same-`x=21/4` local spin-12 primary is not available in the standard Potts/loop spectrum.

What remains conjectural:

- that the measured matching-odd residual couples dominantly to this ordinary `Q4` torus matrix element;
- that its full modulus dependence inherits exactly the same `g2/E4` factor;
- that no lower surviving sector masks the geometric annihilation;
- that total matching-root convergence becomes `L^-6` rather than merely the H4 component.

## References

- E. M. Brehm and I. Runkel, *Lattice models from CFT on surfaces with holes I: Torus partition function via two lattice cells*, J. Phys. A 55 (2022) 235001, arXiv:2112.01563. See the torus descendant Ward recursion and Eq. (4.52).
- P. Roux, S. Ribault and J. L. Jacobsen, *Torus one-point functions in critical loop models*, arXiv:2604.24491 (2026). See the Potts energy one-block solution and critical-loop spectrum.
- `notes/thermal-level4-spin4-candidate.md` and `scripts/virasoro_level4_candidate.py` in this repository.
