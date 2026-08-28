# Elliptic-point spin filter for modular-scalar matching observables

**Status:** structural CFT/topology hypothesis with an exact modular-selection core. This note is deliberately broader than the ordinary-Q4 one-point calculation in `thermal-q4-torus-ward-identity.md`.

The key question is whether the predicted hexagonal H4 zero depends on the matching observable being literally a thermal `Q4` one-point function. It need not.

If the relevant wrapping/homology observable is a **modular scalar**, then the first-order response to a pure spin-`s` irrelevant perturbation obeys an elliptic-fixed-point selection rule. At the hexagonal torus this kills spin 4 and spin 8 but permits spin 12. This gives an H4/H12 filter that is more representation-theoretic than the `E4` ansatz.

## 1. Why homology-class observables are the right scalar variables

On a torus, an FK/percolation configuration determines the image of its first homology in

```text
H1(T^2,Z) = Z x Z.
```

Arguin's torus Potts/FK homology calculation generalizes Pinson's percolation result and explicitly organizes probabilities by the resulting subgroup. The total rank-0, rank-2, and suitable sums over primitive rank-1 sectors are modular invariant: an `SL(2,Z)` basis change relabels primitive winding classes but cannot change homology rank.

This is exactly the distinction the repository should preserve:

- a named **horizontal** winding channel is vector-valued under modular basis changes;
- `cross` / rank-2 wrapping is scalar;
- `no nontrivial winding` / rank-0 is scalar;
- `any nontrivial winding` after summing primitive classes is scalar;
- some `both/either` implementation names require checking their exact topological definition before calling them scalar.

The #134 cross/either erratum makes this typed distinction operationally important as well as theoretical.

## 2. Linear response of a scalar torus observable

Let `O_top` be a dimensionless torus observable that is scalar under modular changes of basis. Perturb the critical theory by a bulk scaling field `Phi_(h,hbar)` of spin

\[
s=h-\bar h
\]

with coupling `lambda_s`. To first order,

\[
\delta\langle O_{top}\rangle
\propto
\lambda_s
\int_{T_\tau} d^2z\,
\langle O_{top}\,\Phi_{h,\bar h}(z)\rangle_\tau.
\]

After the overall area/length power has been factored out, call the remaining shape coefficient `A_s(tau)`.

Under

\[
\gamma=\begin{pmatrix}a&b\\c&d\end{pmatrix}
\in PSL(2,Z),
\qquad
\tau\mapsto\gamma\tau=\frac{a\tau+b}{c\tau+d},
\]

the integrated insertion has modular weights shifted by the measure, but its phase is still controlled by the same spin `s`. At a fixed point `tau_*` for which the induced coordinate map is an isometry, `|c tau_*+d|=1`, the scalar response must satisfy schematically

\[
A_s(\tau_*)
=e^{i s\varphi_\gamma} A_s(\tau_*).
\]

Therefore

\[
\boxed{
e^{i s\varphi_\gamma}\ne1
\quad\Longrightarrow\quad
A_s(\tau_*)=0.
}
\]

This is the same fixed-point logic that makes nontrivial-weight modular forms vanish at elliptic points, but it applies to any modular-scalar observable whose linear response transforms homogeneously in the spin-`s` sector.

## 3. Square and hexagonal elliptic points give complementary filters

### Square point `tau=i`

The `S` stabilizer acts by a quarter turn. Its local phase is

\[
\varphi_S=\pi/2.
\]

A scalar response can therefore survive only when

\[
e^{i s\pi/2}=1,
\qquad
s\in4\mathbb Z.
\]

So the square torus is exactly where H4, H8, H12, ... are all compatible with the stabilizer. This explains why fixed-`tau=i` orientation tomography alone can suffer harmonic aliasing.

### Hexagonal point `rho=e^{2\pi i/3}` or `omega=rho+1`

The `ST` stabilizer of `rho` has

\[
c\rho+d=\rho+1=e^{i\pi/3},
\]

so the local coordinate is rotated by `pi/3` up to orientation convention. A scalar response can survive only when

\[
e^{i s\pi/3}=1,
\qquad
s\in6\mathbb Z.
\]

Combining this with the square-lattice harmonic restriction `s in 4 Z` gives

\[
\boxed{s\in12\mathbb Z.}
\]

Thus for a square-lattice anisotropy measured through a modular-scalar torus observable:

```text
spin 4:   killed at the hexagonal elliptic point
spin 8:   killed at the hexagonal elliptic point
spin 12:  allowed
spin 16:  killed
spin 20:  killed
spin 24:  allowed
...
```

This is a much cleaner reason to use the hexagonal point as an H4/H12 discriminator than merely noting that `E4(omega)=0`.

## 4. Relation to the exact Q4 Ward identity

For the ordinary thermal `Q4` candidate, `thermal-q4-torus-ward-identity.md` derives the stronger explicit formula

\[
\frac{\langle Q_4\phi\rangle}{\langle\phi\rangle}
=\frac{493}{96}g_2(\tau).
\]

Since `g2(omega)=0`, that calculation realizes the general stabilizer selection rule explicitly.

But the stabilizer argument says the zero can survive even when the matching observable is not literally the identity background for a Q4 one-point function. What is required is weaker:

1. the measured topological observable is modular scalar;
2. the leading anisotropic correction transforms homogeneously as spin 4;
3. no singular scalar prefactor cancels the elliptic zero;
4. the relevant finite-size object is not secretly a vector-valued directional channel.

This is the theory bridge that should now be derived for the matching function.

## 5. Logarithmic/Jordan mixing does not automatically evade the zero

A rank-2 logarithmic pair of the **same spin** may transform under the elliptic stabilizer by a Jordan block

\[
R=\lambda(I+N),
\qquad
N^2=0,
\qquad
\lambda=e^{is\varphi}.
\]

If `lambda != 1`, then `I-R` is invertible. A homogeneous scalar fixed-point condition therefore still forces the whole Jordan block to vanish.

So logarithmic mixing by itself is not a loophole for spin 4 at the hexagonal point.

Possible loopholes are more specific:

- the topological observable belongs to a nontrivial vector-valued modular representation;
- the response contains an inhomogeneous/quasimodular connection term mixing with a different sector;
- the lattice anisotropy couples to several spins and a surviving spin-12 component dominates;
- the continuum limit of the chosen channel is not scalar under the elliptic stabilizer.

These alternatives are scientifically useful because they can be distinguished.

## 6. Pell approach: why the extra `1/N` is symmetry-driven

The integer square-lattice Pell periods

```text
u=(2b,0)
v=(b,a)
a^2-3b^2=D
N=2ab
```

approach the hexagonal point with

\[
\tau_N-\omega=O(N^{-1}).
\]

If `A_4(tau)` is analytic and has a symmetry-forced simple zero at `omega`, then generically

\[
A_4(\tau_N)
=A_4'(\omega)(\tau_N-\omega)+\cdots
=O(N^{-1}).
\]

This is enough to predict the extra one power of `N` suppression even without identifying `A_4` with `E4` exactly.

For the observed square-torus H4 law

\[
\Delta M\sim N^{-13/8},
\]

the generic simple-zero prediction is therefore

\[
\boxed{\Delta M_{H4, Pell}\sim N^{-21/8}.}
\]

The H4-induced root component then scales as `N^-3=L^-6` after dividing by the `N^(3/8)` thermal slope.

The stronger `-2` cross-family target follows if the zero is simple and the first derivative is common to the two fixed-defect Pell families: their signed distances satisfy `N(tau-omega)->iD/2` with `D=-2,+1`.

## 7. What should be proved before a Pell production run

### A. Exact channel scalarity

Using the canonical wrapping-channel algebra (#146), classify every candidate matching representation by its `SL(2,Z)` action:

```text
trivial scalar
finite permutation/vector representation
orientation-dependent component
```

The Pell spin-filter test should use a scalar representation, preferably the rank-2/cross homology event if it gives the same matching function under the exact identity.

### B. Matching-pair action

Show that the primal/matching combination used for `DeltaM` remains scalar under the same torus basis change after the complement/matching map. This is weaker than proving an OPE automorphism.

### C. Linear-response transformation

Derive the modular covariance of the first anisotropic correction coefficient. The target is a statement of the form

\[
A_4(\gamma\tau)
=(c\tau+d)^p(c\bar\tau+d)^q A_4(\tau)
\]

with `p-q=4`, possibly vector-valued in a declared topological sector. At the elliptic fixed point only the representation eigenvalue matters for the zero.

### D. Tiny finite-quotient controls

For the first Pell period matrices, verify that apparent sign changes are not caused by channel conventions, Smith-group changes, or orientation canonicalization.

## 8. Decision impact

If norm 5 supports H4, the modular-scalar elliptic filter becomes a high-value independent operator test because it attacks the same spin assignment with a different group action.

If norm 5 instead supports the same-radial H12 adversary, the hexagonal point becomes even more informative: H12 is precisely the first square harmonic that the elliptic stabilizer permits. A nonvanishing residual at the Pell limit would then be qualitatively expected.

If both H4 and H12 mixtures are present, the hexagonal point acts as an **improved geometry** that annihilates the leading H4 piece and enhances relative sensitivity to H12 and scalar sectors.

## References

- L.-P. Arguin, *Homology of Fortuin--Kasteleyn clusters of Potts models on the torus*, J. Stat. Phys. 109 (2002) 301--310, arXiv:hep-th/0111193. The homology-subgroup probabilities are modular invariant.
- P. Roux, S. Ribault and J. L. Jacobsen, *Torus one-point functions in critical loop models*, arXiv:2604.24491 (2026), for modular covariance of critical-loop torus one-point functions.
- E. M. Brehm and I. Runkel, arXiv:2112.01563, for torus descendant Ward identities.
- Repository issues #103, #114, #145, #146 and PR #151.
