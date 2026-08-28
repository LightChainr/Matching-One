# Matching × rotation sector decomposition

Status: high-priority structural hypothesis prompted by the underpowered Pell matching-difference scan.

## 1. The key question is not only whether spin-4 exists

For a fixed wrapping/topological observable `R`, define on the same torus shape

\[
R_G(p,L,\theta)
\]

for the primary square-site lattice and

\[
R_{\hat G}(1-p,L,\theta)
\]

for its site-matching lattice.

The exact finite matching function is

\[
M_L(p,\theta)=R_G(p,L,\theta)-R_{\hat G}(1-p,L,\theta).
\]

Ordinary dimensionless critical observables can carry much larger finite-size corrections than `M` itself. In two-dimensional percolation, the `X_{t2}=4` field gives `L^-2`-type corrections, and Feng-Deng-Blote found a strongly orientation-dependent power-law amplitude on square lattices.

Yet the square-site matching root empirically shifts only near `L^-4`. Since

\[
M'_L(p_c)\sim L^{3/4},
\]

this root behavior corresponds to the very small residual

\[
M_L(p_c)\sim L^{-13/4}.
\]

Therefore any generic `L^-2` correction to the two individual wrapping probabilities must be absent from their **difference**.

This strongly motivates the hypothesis that the leading spin-4 correction is predominantly **matching-even** and cancels in `M`.

## 2. Four-channel parity decomposition

Take two orientations `theta_1,theta_2` at identical `N`, physical scale and modulus. Write

\[
R_{+1}=R_G(\theta_1),\quad
R_{+2}=R_G(\theta_2),
\]

\[
R_{-1}=R_{\hat G}(\theta_1),\quad
R_{-2}=R_{\hat G}(\theta_2),
\]

where probabilities are evaluated at complementary `p` and `1-p`.

Form four orthogonal linear sectors:

### matching-even, orientation-even

\[
E_{++}=\frac14(R_{+1}+R_{+2}+R_{-1}+R_{-2}).
\]

### matching-odd, orientation-even

\[
E_{-+}=\frac14(R_{+1}+R_{+2}-R_{-1}-R_{-2})
=\frac14[M(\theta_1)+M(\theta_2)].
\]

### matching-even, orientation-odd

\[
E_{+-}=\frac14(R_{+1}-R_{+2}+R_{-1}-R_{-2}).
\]

### matching-odd, orientation-odd

\[
E_{--}=\frac14(R_{+1}-R_{+2}-R_{-1}+R_{-2})
=\frac14[M(\theta_1)-M(\theta_2)].
\]

These channels distinguish two questions that the original axis/diamond matching-root scan conflated:

1. Is there an orientation-sensitive correction in the microscopic model?
2. Does matching exchange preserve or reverse that correction amplitude?

## 3. Strong working prediction

At critical complementary probabilities and fixed square-torus modulus, test

\[
E_{+-}\sim C_4L^{-2}
\]

while

\[
E_{--}=o(L^{-2}),
\]

possibly as small as the `L^-13/4` behavior inferred from the matching root.

In words:

> the leading spin-4 anisotropy is visible in the **sum** of primal and matching wrapping corrections but cancels from their **difference**.

If true, this is a structural explanation for why direct Monte Carlo of `M` is statistically difficult and why the matching estimator is unusually convergent.

## 4. Same-N Gaussian tomography makes the test clean

For periods `(a,b),(-b,a)`, every orientation at fixed `N=a^2+b^2` has the same

- site count;
- physical side length;
- torus modulus `tau=i`.

Therefore orientation-even scaling terms cancel directly in `R(theta_1)-R(theta_2)` without any uncertain `p_c` extrapolation.

For multiple angles, fit the two matching-parity combinations separately:

\[
S(\theta)=\frac12(R_G+R_{\hat G}),
\qquad
D(\theta)=\frac12(R_G-R_{\hat G})=\frac12M.
\]

Use

\[
S(\theta)=S_0+S_4\cos4\theta+S_8\cos8\theta+\cdots,
\]

\[
D(\theta)=D_0+D_4\cos4\theta+D_8\cos8\theta+\cdots.
\]

The leading hypothesis is

\[
S_4\sim L^{-2},
\qquad
D_4/S_4\to0.
\]

Do not assume the exact exponent of `D_4` until data resolve it.

## 5. Why this is statistically much easier than the old target

Directly detecting

\[
M(\theta_1)-M(\theta_2)
\]

may require resolving an `L^-13/4` signal. Detecting a single-lattice or matching-even orientation harmonic can instead expose an `L^-2` signal.

At physical size `L`, their ratio scales roughly as

\[
L^{13/4-2}=L^{5/4}.
\]

Thus at `L=8` the discovery signal can be larger by a factor about `13`, and at `L=17` by about `35`, before any covariance advantage.

This is the recommended route for establishing the angular physics. Only after the larger `S_4` sector is understood should production statistics be spent on the small matching-odd residue `D_4`.

## 6. Exact-control hierarchy

### Square bond at p=1/2

Use the same Gaussian tori and wrapping channels. This tests the expected `L^-2 cos4theta` anisotropy with no threshold uncertainty.

### Square site primary alone

At one frozen `p_ref` extremely near `p_c`, same-N orientation differences cancel the common leading thermal displacement. Measure the angular harmonic of `R_G`.

### Matching partner alone

At `1-p_ref`, repeat for `R_hat`.

### Only then form matching parity

Compare `S_4` and `D_4` and their size exponents.

This staged design can tell whether a null result in `M` is due to absence of spin-4 physics or due to successful matching cancellation.

## 7. A possible explanation of the L^-4 root

Suppose the leading ordinary corrections have the schematic form

\[
R_G=R_*+L^{-2}A_{2}^{(+)}(\theta)+\cdots,
\]

\[
R_{\hat G}=R_*+L^{-2}A_{2}^{(-)}(\theta)+\cdots.
\]

If matching forces

\[
A_{2}^{(+)}(\theta)=A_{2}^{(-)}(\theta),
\]

then the whole `L^-2` sector disappears from `M`.

If further lower-order matching-odd sectors vanish, the first surviving correction can be much smaller, producing the observed fast matching-root convergence.

This does not by itself derive the exponent `4`; it converts the empirical superconvergence into an operator-selection question:

> Which irrelevant/analytic sectors are matching-even and therefore projected out, and what is the first matching-odd sector allowed by the exact finite matching relation?

That is a sharper target for issue #12 than a generic list of correction exponents.

## 8. Acceptance criteria

Call the matching-even spin-4 hypothesis supported only if:

1. `S_4` is resolved and has a stable angular `cos4theta` pattern;
2. its size dependence is compatible with a substantially slower decay than `D_4`;
3. `D_4/S_4` decreases with size on held-out orientations/sizes;
4. exact square-bond controls validate the angular implementation;
5. the result persists across more than one wrapping definition or a prespecified optimal equal-mean channel.

A null matching difference by itself is not evidence: the corresponding even sector must be measured simultaneously.
