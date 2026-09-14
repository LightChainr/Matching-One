# Integrable Potts field theory route to the universal charge scaling function

2026-09-14.  Literature-grounded continuum programme for #767, updated after identifying the critical safe sector in the periodic dilute Temperley--Lieb representation.

The problem is now narrower than in the first version of this note.  The UV magnetic twist is no longer arbitrary: critical triangular-site percolation identifies the safe/no-noncontractible-loop sector with the zero-defect periodic dilute-TL module

\[
\boxed{W_{N,0,\omega=i},\qquad \alpha=\omega+\omega^{-1}=0.} \tag{0.1}
\]

The remaining problem is to construct or locate its **off-critical thermal/massive finite-volume continuation** at `q=1`.

## 1. The lattice transfer has identified one continuum magnetic object

The fixed-width square-site safe transfer gives

\[
I^0_{4,w}(p),\qquad I^0_{8,w}(1-p).                           \tag{1.1}
\]

At criticality both converge to the same magnetic primary gap

\[
wI^0\to2\pi x_m,
\qquad x_m=5/48.                                              \tag{1.2}
\]

Near criticality the natural universal organization is one magnetic finite-size energy function sampled at opposite thermal signs:

\[
\boxed{
\mathcal F(X)
=\mathcal E_m(X)-\mathcal E_m(-X).}                          \tag{1.3}
\]

Hence

\[
\mathcal F(-X)=-\mathcal F(X).                               \tag{1.4}
\]

The square NN and complementary matching kernels are two microscopic regularizations of this same magnetic/topological sector.  Their dual-odd difference is the charge free energy.

## 2. The critical periodic dilute-TL sector is now explicit

Morin-Duchesne, Kluemper and Pearce describe critical triangular-site percolation by the Yang--Baxter-solvable dilute `A_2^(2)` model.  Its periodic standard modules are labelled by defect number `d` and twist

\[
\omega=e^{i\gamma}.                                          \tag{2.1}
\]

For the periodic ground states they obtain

\[
(h,\bar h)
=
\left(
\Delta_{\gamma/\pi,d/2},
\Delta_{\gamma/\pi,-d/2}
\right),                                                     \tag{2.2}
\]

where

\[
\Delta_{r,s}=\frac{(3r-2s)^2-1}{24}.                         \tag{2.3}
\]

In the zero-defect module, the noncontractible-loop weight is

\[
\alpha=\omega+\omega^{-1}=2\cos\gamma.                       \tag{2.4}
\]

The safe transfer rejects a state precisely when a horizontal noncontractible occupied loop closes.  The loop-language specialization is therefore

\[
\boxed{\alpha=0,
\qquad \omega=\pm i,
\qquad d=0.}                                                   \tag{2.5}
\]

Choose `omega=i`, so `gamma=pi/2`.  Then

\[
h=\bar h=\Delta_{1/2,0}=5/96,
\qquad x=5/48.                                                \tag{2.6}
\]

Thus the magnetic UV fingerprint is not just inferred from a fitted gap: it is the ground state of the published `W_{N,0,i}` standard module.

Primary source: A. Morin-Duchesne, A. Kluemper, P. A. Pearce, *Critical site percolation on the triangular lattice: From integrability to conformal partition functions*, arXiv:2211.12379v2.

The detailed lattice/transfer interpretation and square-site checks are recorded in `safe-transfer-pdtl-magnetic-sector-20260914.md`.

## 3. Four independent checks of the sector dictionary

The proposed UV dictionary simultaneously explains facts that were obtained independently in the square-site analysis.

### 3.1 Noncontractible-loop semantics

`alpha=0` kills a noncontractible loop exactly, matching the safe transfer's rejection rule.

### 3.2 State-space dimension

The periodic zero-defect dilute-TL module has central-trinomial dimension.  The transparent safe automaton has exactly

\[
1,3,7,19,51,141,393,1107,3139,\ldots                        \tag{3.1}
\]

states.

### 3.3 Magnetic ground gap

Equation (2.6) gives `x_m=5/48`, matching

\[
wI^0_w(p_c)\to2\pi(5/48).                                   \tag{3.2}
\]

### 3.4 Level-one descendants

The first descendants of a scalar magnetic primary have `Delta x=1` and spins `+/-1`.  Direct square-site diagonalization finds a twofold first excited eigenvalue with

\[
w\log(\lambda_0/|\lambda_1|)\to2\pi,                        \tag{3.3}
\]

while the one-column translation eigenvalues on this two-dimensional space are exactly

\[
e^{\pm2\pi i/w}                                               \tag{3.4}
\]

to machine precision for `w=5,...,9`.

See `safe-transfer-momentum-spectrum-w5-w9-20260914.json`.

These checks sharply reduce the risk that `W_{0,i}` is only a numerical weight coincidence.

## 4. Integrable massive Potts theory remains the natural off-critical framework

The scaling `q`-state Potts theory for `q<=4` under thermal perturbation is integrable.  Chim--Zamolodchikov give the kink scattering theory; Dorey--Pocklington--Tateo develop finite-size TBA/NLIE equations for thermal Potts flows and emphasize continuous `q` formulations relevant near `q=1`.  Delfino--Cardy demonstrate that the `q->1` massive Potts continuation yields nontrivial percolation amplitudes and form factors.

What the material inspected here does **not** yet hand us is a ready-made finite-volume equation explicitly carrying the `d=0`, `omega=i`, `alpha=0` magnetic topological sector through the thermal perturbation.

The target is therefore no longer “find some twist with UV `c_eff=-5/4`.”  It is:

> Continue the specific critical module `W_{0,i}` into the massive thermally perturbed Potts theory and compute its finite-volume energy on the two signs of the thermal coupling.

The 2002 Dorey--Pocklington--Tateo finite-size equations remain a plausible continuum framework, but the sector insertion/source implementing (2.5) must be derived or found.

## 5. The UV effective-central-charge check is retained

For a sector with lowest total dimension `x`,

\[
c_{eff}=c-12x.                                                \tag{5.1}
\]

Percolation has `c=0`, so the magnetic module requires

\[
\boxed{c_{eff}^{(m)}=-12(5/48)=-5/4.}                        \tag{5.2}
\]

Any proposed off-critical twisted NLIE must return (5.2) in its ultraviolet limit **and** its lattice/topological twist must reduce to `alpha=0` / `omega=i`.  Matching only one of these checks is insufficient.

## 6. Desired massive finite-volume output

Let

\[
r=m_{phys}R                                                   \tag{6.1}
\]

be the usual dimensionless massive circumference.  We seek

\[
\mathcal E_{W_{0,i}}^{(+)}(r),
\qquad
\mathcal E_{W_{0,i}}^{(-)}(r),                               \tag{6.2}
\]

on the two signs of the thermal perturbation, in a common mass normalization.

Potts duality should relate the two branches.  The universal charge function is

\[
\boxed{
\mathcal F(r)
=\mathcal E_{W_{0,i}}^{(+)}(r)
 -\mathcal E_{W_{0,i}}^{(-)}(r).}                            \tag{6.3}
\]

After fixing the lattice thermal metric,

\[
X\propto(p-p_c)w^{3/4},                                      \tag{6.4}
\]

(6.3) should be the continuum limit of

\[
w\Theta_w(p).                                                \tag{6.5}
\]

## 7. Lattice data now give an entire target curve, not only derivatives

`dual-odd-thermal-metric-20260914.md` defines

\[
\mathcal F_w(X)
=w\Theta_w(h_w+Xw^{-3/4}).                                   \tag{7.1}
\]

The curves for `w=4,...,9` already collapse well over `|X|<=1.5`.

Their leading symmetric contamination scales as `w^-3/4` and is quantitatively consistent with a quadratic analytic thermal-field reparameterization.  After removing that coordinate effect, the odd part is the appropriate lattice target for an integrable calculation.

This improves the comparison protocol:

1. determine the analytic lattice thermal metric from the even residual;
2. extract the antisymmetric continuum curve;
3. compare the full curve, not just `F'(0)`;
4. use the first/third/fifth derivative ratios as local checks.

The current raw-logit local fit predicts a fifth derivative of order

\[
w^{-11/4}\Theta_h^{(5)}\sim-1.5\times10^{-2},                \tag{7.2}
\]

before final metric calibration.

## 8. Mandatory checks on any massive `W_{0,i}` construction

A successful TBA/NLIE or other exact construction should satisfy:

1. **UV module:** `d=0`, `omega=i`, `alpha=0` and `x=5/48`.
2. **UV effective central charge:** `c_eff=-5/4`.
3. **First descendants:** level-one spin `+/-1` above the magnetic primary.
4. **Dual oddness:** `F(-X)=-F(X)` in the correctly normalized thermal coordinate.
5. **Thermal slope:** reproduce the safe pivotal/Perron amplitude after metric fixing.
6. **Nonlinear curve:** agree with the transfer collapse over a finite `X` interval, not just at `X=0`.
7. **Infrared sectors:** match the appropriate kink/order/disorder excitation on the two thermal signs.

The first three are now unusually rigid UV constraints.

## 9. Square-site versus triangular-site scope

For critical triangular-site percolation, `W_{N,0,i}` belongs to the explicit integrable periodic dilute-TL lattice solution.

For square-site NN/matching percolation, we do **not** claim a Yang--Baxter conjugacy.  The square Bernoulli safe transfer appears to act on the same annular/topological module and has the same continuum sector, but its local weights are a nonintegrable regularization.

This is enough for universality-based continuum matching; it is not enough to import finite-width Bethe roots from the triangular model.

## 10. Sharpened unresolved problem

The previous version of this note asked for a generic-q magnetic twist.  The critical twist has now been identified.

The unresolved task is specifically:

\[
\boxed{
\text{massive / thermal finite-volume continuation of }W_{0,i}
\text{ at }q=1.}                                             \tag{10.1}
\]

Promising routes include:

- insert the corresponding noncontractible-loop/twist weight into a continuous-q Potts finite-volume NLIE;
- derive the excited/twisted solution by analytic continuation from the ground-state equation;
- construct an off-critical staggered/dilute-`A_2^(2)` realization whose scaling limit is the thermal Potts perturbation and retain `omega=i`.

A literature search in this pass found the continuous-q thermal Potts finite-size framework and the critical periodic dilute-TL twist, but not a ready-made equation already combining both in the required percolation magnetic sector.

## 11. Claim boundary

The critical periodic dilute-TL standard-module formula, twist label and conformal weights are literature facts.  `alpha=0 <-> omega=+/-i` is exact algebra.  The safe-transfer state-count, momentum and gap checks are deterministic lattice calculations.

The assertion that the square-site near-critical charge curve converges to the massive finite-volume `W_{0,i}` Potts energy difference is a universality/duality synthesis.  The explicit massive `q=1` twisted NLIE remains to be derived or located.
