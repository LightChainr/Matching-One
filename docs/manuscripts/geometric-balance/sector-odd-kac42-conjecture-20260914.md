# A sector-odd Kac-(4,2) conjecture for the n^-4 charge crossing

2026-09-14.  CONJECTURAL CFT interpretation of the fixed-width charge-transfer data.  The numerical eigenvalue sequence itself is not new: the charge-coexistence roots reproduced by `fixed_width_charge_transfer.py` coincide with the `n x infinity` square-site eigenvalue-identity sequence of Jacobsen (2015).  The new question here is why the primal/dual **sector difference** starts at the observed power.

## 1. Probability/free-energy object

For circumference `w`, define

\[
\Theta_w(p)=I^0_{4,w}(p)-I^0_{8,w}(1-p),                     \tag{1.1}
\]

where `I^0=-log lambda^0` is the per-row free-energy cost of remaining in the horizontally trivial homology sector.  Its zero is the fixed-width charge-coexistence point

\[
p_w^{ch}:\quad\Theta_w(p_w^{ch})=0.                           \tag{1.2}
\]

Digital Alexander duality identifies this criterion with equality of the primal/open and dual/closed topological transfer sectors.  Numerically `p_w^{ch}` is exactly the Jacobsen semi-infinite-cylinder sequence to the displayed precision through the checked widths.

Jacobsen's CFT explanation is that the two sectors both propagate the magnetic excitation and therefore have the same leading scaling dimension

\[
x_m=5/48.                                                      \tag{1.3}
\]

For square geometry the per-row excitation cost should therefore satisfy

\[
w I^0_{G,w}(p_c)\longrightarrow 2\pi x_m
=0.6544984694978736\ldots                                    \tag{1.4}
\]

for both sectors.  The safe-transfer control gives

| `w` | `w I^0_4(p_c)` | `w I^0_8(1-p_c)` |
|---:|---:|---:|
| 4 | 0.6808677452 | 0.6677642268 |
| 5 | 0.6702211440 | 0.6643381257 |
| 6 | 0.6651785592 | 0.6620431213 |
| 7 | 0.6622847725 | 0.6604282152 |
| 8 | 0.6604615126 | 0.6592750557 |

The common magnetic gap is visible directly; the much smaller difference is the object of this note.

## 2. Three observed powers

Using only the repository reference `p_c=0.59274605079` as a scaling diagnostic, not as an input to the charge roots, the transparent safe transfer gives

\[
\Theta_w(p_c) w^{17/4}
=1.185995,1.099646,1.059961,1.035803,1.021632                \tag{2.1}
\]

for `w=4,...,8`, while

\[
\Theta_w'(p_w^{ch}) w^{1/4}
=3.485121,3.447892,3.426227,3.412433,3.403145.                \tag{2.2}
\]

Finally

\[
(p_c-p_w^{ch})w^4
=0.340193,0.318892,0.309348,0.303528,0.300197.                \tag{2.3}
\]

The powers are mutually consistent rather than three independent fits:

\[
w^{-17/4}/w^{-1/4}=w^{-4}.                                   \tag{2.4}
\]

Jacobsen independently found the pseudo-critical correction sequence `Delta_1=4`, followed by corrections compatible with `6,8,...`.  The present transfer gives a probability/topology decomposition of the first exponent into a critical sector mismatch and a thermal response.

## 3. Finite-size scaling dictionary

The percolation thermal scaling dimension is

\[
x_t=5/4,\qquad y_t=2-x_t=3/4.                                \tag{3.1}
\]

A cylinder excitation energy per transfer step has the standard scaling form

\[
I_w(p)=w^{-1}\mathcal X((p-p_c)w^{y_t},\{u_iw^{2-x_i}\}).     \tag{3.2}
\]

Therefore

\[
\partial_p I_w(p_c)\asymp w^{y_t-1}=w^{-1/4}.                \tag{3.3}
\]

If the first correction whose matrix element differs between the primal and dual magnetic sectors is a scalar irrelevant field of dimension `x_odd`, then the common lower corrections cancel in the difference and

\[
\Theta_w(p_c)\asymp w^{1-x_{odd}}.                            \tag{3.4}
\]

Linearizing the zero of `Theta` gives

\[
p_w^{ch}-p_c
\asymp w^{-(x_{odd}-x_t)}.                                   \tag{3.5}
\]

The observed exponent four therefore demands

\[
\boxed{x_{odd}=x_t+4=21/4.}                                  \tag{3.6}
\]

This is exactly the exponent independently suggested by (2.1).

## 4. The Kac-(4,2) candidate

For the extended `c=0` percolation Kac table,

\[
h_{r,s}=\frac{(3r-2s)^2-1}{24}.                              \tag{4.1}
\]

Thus

\[
h_{4,2}=21/8,\qquad x_{4,2}=h+\bar h=21/4                  \tag{4.2}
\]

for the scalar diagonal field.

This gives the concrete conjecture:

> **Kac-(4,2) sector-odd conjecture.**  In the square-site primal/dual magnetic cylinder sectors, all lower-dimensional correction fields have equal amplitudes and cancel from `I^0_4-I^0_8`; the first scalar field with nonzero sector-odd matrix element belongs to the `(4,2)` conformal family, with `x=21/4`.

Then (3.4)--(3.6) give the observed `w^-17/4` charge mismatch and `w^-4` pseudo-critical shift.

## 5. Why the full 4,6,8,... sequence becomes natural

Scalar descendants in the same family raise the scaling dimension by even integers,

\[
x_{4,2}^{(j)}=21/4+2j,\qquad j=0,1,2,\ldots.                 \tag{5.1}
\]

Dividing their sector-difference corrections by the same thermal response yields pseudo-critical shift exponents

\[
\boxed{\Delta_j=x_{4,2}^{(j)}-x_t=4+2j,}                     \tag{5.2}
\]

namely

\[
4,6,8,10,\ldots.                                              \tag{5.3}
\]

This reproduces the empirical correction ladder reported by Jacobsen without introducing a separate unrelated irrelevant field for each power.

## 6. What must be proved or falsified

The exponent arithmetic alone is not a field identification.  Two substantive statements remain.

1. **Sector content.** The `(4,2)` family must actually occur in the relevant topological/magnetic cylinder sector with square-lattice scalar symmetry.
2. **Selection/cancellation.** Every lower-dimensional irrelevant family allowed by bulk symmetry must have equal matrix elements in the primal and dual sectors, or vanish separately, so that it drops out of the sector difference.

Either failure kills the proposed identification even if a short-width power fit looks good.

The cleanest numerical/transfer tests are therefore not another fit of `p_w`.  They are:

- compute `Theta_w(p_c)` itself at substantially larger `w` and test the direct exponent `17/4`;
- compute the two sector eigenvectors and project explicit lattice irrelevant perturbations into their difference;
- use a lattice with different rotational symmetry.  A genuine sector-odd CFT family should persist with changed amplitudes subject to symmetry selection, whereas a purely square-lattice analytic artifact need not.

## 7. Literature boundary

Jacobsen, *Critical points of Potts and O(N) models from eigenvalue identities in periodic Temperley-Lieb algebras*, arXiv:1507.03027, establishes the semi-infinite-cylinder eigenvalue criterion, explains that the open/closed sectors share the magnetic exponent, and reports pseudo-critical corrections `4,6,8,...`.  The paper explicitly notes that the shared leading exponent alone does not derive the `n^-4` law.

The general `c=0` logarithmic percolation CFT/Kac framework is standard; see, for example, Mathieu--Ridout, arXiv:0708.0802, and the subsequent percolation LCFT literature.  I did not find a source identifying Jacobsen's `Delta_1=4` specifically with the scalar `(4,2)` family.  Accordingly sections 4--6 are a research conjecture, not a literature claim.
