# Five scaling layers in the fixed-width charge transfer

2026-09-14.  Integration note.  The purpose is to keep several numerically visible exponents from being assigned to the same physical mechanism.

The first and fourth layers below use standard percolation CFT/thermal scaling.  The second and third are strong transfer diagnostics.  The fifth is the separate Kac-(4,2) conjecture.

## 1. Magnetic primary: x_m=5/48

At criticality the NN-safe and complementary-matching-safe sectors both propagate the magnetic/topological excitation.  Thus

\[
I^0_{G,w}(p_c)\sim\frac{2\pi x_m}{w},
\qquad x_m=5/48.                                               \tag{1.1}
\]

The transparent transfer gives

| `w` | `w I4` | `w I8` |
|---:|---:|---:|
| 4 | 0.6808677452 | 0.6677642268 |
| 5 | 0.6702211440 | 0.6643381257 |
| 6 | 0.6651785592 | 0.6620431213 |
| 7 | 0.6622847725 | 0.6604282152 |
| 8 | 0.6604615126 | 0.6592750557 |

against

\[
2\pi x_m=0.6544984694978736\ldots.                            \tag{1.2}
\]

This is the common primary sector and is the reason the eigenvalue identity can converge unusually fast.

## 2. Common sector-even irrelevant correction: candidate x=4

Average the two magnetic gaps:

\[
\bar I_w=\frac12\left[I^0_{4,w}(p_c)+I^0_{8,w}(1-p_c)\right].\tag{2.1}
\]

The diagnostic

\[
w^2\left[w\bar I_w-2\pi x_m\right]                         \tag{2.2}
\]

has values

```text
w=4  0.31708
w=5  0.31953
w=6  0.32805
w=7  0.33604
w=8  0.34367
```

so a correction

\[
I_w=\frac{2\pi x_m}{w}+C_{even}w^{-3}+\cdots                 \tag{2.3}
\]

is natural.  In CFT language a per-row correction `w^(1-x)` with exponent `-3` corresponds to a scalar irrelevant dimension

\[
\boxed{x_{even}=4.}                                          \tag{2.4}
\]

The obvious candidate is the usual identity-family/T-bar-T-type scalar correction.  The important point for the matching method is not its precise operator name: this lower-dimensional correction appears **common to both sectors** and therefore cancels to high order from their difference.

The `x=4` identification is a scaling interpretation, not an operator matrix-element proof.

## 3. Within-sector relaxation: Delta x=1

The second eigenvalue inside each safe magnetic kernel gives

\[
g_{G,w}=w\log(\lambda_1/|\lambda_2|).                        \tag{3.1}
\]

For `w=4,...,8`, the NN values decrease toward `2 pi` and the matching values increase toward `2 pi`.  This suggests

\[
\boxed{\Delta x_{relax}=1,}                                  \tag{3.2}
\]

consistent with a level-one descendant in the same magnetic module.

This gap controls convergence in the **longitudinal aspect ratio** and should not be confused with the irrelevant correction that shifts the pseudo-critical width sequence.

## 4. Thermal response: x_t=5/4

The charge-sector difference is

\[
\Theta_w(p)=I^0_{4,w}(p)-I^0_{8,w}(1-p).                     \tag{4.1}
\]

Differentiation in the thermal parameter couples to the percolation thermal field

\[
x_t=5/4,\qquad y_t=2-x_t=3/4.                                \tag{4.2}
\]

Hence an excitation energy per row has derivative

\[
\boxed{\Theta'_w\asymp w^{y_t-1}=w^{-1/4}.}                  \tag{4.3}
\]

The same power follows from the exact safe-sector pivotal-density identity plus the four-arm scaling mechanism.

## 5. Sector-odd correction: conjectural x=21/4

The difference at the critical point is far smaller than the average correction:

\[
\Theta_w(p_c)\asymp w^{-17/4}.                               \tag{5.1}
\]

A per-row correction `w^(1-x)` therefore points to

\[
\boxed{x_{odd}=21/4.}                                        \tag{5.2}
\]

The extended `c=0` Kac table contains exactly the scalar `(4,2)` value

\[
2h_{4,2}=21/4.                                                \tag{5.3}
\]

This motivates, but does not prove, the sector-odd Kac-(4,2) conjecture recorded separately.

Dividing the sector-odd mismatch by the thermal response gives

\[
\frac{w^{-17/4}}{w^{-1/4}}=w^{-4},                            \tag{5.4}
\]

which is the Jacobsen pseudo-critical shift.

## 6. Why lower irrelevant dimensions do not contradict the n^-4 root shift

A common source of confusion is to say: if an `x=4` irrelevant field exists, why does the pseudo-critical root not shift with exponent

\[
x=4\quad\Rightarrow\quad 4-x_t=11/4 ?                        \tag{6.1}
\]

The answer is that the root uses the **difference** of the primal and dual magnetic-sector energies.  A sector-even correction can be large in each energy and still disappear from

\[
I^0_4-I^0_8.                                                   \tag{6.2}
\]

The transfer data show exactly this hierarchy:

```text
individual/common correction  >>  sector difference.
```

Thus the relevant question for the matching root is not “what is the leading irrelevant field of percolation?” but

> what is the lowest-dimensional field whose matrix element is ODD under exchange of the two topological magnetic sectors?

The Kac-(4,2) conjecture is a proposed answer to that narrower question.

## 7. Summary table

| phenomenon | scale in per-row transfer | proposed dimension/spacing | role |
|---|---:|---:|---|
| magnetic primary | `w^-1` | `x_m=5/48` | common leading sector |
| common irrelevant | `w^-3` | `x=4` | cancels in sector difference |
| longitudinal relaxation | gap `2pi/w` | `Delta x=1` | aspect convergence |
| thermal derivative | `w^-1/4` | `x_t=5/4` | root susceptibility |
| sector-odd mismatch | `w^-17/4` | conjectural `x=21/4` | produces `w^-4` root shift |

This table should be used when interpreting future transfer fits; assigning all five scales to one “irrelevant exponent” would erase the mechanism.

## 8. Claim boundary

The finite transfer numbers are direct controls.  The primary magnetic interpretation is the same CFT argument already used by Jacobsen.  `x=4`, `Delta x=1`, and `x=21/4` are increasingly interpretive statements: the first two are strongly suggested by direct scaling diagnostics, while the Kac-(4,2) sector-odd identification remains a targeted conjecture requiring a sector-selection proof.
