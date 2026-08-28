# Thermal level-4 spin-4 candidate at c=0

The current same-N data motivate testing a correction with

\[
\Delta M\propto \cos(4\theta)L^{-13/4}.
\]

Because percolation has `y_t=1/nu=3/4`, the bulk thermal field has

\[
x_t=2-y_t=5/4,\qquad h=\bar h=5/8.
\]

At `c=0,h=5/8` the chiral module is level-2 degenerate:

\[
\chi_2=(L_{-2}-\tfrac23L_{-1}^2)|h\rangle=0.
\]

The key algebraic question is whether quotienting this null submodule removes every level-4 spin-4 descendant. It does not.

In the level-4 PBW basis

`[L_-4, L_-3 L_-1, L_-2^2, L_-2 L_-1^2, L_-1^4]`,

the quotient contains one non-null quasiprimary direction. A convenient representative is

\[
\boxed{Q_4=(40L_{-2}^2-60L_{-3}L_{-1}-9L_{-4})|h\rangle},
\]

with

\[
L_1Q_4=0.
\]

It is not in the span of `L_-2 chi2` and `L_-1^2 chi2`. The exact rational check is in `scripts/virasoro_level4_candidate.py`.

Tensoring with the antiholomorphic thermal primary gives weights

\[
(37/8,5/8)\quad\text{or}\quad(5/8,37/8),
\]

hence

\[
\boxed{x=21/4,\qquad s=\pm4.}
\]

Thus the candidate required by the observed scaling exists algebraically; `5/4+4` is not killed by the level-2 null relation.

## Working identification hypothesis

The root-moving matching-odd, rotation-spin-4, thermal-even sector may contain

\[
D_4^{\rm even}(L,0)=L^{-13/4}(A+B\log L+\cdots).
\]

Since `M'_L(pc)~L^(3/4)`, this gives

\[
p_L^*-p_c\sim L^{-4}
\]

(up to logarithmic/subleading corrections).

This is a candidate identification, not a result. Bulk percolation is a `c=0` LCFT, and the physical energy/Kac operator can participate in zero-norm/logarithmic multiplets. Matching parity also does not follow from Virasoro algebra; it must be measured.

## Required numerical discriminator

After threshold-rank Newman-Ziff reconstructs the full curve, project independently onto:

1. matching even/odd;
2. rotation spin 4;
3. thermal even/odd;
4. size exponent.

For two same-N orientations define

\[
P_4[X]=\frac{X(\theta_1)-X(\theta_2)}{\cos4\theta_1-\cos4\theta_2}.
\]

Define intrinsic thermal reflection points using the direction-averaged matching function:

\[
\bar M(p_-)=-u,\qquad \bar M(p_+)=+u.
\]

Then measure

\[
X_4^{even}(u)=\tfrac12[P_4X(p_+)+P_4X(p_-)],
\]

\[
X_4^{odd}(u)=\tfrac12[P_4X(p_+)-P_4X(p_-)].
\]

Compare held-out predictions for `L^-13/4`, `L^-13/4(A+B log L)`, and a free exponent. This avoids using disputed last digits of `pc` to define thermal parity.

## Magic-angle design correction

`N=169: (13,0) vs (12,5)` has equal area and modulus but not the same finite translation group: the axis quotient is `Z_13 x Z_13`, while primitive `(12,5)` is cyclic `Z_169`. It is not a canonical shared-cyclic-field CRN pair.

Cleaner primitive same-N confirmation candidates are:

- `N=565`: `(23,6)` vs `(22,9)`, with the latter near `pi/8` (`cos4theta ~ 0.0175`);
- `N=985`: `(27,16)` vs `(29,12)`, with `(29,12)` extremely near the spin-4 zero (`cos4theta ~ 0.00144`).

These are confirmation sizes. The immediate priority remains full-curve threshold ranks at `N=65,85,145`.

Background: bulk `c=0` percolation is logarithmic (Vasseur et al., arXiv:1110.1327). Modern work on Kac/energy operators at `c=0` emphasizes zero-norm/logarithmic multiplets (Y. He, arXiv:2411.18696).
