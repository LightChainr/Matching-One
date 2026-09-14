# A two-free-energy phase diagram for the homology rank at fixed p

2026-09-14.  Exact digital Alexander duality plus the arbitrary-shape fixed-`p` homological-free-energy theorem reduce the asymptotic rank law to two complementary scalar free energies.

## 1. Complementary correlation-norm minima

For an honest period lattice `Lambda_n` with `N_n` sites, fix black NN density `p` and complementary white matching density

\[
q=1-p.                                                        \tag{1.1}
\]

Define

\[
\rho_{4,n}(p)
=\min_{\lambda\in\Lambda_n\setminus0}\tau_{4,p}(\lambda),   \tag{1.2}
\]

\[
\rho_{8,n}(q)
=\min_{\lambda\in\Lambda_n\setminus0}\tau_{8,q}(\lambda).   \tag{1.3}
\]

Assume both minima diverge.  Introduce the opportunity ratios

\[
\alpha_{4,n}(p)=\frac{\log N_n}{\rho_{4,n}(p)},
\qquad
\alpha_{8,n}(q)=\frac{\log N_n}{\rho_{8,n}(q)}.               \tag{1.4}
\]

Equivalently define signed free energies

\[
\Phi_{4,n}(p)=\rho_{4,n}(p)-\log N_n,
\qquad
\Phi_{8,n}(q)=\rho_{8,n}(q)-\log N_n.                         \tag{1.5}
\]

## 2. One-colour theorem used twice

`general-period-homological-free-energy-20260914.md` gives for either graph `G`

\[
\frac1{\rho_{G,n}}
\log P(r_G>0)
\to-\max\{1-\alpha_G,0\}                                    \tag{2.1}
\]

when `alpha_G` has a limit.

Therefore:

- if `alpha_4<1`, then black positive rank is exponentially rare and
  \[
  P(r_4=0)\to1;                                               \tag{2.2}
  \]
- if `alpha_4>1`, then
  \[
  P(r_4>0)\to1;                                               \tag{2.3}
  \]
- if `alpha_8<1`, then white matching positive rank is exponentially rare;
- if `alpha_8>1`, then white matching positive rank occurs with probability tending to one.

## 3. Digital Alexander converts the white statement into the black upper-rank statement

Configuration by configuration under complementary labels,

\[
\boxed{r_4(\omega)+r_8(\omega^c)=2.}                         \tag{3.1}
\]

Hence in distribution at black density `p` / white density `q`,

\[
P_p(r_4=2)=P_q(r_8=0),                                       \tag{3.2}
\]

\[
P_p(r_4=0)=P_q(r_8=2).                                       \tag{3.3}
\]

Thus the matching free energy controls the **second** NN rank birth.

## 4. Three open phases away from the free-energy boundaries

Suppose along a sequence

\[
\alpha_{4,n}(p)\to\alpha_4,
\qquad
\alpha_{8,n}(q)\to\alpha_8,                                 \tag{4.1}
\]

with neither limit equal to one.

### Phase 0: no black homology

If

\[
\boxed{\alpha_4<1,}                                           \tag{4.2}
\]

then

\[
\boxed{r_4\xrightarrow{P}0.}                                 \tag{4.3}
\]

Consequently Alexander duality forces

\[
r_8(\omega^c)\xrightarrow{P}2.                              \tag{4.4}
\]

### Phase 2: full black homology

If

\[
\boxed{\alpha_8<1,}                                           \tag{4.5}
\]

then white matching rank is zero with probability tending to one, so

\[
\boxed{r_4\xrightarrow{P}2.}                                 \tag{4.6}
\]

### Phase 1: neutral rank-one plateau

If

\[
\boxed{\alpha_4>1\quad\text{and}\quad\alpha_8>1,}            \tag{4.7}
\]

then

\[
P(r_4=0)\to0                                                  \tag{4.8}
\]

from the NN theorem, while

\[
P(r_4=2)=P(r_8=0)\to0                                        \tag{4.9}
\]

from the matching theorem.  Hence

\[
\boxed{r_4\xrightarrow{P}1.}                                 \tag{4.10}
\]

So the macroscopic rank-one phase is exactly the region where **both colours have enough opportunities to create positive homology**.

## 5. The apparently fourth phase is forbidden

Could both

\[
\alpha_4<1
\quad\text{and}\quad
\alpha_8<1                                                     \tag{5.1}
\]

hold along an asymptotic sequence away from the boundary?

The one-colour free-energy theorem would then give simultaneously

\[
P(r_4=0)\to1,
\qquad
P(r_8=0)\to1.                                                 \tag{5.2}
\]

But from (3.1), black rank zero forces complementary white rank two on the same configuration.  In particular the exact marginal inequality

\[
P_p(r_4=0)+P_q(r_8=0)
=P_p(r_4=0)+P_p(r_4=2)
\le1                                                         \tag{5.3}
\]

holds at every finite size.

Therefore the double-rank-zero phase is impossible:

\[
\boxed{(\alpha_4,\alpha_8)\in(-\infty,1)^2
\text{ cannot be an off-boundary asymptotic limit.}}          \tag{5.4}
\]

This is a nontrivial compatibility condition between the two complementary correlation-norm minima and the period-lattice entropy.

## 6. Matching observable in the three phases

Recall

\[
M(p)=P(r_4=2)-P(r_4=0).                                       \tag{6.1}
\]

Therefore

\[
\boxed{
M\to
\begin{cases}
-1,&\alpha_4<1,\\
0,&\alpha_4>1,\ \alpha_8>1,\\
+1,&\alpha_8<1.
\end{cases}}                                                  \tag{6.2}
\]

The balance root in a moving-parameter problem is the finite-size location where the two charged endpoint sectors exchange dominance through the rank-one region.

This phase diagram explains why the balance observable may have a sharp zero even when a wide region has `M≈0`: the bulk neutral phase is controlled by two negative free energies, while the actual root is an `O(1)` topological-charge balance inside it.

## 7. Free-energy sign formulation

In terms of (1.5):

- `Phi_4 >> +1` -> rank 0;
- `Phi_8 >> +1` -> rank 2;
- `Phi_4 << -1` and `Phi_8 << -1` -> rank 1.

The exact topology forbids `Phi_4` and `Phi_8` from both being macroscopically positive in a regime where the fixed-`p` theorem applies.

The codimension-one surfaces

\[
\Phi_4\approx0,
\qquad
\Phi_8\approx0                                                \tag{7.1}
\]

are the lower and upper homology birth boundaries.  Their separation or merger is a geometry/parameter question; the state labels between them are already fixed by topology.

## 8. Exponential-aspect special case

For a short direction `e` and transverse entropy `d`,

\[
\rho_4\sim\ell\tau_{4,p}(e),
\qquad
\rho_8\sim\ell\tau_{8,1-p}(e),
\qquad
\log N\sim d\ell.                                             \tag{8.1}
\]

Thus the phase boundaries are

\[
\tau_{4,p}(e)=d,                                              \tag{8.2}
\]

and

\[
\tau_{8,1-p}(e)=d,                                            \tag{8.3}
\]

which are exactly the directional birth-centre equations already proved on this branch.

The general phase diagram therefore contains the axial and varying-direction exponential tori as one special geometry, rather than being derived from them.

## 9. Boundary regimes

If one or both `alpha` values tend to one, the normalized free-energy theorem only fixes exponential rates.  The actual endpoint probabilities can have nontrivial limits determined by subexponential component activities, prefactors and correlations.

That is precisely the regime where the Poisson-window and common-crossover analyses become necessary.  The phase diagram should not be used to replace those finer results at `alpha=1`.

## 10. Claim boundary

The three-phase classification uses only the fixed-`p` arbitrary-shape free-energy theorem and exact finite digital Alexander duality.  It does not assume independence of black and white events and does not posit a continuum field.