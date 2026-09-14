# One-sided free energies and the three homology-rank phases

2026-09-14.  Corrected formulation after self-audit.  The NN and complementary matching inverse-correlation norms are **not simultaneously subcritical at one common p**, because

\[
p_c(G8)=1-p_c(G4).                                            \tag{0.1}
\]

Therefore the lower and upper rank boundaries are controlled by subcritical free energies on **opposite sides of the critical point**, not by two same-p subcritical masses.

The correct phase diagram is:

- below `p_c`: the black NN free energy separates rank 0 from rank 1;
- above `p_c`: the complementary white matching free energy separates rank 1 from rank 2;
- at/near the boundaries: component-Poisson/crossover information is needed.

## 1. Lower-side NN free energy, p < pc(G4)

Fix

\[
0<p<p_c(G4).                                                   \tag{1.1}
\]

For an honest period lattice `Lambda_n` with `N_n` sites define

\[
\rho_{4,n}(p)
=\min_{\lambda\in\Lambda_n\setminus0}\tau_{4,p}(\lambda).   \tag{1.2}
\]

Assume

\[
\rho_{4,n}(p)\to\infty,
\qquad
\alpha_{4,n}(p)=\frac{\log N_n}{\rho_{4,n}(p)}\to\alpha_4.   \tag{1.3}
\]

The arbitrary-shape free-energy theorem gives

\[
\frac1{\rho_{4,n}}
\log P_p(r_4>0)
\to-\max\{1-\alpha_4,0\}.                                   \tag{1.4}
\]

### 1.1 If alpha4 < 1: rank zero

Then positive black homology is exponentially rare:

\[
\boxed{r_4\xrightarrow P0.}                                  \tag{1.5}
\]

### 1.2 If alpha4 > 1: positive rank appears

Then

\[
P_p(r_4>0)\to1.                                               \tag{1.6}
\]

To identify the rank as one rather than two, use the arbitrary-period subcritical rank comparison already proved in the geometric manuscript: for every fixed `p<pc(G4)` and sufficiently large Euclidean shortest period `ell_n`,

\[
\boxed{
\frac{P_p(r_4=2)}{P_p(r_4=0)}
\le
\exp[-c_p N_n/\ell_n].}                                      \tag{1.7}
\]

Because `rho_4->infinity`, norm equivalence implies `ell_n->infinity`.  Equation (1.6) gives `P_0->0`; multiplying by the exponentially small ratio in (1.7) yields

\[
P_2\to0.                                                       \tag{1.8}
\]

Hence

\[
\boxed{r_4\xrightarrow P1\qquad(p<p_c,\ \alpha_4>1).}        \tag{1.9}
\]

Thus below criticality the single black NN free energy controls the transition

\[
\boxed{0\longleftrightarrow1.}                               \tag{1.10}
\]

## 2. Upper-side matching free energy, p > pc(G4)

Now fix

\[
p>p_c(G4),
\qquad
q=1-p<p_c(G8).                                                \tag{2.1}
\]

Define the **subcritical white matching** period cost

\[
\rho_{8,n}(q)
=\min_{\lambda\in\Lambda_n\setminus0}\tau_{8,q}(\lambda),   \tag{2.2}
\]

and suppose

\[
\rho_{8,n}(q)\to\infty,
\qquad
\alpha_{8,n}(q)=\frac{\log N_n}{\rho_{8,n}(q)}\to\alpha_8.   \tag{2.3}
\]

Digital Alexander duality is exact:

\[
\boxed{r_4(\omega)+r_8(\omega^c)=2.}                         \tag{2.4}
\]

### 2.1 If alpha8 < 1: white rank zero, black rank two

The matching free-energy theorem gives

\[
r_8(\omega^c)\xrightarrow P0.                               \tag{2.5}
\]

Therefore

\[
\boxed{r_4\xrightarrow P2.}                                  \tag{2.6}
\]

### 2.2 If alpha8 > 1: white positive rank is rank one

The matching theorem gives

\[
P_q(r_8>0)\to1.                                               \tag{2.7}
\]

Apply the same subcritical rank comparison (1.7), now to graph `G8` at `q<pc(G8)`.  Its rank-two probability is negligible relative to its rank-zero probability.  Since rank zero itself tends to zero under `alpha8>1`,

\[
P_q(r_8=2)\to0.                                               \tag{2.8}
\]

Hence

\[
r_8(\omega^c)\xrightarrow P1.                               \tag{2.9}
\]

and Alexander duality gives

\[
\boxed{r_4\xrightarrow P1\qquad(p>p_c,\ \alpha_8>1).}        \tag{2.10}
\]

Thus above criticality the single **white matching** free energy controls

\[
\boxed{1\longleftrightarrow2.}                               \tag{2.11}
\]

## 3. Correct three-phase picture

Away from the free-energy boundaries and away from the critical point:

\[
\boxed{
\begin{array}{c|c|c}
\text{parameter side}&\text{free-energy condition}&\text{black NN rank}\\
\hline
p<p_c&\log N/\rho_4(p)<1&0\\
p<p_c&\log N/\rho_4(p)>1&1\\
p>p_c&\log N/\rho_8(1-p)>1&1\\
p>p_c&\log N/\rho_8(1-p)<1&2
\end{array}}                                                  \tag{3.1}
\]

The rank-one phase is the geometry-induced bridge between the lower NN birth and the upper complementary-matching birth.

There is no regime in which one should simultaneously plug `p` into a subcritical NN mass and `1-p` into a subcritical matching mass: by (0.1), exactly one of those two complementary parameters is subcritical away from `p_c`.

## 4. Matching observable

Recall

\[
M(p)=P_2(p)-P_0(p).                                           \tag{4.1}
\]

The phase limits are

\[
\boxed{
M\to
\begin{cases}
-1,&p<p_c,\ \log N/\rho_4(p)<1,\\
0,&p<p_c,\ \log N/\rho_4(p)>1,\\
0,&p>p_c,\ \log N/\rho_8(1-p)>1,\\
+1,&p>p_c,\ \log N/\rho_8(1-p)<1.
\end{cases}}                                                  \tag{4.2}
\]

This makes the plateau mechanism explicit.  The lower free-energy boundary drives `M` from `-1` to `0`; the upper complementary boundary drives it from `0` to `+1`.  The finite matching root lies inside the central topological crossover and is not the same object as either fixed-p free-energy threshold.

## 5. Exponential-aspect special case

Let the shortest-period direction converge to `e` and

\[
\frac{\log(N/\ell)}{\ell}\to d>0.                            \tag{5.1}
\]

Below `p_c`,

\[
\rho_4(p)=\ell\tau_{4,p}(e)+o(\ell),
\qquad
\log N=d\ell+o(\ell).                                        \tag{5.2}
\]

The lower boundary is therefore

\[
\boxed{\tau_{4,p}(e)=d,}                                     \tag{5.3}
\]

which gives `p=a_e(d)`.

Above `p_c`, write `q=1-p<pc(G8)`.  Then

\[
\rho_8(q)=\ell\tau_{8,q}(e)+o(\ell),                         \tag{5.4}
\]

and the upper boundary is

\[
\boxed{\tau_{8,1-p}(e)=d,}                                   \tag{5.5}
\]

which gives `p=b_e(d)`.

So the directional two-centre theorem is exactly the one-dimensional section of this corrected lower/upper free-energy picture.

## 6. Boundary alpha = 1

At

\[
\log N/\rho\to1,                                              \tag{6.1}
\]

the arbitrary-shape theorem only gives

\[
\log P(r>0)=o(\rho).                                         \tag{6.2}
\]

The actual endpoint-sector probability can be nontrivial and depends on subexponential component activities.  This is precisely where the Poisson birth-window theory enters.

Thus the hierarchy is:

1. correlation-norm free energy locates the rank boundary at exponential scale;
2. component intensity resolves the boundary law;
3. its derivative/semiconvexity resolves the affine Gumbel window at regular points.

## 7. Critical p itself

At `p=p_c(G4)`, both complementary graphs are critical and neither subcritical correlation norm supplies a positive free energy.  This note deliberately makes no rank-law claim at the critical parameter for arbitrary aspect sequences.

Critical fixed-aspect wrapping laws and the geometric balance-root theorem are separate inputs/problems.

## 8. Claim boundary

This corrected phase diagram uses the arbitrary-shape subcritical free-energy theorem, the geometric manuscript's subcritical `P2/P0` suppression, and exact digital Alexander duality.  It does not analytically continue an inverse-correlation norm through criticality.