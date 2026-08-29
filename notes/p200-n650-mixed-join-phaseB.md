# N650 Phase B: a scoreable mixed-factor connectivity defect

The exact CRT result rules out an antisymmetric linear path commutator.  This
note replaces it with a quantity that is nonlinear, canonical on every final
configuration, and explicit about what a nonzero result means.

## 1. Work on typed final-lift connectivity, not a binary quotient

A generic `N650` binary configuration does not descend to either intermediate
torus.  Therefore no OR, majority, or representative rule is silently used.
For color (c\in\{B,W\}), let (\Pi_c(\omega)) be the complete connectivity
partition of the final configuration:

- black uses the square NN graph;
- white uses the NN+NNN matching graph;
- the implementation retains universal-cover displacement potentials even
  though the primary partition-rank statistic does not need them.

Within each of the 65 source fibers, label the ten lifts by
(D=C_2\times C_5).  Let (R_2^c(\omega)) relate same-color lifts that differ
only in the (C_2) coordinate, and let (R_5^c(\omega)) relate same-color
lifts that differ only in the (C_5) coordinate.  These are artificial join
edges with their honest deck displacements, not a new percolation lattice.

For a partition (P) of the colored vertices, put

\[
 r(P)=|V_c|-\#\operatorname{blocks}(P).
\]

The full mixed redundancy is

\[
 J_c^{\rm full}=
 r(\Pi_c\vee R_2^c)+r(\Pi_c\vee R_5^c)
 -r(\Pi_c\vee R_2^c\vee R_5^c)-r(\Pi_c).
\]

Partition-rank submodularity proves (J_c^{\rm full}\ge0).  This is not an
order commutator: both join orders still give the same final partition.

## 2. Remove the exact isolated-fiber interaction

Even ten disconnected Bernoulli bits have a local mixed join.  For selected
cells (S\subset C_2\times C_5), form the bipartite incidence graph whose row
vertices are occupied (C_2) classes, whose column vertices are occupied
(C_5) classes, and whose edges are the selected cells.  Direct rank counting
gives

\[
 J_c^{\rm iso}(S)=b_1(\operatorname{Inc}(S)).
\]

This definition is invariant under separate relabeling of both CRT factors.
The closed form looks asymmetric only because (C_2) has two rows: if
(k_c) is the number of columns containing both row lifts, then

\[
 J_c^{\rm iso}=\max(k_c-1,0).
\]

The production observable is the connected residual

\[
 R_c(\omega)=J_c^{\rm full}(\omega)
 -\sum_{x=1}^{65}J_{c,x}^{\rm iso}(\omega).
\]

It vanishes configurationwise when the source fibers have no inter-fiber
connectivity.  A nonzero mean therefore identifies a **nonlocal mixed
(C_2\)-by-(C_5) connectivity interaction**.  It does not, by itself,
identify Jordan structure, RG memory, or path noncommutativity.

## 3. Exact null and the boundary of the radial clocks

For any four-corner response with factor-additive form

\[
 h_{ab}=h_{00}+a u_2+b u_5,
\]

the mixed difference is exactly

\[
 h_{11}-h_{10}-h_{01}+h_{00}=0.
\]

Adding a bilinear term (\lambda ab) makes the same difference exactly
(\lambda).  Thus the score is a direct gate for a mixed factor sector.

The ordinary `q2` and Jordan clocks constrain unmarked endpoint states; alone
they say nothing about this newly marked topology functional.  Their
**minimal factor-additive topology extension** predicts (E[R_c]=0).  A
rejection falsifies that bridge, not automatically the radial clock selected
by N580.  This distinction prevents a finite-lattice connectivity effect from
being renamed “morphism memory” without a scaling test.

## 4. Frozen four-vector score

For each orientation (o\in\{(23,11),(17,19)\}), define

\[
 E_o=R_B^o+R_W^o,\qquad O_o=R_B^o-R_W^o.
\]

The primary state is

\[
 y=(E_S,E_D,O_S,O_D)
\]

with (S=(o_1+o_2)/2) and (D=(o_1-o_2)/2).  The factor-additive null is the
fixed zero vector.  Use one synchronized counter stream for both orientations
and compute all four joins on the same configuration.  At the frozen

\[
 p_{\rm ref}=0.592746050790,
\]

draw (K\sim\operatorname{Binomial}(650,p_{\rm ref})) for each permutation
and evaluate its prefix; use the same (K) for both orientations.  At least
100 synchronized batches give the full delete-one (4\times4) covariance.
The primary score is the joint GLS statistic (y^T\Sigma^+y), with degrees
of freedom equal to the numerical covariance rank.  Marginal signed z values
classify a rejection:

- `ES`: scalar/color-even interaction;
- `ED`: orientation-sensitive color-even interaction;
- `OS`: matching-odd orientation average;
- `OD`: chiral matching-odd interaction, the only channel worth comparing
  prospectively with the archived P57 odd direction.

Do not fit the P57 template in the primary score.  The new statistic is a
scalar-at-(p_{\rm ref}) topology gate, not yet an `r=2..6` curve jet.

## 5. Tiny exact normalization

All (2^{10}=1024) fiber colorings are exhausted by the oracle.  At (p=1/2),

\[
 E[J_B]=E[J_W]=\frac{499}{1024},\qquad
 \operatorname{Var}(J_B-J_W)=\frac{681}{512}.
\]

The odd-count distribution from `-4` through `4` is

```text
1, 15, 80, 210, 412, 210, 80, 15, 1.
```

The artifact stores decimal moments obtained by exact rational enumeration at
the frozen (p_{\rm ref}); the rational input and executable derivation are
retained without serializing enormous power-of-ten fractions.  For a readable
dimensionless diagnostic, divide one
orientation's odd connected residual by

\[
 \sqrt{65\,\operatorname{Var}_{p_{\rm ref}}(J_B-J_W)}.
\]

This is only an exact reference scale; the empirical synchronized jackknife
covariance remains primary.

The included two-fiber witness has zero isolated interaction but one mixed
cycle after ordinary connectivity joins the two fibers, so (R=1).  It proves
the subtraction does not erase the desired nonlocal mechanism.

## Secondary lift-aware archive

For each of the four corners, also archive

\[
 \rho_c^A=\operatorname{rank}\operatorname{im}
 [H_1(G_c^A)\to H_1(T^2)].
\]

Its mixed four-corner difference is a typed secondary observable.  It is not
folded into the primary score until a tiny real-HNF lift oracle freezes its
sign and displacement conventions.

## Reproduce

```bash
python3 scripts/gaussian_crt_mixed_join_phaseb.py \
  --output predictions/p200_n650_mixed_join_phaseB_20260829.json
python3 -m unittest discover -s tests -p 'test_gaussian_crt_mixed_join_phaseb.py'
```
