# Near-critical loop versus endpoint insertion: a decisive prefactor diagnostic

2026-09-14.  Conjectural synthesis of three rigorous ingredients that currently live in different models/levels:

1. the finite-state cyclic matrix-sewing theorem on this branch, whose simple Perron-band log-determinant residue is exactly one;
2. the 2026 D'Alimonte--Manolescu near-critical OZ theorem for square-lattice **bond** random-cluster connectivity;
3. their uniform Brownian/local-CLT scale, which gives transverse variance of order `w xi(p)`.

The goal is not to claim a SITE theorem.  It is to isolate a sharp observable that distinguishes whether a complete winding component behaves like a genuinely cyclic renewal object or still carries point-to-point endpoint insertions near criticality.

## 1. Correlation-length variables

Let

\[
\xi=\kappa^{-1},
\qquad
s=w/\xi=\kappa w.                                             \tag{1.1}
\]

Here `w` is the physical winding length and `s` is the number of correlation-length units around the loop.

The near-critical OZ theorem for the bond-FK two-point function has the form

\[
\boxed{
G_{2pt}(w)
\asymp
\pi_1(\xi)^2s^{-1/2}e^{-s}.}                                 \tag{1.2}
\]

The conditioned-cluster transverse variance is of order

\[
\operatorname{Var}X_\perp\asymp w\xi.                        \tag{1.3}
\]

Thus a Brownian diffusion coefficient per unit physical longitudinal distance satisfies

\[
D\asymp\xi.                                                    \tag{1.4}
\]

## 2. What a pure cyclic renewal loop predicts

For a cyclic Markov-additive renewal kernel with one simple Perron band, the branch theorem `matrix-sewing-unit-residue-20260914.md` gives

\[
\nu_w
\sim
\frac{\zeta_{loop}}{\sqrt{2\pi D w}}e^{-\kappa w}.            \tag{2.1}
\]

If the cyclic closure has unit insertion residue at correlation-length scale,

\[
\zeta_{loop}=O(1)                                             \tag{2.2}
\]

as `p->pc`, then using `D~xi` and `w=xi s`,

\[
\boxed{
\nu_w^{cyclic}
\asymp
\xi^{-1}s^{-1/2}e^{-s}.}                                     \tag{2.3}
\]

The factor `xi^{-1}` has a simple local-CLT meaning.  A transverse bridge after `s` correlation-length steps has physical standard deviation

\[
\xi\sqrt s.                                                   \tag{2.4}
\]

Returning to one specified **microscopic row** therefore costs

\[
(\xi\sqrt s)^{-1}.                                            \tag{2.5}
\]

This is exactly the physical-row version of the Gaussian closure factor.

## 3. Why the two-point function has a different insertion weight

The two-point probability (1.2) is dimensionless and pins two microscopic endpoint vertices.  Compare it with the raw microscopic endpoint local-CLT density `(xi sqrt s)^{-1}e^{-s}` suggested by the renewal walk.

The ratio is

\[
\boxed{
J_{2pt}(\xi)
\asymp
\pi_1(\xi)^2\xi.}                                            \tag{3.1}
\]

So in the near-critical bond-FK theorem, the two endpoint insertions together contribute an effective factor of order

\[
\pi_1(\xi)^2\xi.                                              \tag{3.2}
\]

Heuristically one may think of each microscopic endpoint insertion as having scale

\[
\pi_1(\xi)\sqrt\xi,                                          \tag{3.3}
\]

though only the product is relevant here.

A genuine closed loop has no prescribed endpoints.  Therefore there is no reason for (3.2) to survive unchanged after cyclic closure and unrooting.

## 4. Three competing complete-component insertion hypotheses

Write the general near-critical complete-component ansatz as

\[
\boxed{
\nu_w(p)
\asymp
\xi^{-1}J_{comp}(\xi,s)
 s^{-1/2}e^{-s}.}                                             \tag{4.1}
\]

The unknown is now isolated in one dimensionless insertion factor `J_comp`.

### H0: pure cyclic closure

\[
\boxed{J_{comp}\asymp1.}                                     \tag{4.2}
\]

Then

\[
\nu_w\asymp\xi^{-1}s^{-1/2}e^{-s}.                           \tag{4.3}
\]

This is the natural continuation of the unit-residue logdet mechanism.

### H2: two endpoint insertions survive

If cutting/opening the complete component effectively introduces the same microscopic endpoint factors as a two-point connection, then

\[
J_{comp}\asymp\pi_1(\xi)^2\xi,                               \tag{4.4}
\]

and

\[
\boxed{
\nu_w\asymp\pi_1(\xi)^2s^{-1/2}e^{-s},}                     \tag{4.5}
\]

at comparability level.

### H1 / marked closure

A one-mark or asymmetric anchor construction could produce an intermediate factor.  More generally write

\[
J_{comp}(\xi)\asymp\xi^\beta\pi_1(\xi)^\gamma               \tag{4.6}
\]

as a diagnostic parameterization only, not a claimed power law.

The important point is that `beta,gamma` describe **insertion semantics**, while the universal Gaussian closure remains `s^{-1/2}`.

## 5. The three hypotheses are invisible at fixed p

At any fixed subcritical `p`,

\[
\xi(p)<\infty,
\qquad
\pi_1(\xi(p))>0                                               \tag{5.1}
\]

are constants.  Every hypothesis therefore reduces to

\[
\nu_w=C(p)w^{-1/2}e^{-\kappa(p)w}(1+o(1))                    \tag{5.2}
\]

with a different `C(p)`.

So a fixed-p width ladder, even a mathematically perfect one, cannot determine whether the prefactor is a pure cyclic residue or an endpoint-dressed insertion mechanism.

The hypotheses separate only when `p=p_w->pc` and `xi(p_w)->infinity`.

## 6. Distinct extreme-value centre corrections

Suppose a longitudinal opportunity count `m` creates complete components with mean

\[
\lambda=m\nu_w.                                               \tag{6.1}
\]

At first birth, `lambda=O(1)`.  From (4.1), the general logarithmic balance is

\[
\boxed{
s+\tfrac12\log s+\log\xi-\log J_{comp}(\xi,s)
=\log m+O(1).}                                                \tag{6.2}
\]

### Pure cyclic H0

\[
\boxed{
s+\tfrac12\log s+\log\xi
=\log m+O(1).}                                                \tag{6.3}
\]

### Endpoint-dressed H2

Insert (4.4):

\[
\boxed{
s+\tfrac12\log s-2\log\pi_1(\xi)
=\log m+O(1),}                                                \tag{6.4}
\]

which is exactly the two-point OZ balance.

The difference between (6.3) and (6.4) is

\[
\boxed{
\log\xi+2\log\pi_1(\xi).}                                   \tag{6.5}
\]

This grows logarithmically/polynomial-arm scale rather than staying `O(1)`.  Near criticality it is therefore a strong discriminator.

## 7. A direct dimensionless ratio diagnostic

Define

\[
\boxed{
R_{loop/2pt}(p,w)
=\frac{\xi(p)\nu_w(p)}{G_{2pt}(w;p)}.}                        \tag{7.1}
\]

Using the common Gaussian/exponential factors:

### H0 predicts

\[
\boxed{
R_{loop/2pt}\asymp\pi_1(\xi)^{-2}.}                          \tag{7.2}
\]

### H2 predicts

\[
\boxed{R_{loop/2pt}\asymp\xi.}                               \tag{7.3}
\]

depending on the precise normalization chosen for the two-point denominator.  An even cleaner comparison is to strip the common skeleton directly:

\[
\boxed{
J_{emp}(p,w)
:=\xi\sqrt{s}\,e^s\nu_w(p).}                                 \tag{7.4}
\]

Then

\[
J_{emp}\asymp1                                                \tag{7.5}
\]

under pure cyclic closure, while

\[
J_{emp}\asymp\pi_1(\xi)^2\xi                                \tag{7.6}
\]

under two-endpoint dressing.

This is the preferred statistic because it does not require a separate two-point simulation if `xi` and the critical one-arm baseline are available.

## 8. Interaction with the surface-excess identity

The branch already gives an independent component-Palm route to the mass slope through

\[
qE_WB-pE_WN=pq\,\partial_p\log\nu_w.                          \tag{8.1}
\]

If (4.1) holds, then

\[
\partial_p\log\nu_w
= -\partial_p\log\xi
  +\partial_p\log J_{comp}
  -\frac12\partial_p\log s
  -\partial_p s.                                              \tag{8.2}
\]

Near a crossover centre, the dominant term may still be `-partial_p s`, but the insertion derivative can enter at logarithmic order.

Thus morphology `(N,B)` gives a second route to detecting whether the near-critical prefactor is becoming singular, without fitting the activity alone.

## 9. Relation to matrix sewing

The finite-state matrix theorem says:

- a simple Perron band in a pure cyclic logdet contributes unit logarithmic residue;
- finite local memory changes `D` but not the residue;
- a nontrivial continuously varying residue must arise from insertion/mark/unrooting semantics, multiple bands or an infinite-state limit.

Near-critical divergence of the correlation scale naturally sends a fixed microscopic transfer description towards an effectively infinite-state object.  Therefore observing `J_comp` drift does not refute Gaussian sewing.  It diagnoses which microscopic insertion survives the scaling limit.

This is exactly the distinction that fixed-p data could not make.

## 10. A concrete numerical/theoretical gate

Choose a sequence `p_j->pc` and widths

\[
w_j=s_j\xi(p_j),
\qquad
s_j\to\infty\text{ slowly},                                  \tag{10.1}
\]

so the system is safely on the OZ side but the correlation length grows.

For each `j`, estimate/certify

\[
J_{emp,j}=\xi_j\sqrt{s_j}e^{s_j}\nu_{w_j}(p_j).              \tag{10.2}
\]

Interpretation:

- bounded nonzero `J_emp` -> supports pure cyclic closure;
- tracking `pi_1(xi)^2 xi` -> supports endpoint-dressed closure;
- another systematic scale -> identifies a different insertion class;
- extra power of `s` -> refutes the simple one-soft-band Gaussian sewing hypothesis itself.

The existing tagged complete-component resolvent is exact at finite width but cannot reach `w~xi->infinity` by brute force.  A correlation-length-block transfer/renewal construction is therefore the right next theorem engine, not simply a larger microscopic width ladder.

## 11. Claim boundary

All SITE formulas in this note are conjectural diagnostics.  The bond-FK two-point formula and Brownian scale are rigorous in D'Alimonte--Manolescu; the unit cyclic residue is rigorous for the finite-state matrix class on this branch.  The scientific question is how the actual square-site complete-component observable interpolates between those structures.