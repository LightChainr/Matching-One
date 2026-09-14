# Fixed-width charge-sector free energy and the infinite-length balance root

2026-09-14.  A one-dimensional transfer/free-energy explanation of “balance without concentration.”  The order of limits in the main theorem is different, but the fixed-width `m->infinity` problem cleanly identifies what the matching root is balancing.

## 1. Trivial-homology survival in an open strip

Fix a circumference `w` and graph `G=G4` or `G8`.  On the open vertical strip

\[
C_w\times\{1,\ldots,m\},                                     \tag{1.1}
\]

let

\[
Q^G_{w,m}(p)
=P_p(\text{there is no horizontally essential occupied component}).\tag{1.2}
\]

There is no vertical periodic identification in this strip; the only ambient homology to forbid is horizontal winding around `C_w`.

### 1.1 Submultiplicative upper inequality

If a strip of height `m+n` has no horizontal essential component, then neither its first `m` rows nor its last `n` rows has one.  These two subevents use disjoint row variables.  Therefore

\[
\boxed{
Q_{w,m+n}^G(p)
\le Q_{w,m}^G(p)Q_{w,n}^G(p).}                               \tag{1.3}
\]

Thus

\[
a_m=-\log Q_{w,m}^G(p)                                      \tag{1.4}
\]

is superadditive and Fekete gives the limit

\[
\boxed{
I^0_{G,w}(p)
:=\lim_{m\to\infty}-\frac1m\log Q_{w,m}^G(p)
=\sup_m-\frac1m\log Q_{w,m}^G(p).}                           \tag{1.5}
\]

### 1.2 Reverse inequality up to one empty separator row

Let

\[
\delta_{w,p}=(1-p)^w>0.                                      \tag{1.6}
\]

Force one whole row between two strips to be empty.  For both NN and matching adjacency, one empty row separates occupied components above and below because all vertical jumps have size one.

Hence

\[
\boxed{
Q_{w,m+n+1}^G(p)
\ge\delta_{w,p}Q_{w,m}^G(p)Q_{w,n}^G(p).}                    \tag{1.7}
\]

So the strip survival probability has a genuine one-dimensional free energy with only `O(1)` concatenation cost.

## 2. The torus rank-zero probability has the same exponential rate

Let

\[
P^G_{0;w,m}(p)
=P_p^{C_w\times C_m}(r_G=0).                                 \tag{2.1}
\]

### Upper bound

A rank-zero torus configuration has no horizontal essential component.  Cut the torus between two rows and forget the vertical periodic edges.  The resulting open strip still has no horizontal winding.  Therefore

\[
\boxed{P^G_{0;w,m}(p)\le Q^G_{w,m}(p).}                      \tag{2.2}
\]

### Lower bound

Force one specified torus row to be completely empty.  Cut there.  If the remaining open strip has no horizontal essential component, the torus has no horizontal homology, and the empty row also destroys every possible vertical homology cycle.

Thus

\[
\boxed{
P^G_{0;w,m}(p)
\ge\delta_{w,p}Q^G_{w,m-1}(p).}                              \tag{2.3}
\]

Combining (1.5), (2.2), and (2.3),

\[
\boxed{
\lim_{m\to\infty}
-\frac1m\log P^G_{0;w,m}(p)
=I^0_{G,w}(p).}                                               \tag{2.4}
\]

This statement is exact for every fixed `w` and interior `p`.

## 3. Finite-state Perron representation

For fixed `w`, horizontal-winding avoidance is recognized by the standard finite frontier connectivity/homology state used throughout the repository.

Reject a transition as soon as horizontal homology is created.  This gives a finite nonnegative substochastic transfer matrix

\[
T^0_{G,w}(p).                                                  \tag{3.1}
\]

Its entries are polynomials in `p` and `1-p`.

The all-empty row gives a reset state accessible with positive probability from every safe frontier state.  Restricting to reachable/co-reachable safe states therefore produces a primitive Perron block.  Hence

\[
\boxed{
I^0_{G,w}(p)=-\log\lambda^0_{G,w}(p),}                        \tag{3.2}
\]

where `lambda^0` is the Perron root of that block.

For `0<p<1`, the Perron root is real analytic in `p`.  Monotone coupling makes it nonincreasing; strict finite-energy enhancement of horizontal winding makes

\[
\boxed{I^0_{G,w}(p)\text{ strictly increasing in }p.}         \tag{3.3}
\]

A fully machine-checked version can be obtained by constructing exactly the same safe-state automaton already used for tagged/cylinder calculations; no large-`w` claim is required for the fixed-width theorem.

## 4. Charge fugacity has a fixed-width thermodynamic rate

For black NN at density `p`, digital Alexander gives

\[
P^{4}_{2;w,m}(p)
=P^{8}_{0;w,m}(1-p).                                          \tag{4.1}
\]

The charge fugacity is

\[
\theta_{w,m}(p)
=\log\frac{P^4_{2;w,m}(p)}{P^4_{0;w,m}(p)}.                  \tag{4.2}
\]

Using (2.4) for the two graphs,

\[
\boxed{
\lim_{m\to\infty}\frac1m\theta_{w,m}(p)
=\Theta_w(p),}                                                \tag{4.3}
\]

with

\[
\boxed{
\Theta_w(p)
=I^0_{4,w}(p)-I^0_{8,w}(1-p).}                               \tag{4.4}
\]

Thus the large-`m` matching balance is a **difference of two topological void free energies**, not a difference of the black/white essential-component intensities.

That distinction matters because the infinite-cylinder black and complementary-white essential-component intensities are equal by alternation, while their long-interval void probabilities need not have the same exponential rate.

## 5. Unique fixed-width charge-coexistence point

`I^0_{4,w}(p)` is strictly increasing in `p`, while

\[
I^0_{8,w}(1-p)                                                \tag{5.1}
\]

is strictly decreasing in `p`.  Therefore

\[
\boxed{\Theta_w(p)\text{ is strictly increasing}.}           \tag{5.2}
\]

At `p downarrow0`, black horizontal winding is extremely unlikely, so

\[
I^0_{4,w}(p)\to0,                                             \tag{5.3}
\]

while at matching density `1-p up to1`, rank-zero survival is strongly suppressed, giving

\[
\Theta_w(p)<0                                                 \tag{5.4}
\]

for sufficiently small `p`.

The opposite sign holds near `p=1`.  Hence there is a unique

\[
\boxed{p^{ch}_w\in(0,1):\Theta_w(p^{ch}_w)=0.}                \tag{5.5}
\]

Equivalently,

\[
\boxed{
I^0_{4,w}(p^{ch}_w)
=I^0_{8,w}(1-p^{ch}_w).}                                     \tag{5.6}
\]

This is the infinite-length fixed-width **charge coexistence point**.

The finite reflection-dominance theorem implies

\[
\boxed{p^{ch}_w\ge1/2,}                                      \tag{5.7}
\]

with strict inequality for ordinary square cylinders once strict matching enhancement is reflected in the transfer sectors.

## 6. Finite-torus roots converge to the charge-coexistence point as m -> infinity

Let

\[
p^*_{w,m}                                                     \tag{6.1}
\]

be the unique finite-torus root where

\[
P_2=P_0.                                                       \tag{6.2}
\]

Fix any `epsilon>0`.  By strict monotonicity of `Theta_w`,

\[
\Theta_w(p^{ch}_w-\epsilon)<0,
\qquad
\Theta_w(p^{ch}_w+\epsilon)>0.                                \tag{6.3}
\]

Equation (4.3) then gives, for sufficiently large `m`,

\[
\theta_{w,m}(p^{ch}_w-\epsilon)<0,
\qquad
\theta_{w,m}(p^{ch}_w+\epsilon)>0.                            \tag{6.4}
\]

Since the finite fugacity is strictly increasing and its root is `p^*_{w,m}`,

\[
\boxed{
p^*_{w,m}\longrightarrow p^{ch}_w
\qquad(m\to\infty).}                                        \tag{6.5}
\]

If one additionally extracts the Perron amplitudes and `Theta_w'(p^{ch}_w)>0`, standard analytic transfer asymptotics predict an `O(1/m)` root correction.  That quantitative rate is not needed for (6.5).

## 7. Then the charge-coexistence point tends to the planar pc

The geometric-balance manuscript proves balance-root consistency for any honest torus sequence whose Euclidean shortest period tends to infinity.

Choose widths

\[
w_j\to\infty.                                                 \tag{7.1}
\]

For each `w_j`, choose `m_j` large enough that

\[
|p^*_{w_j,m_j}-p^{ch}_{w_j}|<1/j                              \tag{7.2}
\]

and `m_j>=w_j`.

The tori `C_{w_j}xC_{m_j}` have shortest period `w_j->infinity`, so root consistency gives

\[
p^*_{w_j,m_j}\to p_c(G4).                                   \tag{7.3}
\]

Therefore

\[
\boxed{p^{ch}_w\longrightarrow p_c(G4)\qquad(w\to\infty).}   \tag{7.4}
\]

This is an iterated-limit theorem requiring no critical exponent or birth-law concentration.

## 8. Explanation of balance without concentration

In an exponentially elongated torus, the two rank births can converge to distinct constants

\[
a(d)<b(d),                                                     \tag{8.1}
\]

so the fair birth distribution does not concentrate near the root.

The root is nevertheless stable because it is controlled by a different object:

\[
\boxed{\text{root} = \text{equality of two rare charge-sector void free energies}.}\tag{8.2}
\]

At fixed large `w`, making `m` enormous suppresses both charged sectors exponentially in `m`, but the **sign of their exponent difference** remains and selects `p^{ch}_w` uniquely.  Sending `w->infinity` then moves that charge-coexistence point to the planar critical parameter.

So there is no contradiction between:

- a broad rank-one plateau / nonconcentrated birth law;
- a sharply defined matching root.

They live in different spectral sectors.

## 9. Relation to the exact charge coordinates

Recall

\[
M=\chi\tanh(\theta/2).                                       \tag{9.1}
\]

At fixed `w` and large `m`,

\[
\theta\sim m\Theta_w(p).                                     \tag{9.2}
\]

while the charged susceptibility `chi=P0+P2` is exponentially small away from the zero of `Theta_w`.

Near `p^{ch}_w`, if `Theta'_w>0`, the natural charge crossover variable is

\[
\boxed{x=m\Theta'_w(p^{ch}_w)(p-p^{ch}_w).}                  \tag{9.3}
\]

The actual limiting `chi` and neutral count law require transfer/Perron amplitudes; the **location scale** of the charge sign change is already visible from the free-energy crossing.

This is a different window from the lower/upper winding-component Gumbel windows in `w`.

## 10. A new exact-computation opportunity

For modest fixed widths, `lambda^0_{G,w}(p)` can be computed from a safe homology transfer matrix without tracking component ages.  One can therefore obtain `p^{ch}_w` as the root of

\[
\boxed{
\lambda^0_{4,w}(p)
=\lambda^0_{8,w}(1-p).}                                      \tag{10.1}
\]

This is numerically far cleaner than finding a finite-`m` root from two extremely small probabilities.

It also gives a new width sequence approaching `p_c` that is independent of the existing threshold-rank finite-torus root estimator.  Agreement of the two sequences after taking `m` large would be a strong transfer/topology cross-check.

## 11. Claim boundary

Existence of the strip/torus void free energy follows from the elementary concatenation bounds.  Analytic Perron representation and strictness are standard finite-state transfer consequences but should receive code/proof review if promoted as a production estimator.  The limit `p^{ch}_w->pc` uses the parent manuscript's author-level root-consistency theorem via a diagonal sequence; no claim of a uniform-in-`w` rate for `p^*_{w,m}->p^{ch}_w` is made.