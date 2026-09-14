# A square-SITE near-critical renewal programme: smaller than full universality

2026-09-14.  Proof programme distilled after auditing D'Alimonte--Manolescu's 2026 near-critical OZ construction and the classical Bernoulli near-critical toolbox.

The strategic conclusion is that Matching-One does **not** need a proof of conformal universality for square-site percolation in order to obtain a useful near-critical OZ/renewal theorem.  The independent-SITE model removes some of the hardest random-cluster boundary-condition issues.  A comparability-level renewal theorem appears to require a much smaller set of inputs.

Nothing in this note is promoted to an accepted theorem.  Each missing gate is isolated explicitly.

## 1. Why SITE may be simpler than the bond-FK proof

In a random-cluster exploration, conditioning on the explored past changes the law in the future through induced wired/free boundary conditions.  D'Alimonte--Manolescu therefore need robust RSW under boundary conditions and a nontrivial coupling/mixing construction to obtain a killed Markov renewal process with uniform mass gap.

For Bernoulli SITE percolation:

- site states are independent;
- after revealing an arbitrary explored set, all unrevealed sites remain iid Bernoulli(`p`);
- no wiring information propagates through the unexplored region;
- finite-energy ratios are explicit powers of `p/(1-p)`.

Thus the difficult “future forgets the past” step should reduce to a **geometric separation/barrier** statement rather than a measure-comparison theorem.

The matching 4/8 convention provides exactly the site-dual geometry needed to construct white barriers around black NN clusters.

## 2. Characteristic length without critical exponents

For black NN site percolation define a characteristic length `L(p)` below criticality by a fixed rectangle-crossing threshold.  One possible convention is the smallest `L` such that the long-direction black crossing probability of a fixed-aspect rectangle is below a small constant `epsilon_0`.

The desired properties are only up to bounded factors:

\[
L(p)\to\infty\quad(p\uparrow p_c),                            \tag{2.1}
\]

and for rectangles/quads of diameter at most `cL(p)`, black and complementary white-matching crossing probabilities are uniformly nondegenerate.

Critical square-site RSW on `Z^2` is available in the literature despite the lack of self-duality.  The remaining near-critical extension is a finite-size criterion/RSW stability question, not an exponent-identification question.

No use of `nu=4/3`, `5/48`, SLE, or conformal invariance is needed for the renewal architecture.

## 3. Gate A: sub-characteristic RSW for the 4/8 pair

**Target A.**  There exist fixed constants `c,C,epsilon>0` such that, uniformly for `p<pc`, every topological rectangle of bounded aspect ratio and diameter at most `cL(p)` satisfies

\[
\epsilon
\le P_p(\text{black NN crossing})\le1-\epsilon,       \tag{3.1}
\]

with the complementary statement for white matching crossings.

A version stable under conditioning on already revealed sites outside the rectangle is automatic once the rectangle's site set is disjoint from the revealed set, because SITE variables inside remain iid.  This is much simpler than FK arbitrary-boundary-condition RSW.

The exact finite rectangle matching dichotomy should be used to transfer missing black crossings to white matching transverse crossings.

**Acceptance gate:** prove A using the site's critical RSW plus the chosen definition of `L(p)` and standard gluing.  If a particular definition of `L` makes one side awkward, change the definition; only bounded equivalence to the correlation length is needed.

## 4. Gate B: one-arm stability to L(p)

D'Alimonte--Manolescu's endpoint factor is the critical one-arm probability at correlation scale.  The SITE analogue needs

\[
\boxed{
P_p(0\leftrightarrow\partial B_r)
\asymp
P_{p_c}(0\leftrightarrow\partial B_r),
\qquad r\le cL(p),}                                          \tag{4.1}
\]

uniformly near criticality.

Classical Kesten near-critical theory proves this type of arm stability for planar Bernoulli percolation in standard settings; the exact square-site/matching formulation should be checked line by line rather than cited by universality.

This is likely the most literature-sensitive gate of the programme.

**Minimal substitute:** for #767 comparability, exact equality to the critical arm is not essential.  It would already suffice to define

\[
\pi_{1,site}(p):=P_p(0\leftrightarrow\partial B_{L(p)})       \tag{4.2}
\]

and retain it as a measured/controlled endpoint insertion.  The renewal theorem can be stated with `pi_{1,site}(p)` directly.  Critical replacement can be a later corollary.

## 5. Gate C: correlation-length slab barriers and killing

Fix a direction `w`.  Explore the black cluster through slabs of longitudinal thickness

\[
M L(p)                                                        \tag{5.1}
\]

for one large universal `M`.

Using Gate A, require in each coarse block a finite collection of complementary white-matching crossings that isolates the surviving black channel and prevents long memory around the block boundary.

Because each such white barrier has probability bounded below uniformly, one obtains:

- a uniformly positive probability that the black exploration dies in the next coarse block;
- a uniformly positive probability of one clean surviving channel;
- exponentially small probability of an anomalously large irreducible piece after iterating blocks.

This is the SITE analogue of a uniform killing rate and mass gap.

Unlike FK, once the barrier site's states are revealed, the unexplored future is still iid.  No boundary-condition coupling is required.

**Target C.**  Construct a coarse exploration state `Y_k` such that the longitudinal/transverse increments form an aperiodic killed Markov-additive process with

\[
0<c\le\kappa_{kill}\le1-c,                                   \tag{5.2}
\]

and exponential tails of the irreducible block length in units of `L(p)`, uniformly in `p` and direction.

## 6. Gate D: nondegenerate transverse step variance

A local CLT needs the transverse increment variance to be bounded above and below in `L(p)` units.

The upper bound follows from exponential tails in Gate C.

For the lower bound, build two finite crossing patterns of uniformly positive probability that move the clean channel transversely by distinct `O(L(p))` amounts over one or a bounded number of renewal blocks.

Gate A supplies these patterns.

Thus the expected target is

\[
\boxed{
c\le \operatorname{Var}(\Delta X/L(p))\le C}               \tag{6.1}
\]

uniformly near criticality and in direction.

This is geometric and should not require arm exponents.

## 7. Consequence of A--D: a SITE killed-renewal OZ comparability

Classical killed-renewal theory would then give, for point-to-point connections at physical distance `r>=CL(p)`,

\[
P_p(0\leftrightarrow re)
\asymp
J_{left}(p,e)J_{right}(p,e)
\left(\frac r{L(p)}\right)^{-1/2}
\exp[-r/\xi_p(e)],                                            \tag{7.1}
\]

with

\[
\xi_p(e)\asymp L(p).                                         \tag{7.2}
\]

Gate B would identify the endpoint insertions up to constants as

\[
J_{left}J_{right}\asymp\pi_{1,site}(L(p))^2.                 \tag{7.3}
\]

The same renewal process gives

\[
\boxed{D_p(e)\asymp\xi_p(e)}                                 \tag{7.4}
\]

and a Brownian bridge for the conditioned core.

This already supplies nearly every near-critical structural input needed by #740/#758/#767, without any exact amplitude.

## 8. Gate E: from an open connection to a complete winding component

This is the genuinely Matching-One-specific step and should be kept separate from point-to-point OZ.

On a cylinder of circumference `w=sL(p)`, define a clean renewal core that returns to its starting transverse row after one horizontal period.  A cyclic renewal calculation predicts a closure density

\[
\asymp L(p)^{-1}s^{-1/2}e^{-cs}.                              \tag{8.1}
\]

What remains is to compare this clean cyclic core to the **actual complete SITE component anchor activity**, including branches and the unique-anchor convention.

A realistic first theorem is only comparability:

\[
\boxed{
\nu_w(p)
\asymp
L(p)^{-1}J_{comp}(p,s)s^{-1/2}e^{-w/\xi_p},}                 \tag{8.2}
\]

with `J_comp` bounded above/below on a declared range of `s`, or with its endpoint/mark dependence left explicit.

This is strictly weaker than identifying `zeta(p)`, but it is enough to determine which logarithmic terms enter the crossover centre.

## 9. A possible direct route through the existing tagged component construction

The repository already has an exact microscopic tagged complete-component resolvent at each fixed width.  Its weakness near criticality is state explosion when `w~xi->infinity`.

The renewal programme suggests a renormalized version:

1. group `O(L(p)) x O(L(p))` microscopic boxes into coarse blocks;
2. classify only the finite set of clean crossing/renewal states needed by Gate C;
3. carry the component's winding and unique-anchor mark at coarse scale;
4. integrate all microscopic branches inside each block into transition weights;
5. apply the existing matrix-sewing/logdet analysis to the coarse operator.

If successful, the number of coarse states can remain `O(1)` as `p->pc`, while one winding uses `s=w/L(p)` coarse steps.

This is a more plausible near-critical theorem engine than increasing the microscopic transfer width.

## 10. Why complete conformal universality is unnecessary

The programme uses only:

- critical/sub-characteristic RSW;
- finite-size correlation length;
- arm stability or an explicit near-critical endpoint arm observable;
- independence of unexplored SITE variables;
- BK/Harris;
- classical local CLT for a killed Markov-additive process.

It does not need:

- identification with SLE6;
- exact critical exponents;
- conformal covariance of square-site scaling limits;
- a proof that square-site and bond-FK amplitudes agree.

So the correct research question is not “can we prove universality for square-site percolation?” but rather “can we build a uniform correlation-length-scale renewal decomposition for independent SITE?”

## 11. Suggested proof order

The shortest dependency chain is:

1. freeze a precise `L_site(p)` convention and prove Gate A;
2. build the clean white-matching barrier / black survivor block and prove Gate C;
3. prove Gate D and invoke killed-renewal local CLT;
4. first state the point-to-point theorem with the near-critical endpoint insertion left as `pi_{1,site}(p)`;
5. only then prove/quote Gate B to replace it by the critical one-arm probability;
6. finally attack Gate E for complete component activity.

This order avoids making the hardest SITE-specific insertion issue a prerequisite for the renewal skeleton itself.

## 12. Claim boundary

Critical RSW for square-site percolation on `Z^2` is a published theorem.  Classical near-critical arm/correlation-length stability exists for planar Bernoulli percolation, but the exact square-site/matching formulations required above have **not** been independently certified in this note.  Gates A--E remain a proof programme until those inputs and constructions are written out.