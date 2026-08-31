# P334: uniform random blockade is an exactly closed clock semigroup

An entire unmarked birth clock already determines the mean response to
**every strength of independent uniform site blockade**, and to a uniformly
chosen fixed number of blocked sites. Such averaged interventions do not add
geometric information to that clock. Spatially marked responses and their
replica correlations can add information: the double-star/C4 counterexample
separates these two kinds of observables exactly.

This is an elementary finite-state result for any monotone event, not a
large-size approximation, a torus-specific assumption, or a new Monte Carlo
measurement. It explains both a useful analytic closure and a precise
identifiability boundary.

## The clock and the intervention

Let the fixed prefix have d selectable sites. Write f_j for the number of
safe j-subsets, F(z)=sum_j f_j z^j, and

\[
 S(j)=P(T>j)=f_j/\binom dj.
\]

Permanently block each site independently with probability q=1-a, but retain
its original insertion label as an inert dummy. Then draw an independent
uniform order on all d original labels. Expectations below average both the
blockade mask and this order; a specified, fixed spatial mask need not obey
the averaged formulas.

After k labels, the number J of active sites is Binomial(k,a). Conditional
on J=j, their set is a uniform j-subset of the original d sites. Hence

\[
 \boxed{S_a(k)=\sum_{j=0}^k\binom kj a^j(1-a)^{k-j}S(j).}
\]

Multiplying by binomial(d,k) and collecting powers gives the equivalent
polynomial identity

\[
 \boxed{F_a(z)=[1+(1-a)z]^d
 F\!\left(\frac{az}{1+(1-a)z}\right).}
\]

The prefactor uses the **original number d of selectable labels**, including
inert or redundant ones; it is not the degree of F. Coefficients of F_a are
averaged safe counts and can be nonintegral.

In an independent remaining-site occupancy experiment with occupation
probability u, survival becomes simply S_cont(a u), where

\[
 S_{\rm cont}(u)=\sum_j f_j u^j(1-u)^{d-j}.
\]

Here u concerns the d remaining sites. It must not be silently identified
with the project's full-N canonical p after a prefix at nonzero k0. That
canonical response is obtained by applying its original binomial kernel to
the transformed clock.

## A pure-death semigroup on the complete clock vector

Let B_a denote the binomial transform above. Independent successive
blockades multiply active probabilities:

\[
 B_a B_b=B_{ab},\qquad B_1=I.
\]

With a=exp(-t), the exact generator acting on a clock vector is

\[
 (L S)(k)=k[S(k-1)-S(k)],\qquad (LS)(0)=0.
\]

This is the usual finite pure-death generator on the count of retained
labels. It is a linear closure for **the intervention-averaged complete
clock**, not a statement that the scalar state `(k,H2,b2)` closes the physical
growth process, or that individual spatial interventions have no memory.

At q=0 its first derivative is

\[
 \partial_q S_{1-q}(k)|_{q=0}=kP(T=k)
 =\sum_v\Delta S_v(k),
\]

recovering the marked singleton identity in
`notes/p334-exact-knockout-response-identity.md`. If any blockade prevents
eventual birth, the true mean waiting time is infinite. All formulas for
survival, or the clock censored at d+1, remain finite and valid; no finite
unconditional mean is inferred in that case.

For exactly b uniformly chosen blocked sites instead of an independent mask,

\[
 S^{[b]}(k)=\sum_j
 \frac{\binom{d-b}{j}\binom b{k-j}}{\binom dk}S(j),
\]

where out-of-range binomial coefficients vanish. The one-block case is
`S^[1](k)=((d-k)S(k)+kS(k-1))/d`.

## The five-site example survives every uniform blockade mean probe

The previously constructed double-star and C4 plus one inert site have

\[
 F(z)=1+5z+6z^2+2z^3,\quad
 S=(1,1,3/5,1/5,0,0).
\]

They therefore remain identical under every B_a and every fixed-b average,
not merely without an intervention. For instance, at a=q=1/2 both have

\[
 S_{1/2}=(1,1,9/10,3/4,47/80,7/16).
\]

The final survival probability7/16 is the chance that the retained active
set never contains a triggering pair. Their finite censored mean is187/40;
their true mean is infinite. Neither the blockage dose-response nor this
averaged terminal probability distinguishes the two geometries.

In contrast, for two independent continuations of the **same fixed state**,
the final-site collision probabilities are5/18 and1/4. The distributions of
single-site mean-knockout effects differ as well: their normalized squared
concentrations are113/392 and1/4. These are spatial response statistics,
not another linear transform of the unmarked clock.

## Consequence for the next data stream

Averaging more uniformly chosen blockade masks cannot distinguish
clock-equivalent mechanisms at the mean level, even with the entire dosage
curve. Three observations can move outside that closed family:

1. a named site's response or a prescribed spatially nonuniform blockade;
2. the final insertion site, with conditional independent replica overlap;
3. the variation of conditional responses across blockade masks, retaining
   the same-mask dependency instead of averaging it away first.

The scope is exact but specific: equal clocks imply equal **mean** uniform
blockade responses. Higher mask-conditioned moments are not claimed equal.
The147 existing real-prefix clocks are all distinct (a11d6499); the
five-site counterexample is a constructed mechanism example, not a claim
that an unmodified equal-clock pair was observed in that archive.

## Source and lifecycle

- Source event class: arbitrary finite monotone birth event; original d clock.
- Concrete example:250c589958fc09c52380feb4c99276c8e9c9455b,
  `notes/p334-isoclock-marked-geometry-counterexample.md`.
- Prior singleton relation:06845ff6b651603e36412f774703483e43c958f0.
- Result type: exact algebraic derivation and rational example, zero new
  samples, zero network solves, no new independent empirical block.
- Next discriminating observation: a spatial mark or same-prefix replica
  correlation, rather than a finer mean uniform-blockade curve.
