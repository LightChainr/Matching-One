# Completion-pivotal geometry transmits to the two-birth kernel

2026-09-29. Bounded author derivation, continuing the
[birth-kernel Markov characterization](birth-kernel-markov-characterization-20260929.md).
This note supplies an exact finite-model link

    occupied geometry -> completion-pivotal count -> exit hazard -> H.

It does not certify a large-size limit, independent births, an original-U
operator mechanism, or a recursively sufficient compressed geometry state.
The last section quotes the [L3 geometric calculation](birth-completion-geometry-20260929.md);
no archive, subset enumeration, simulation, or test is run here.

## 1. Finite geometry and the two clocks

Let V be a fixed set of N sites, and let
`r:2^V -> {0,1,2}` be a monotone rank function, with `r(empty)=0` and
`r(V)=2`. In the intended application this is the occupied graph's ambient
homology rank on a fixed torus. The argument uses monotonicity and the
specified sampling law, not a critical exponent or a topology-specific
formula. A single insertion may increase rank directly from zero to two.

For an occupied set A of rank one, define

    nu_2(A) = #{v in V\A : r(A union {v})=2}.           (1.1)

These are the currently empty sites whose individual insertion completes
rank two. Set nu_2=0 outside rank one when writing indicator-weighted
expectations; no exit interpretation is assigned to that extension.

We distinguish two experiments throughout:

- **Count clock.** A uniformly random site permutation gives occupied
  prefixes A_k, |A_k|=k, and births `J_i=min{k:r(A_k)>=i}`.
- **Label clock.** Independent Uniform[0,1] site labels give
  `A_t={v:U_v<=t}` and `T_i=inf{t:r(A_t)>=i}`.

Write

    H_count(a,k)=P(J1<=a, J2>k),     a<=k,
    H_label(s,t)=P(T1<=s, T2>t),     s<=t.

We use H without its subscript only when the clock is explicit. The count
denominator `N-k` is deterministic. In label time the number of remaining
sites is random, and each has its own conditional hazard `1/(1-t)`.
One does not obtain the label formulas by substituting `k=Nt`.

## 2. Count-clock exit and the exact H recursion

Let G_k contain the complete revealed permutation prefix, not merely its
rank. Conditional on G_k, the next site is uniform on the N-k empty sites.
For k<N, on `{r(A_k)=1}`,

    P(J2=k+1 | G_k) = nu_2(A_k)/(N-k).                (2.1)

Fix a<=k. The event `C_a(k)={J1<=a,J2>k}` is G_k-measurable and lies
inside rank one. It loses mass at the next step exactly when a completion
site is chosen. Therefore

    H_count(a,k+1) = H_count(a,k)
      - E[1{J1<=a,J2>k} nu_2(A_k)]/(N-k).            (2.2)

This remains an undivided identity when the cohort has zero probability.
For positive cohort risk, put

    m_a(k)=E[nu_2(A_k) | J1<=a,J2>k].

Its one-step conditional survival is

    H_count(a,k+1)/H_count(a,k) = 1-m_a(k)/(N-k).      (2.3)

Thus two entry cohorts at the same k can have different exit probabilities
only through their different conditional averages of the current
completion count. The full geometry determines that count pointwise;
the entry history changes the distribution of geometries in the risk set.

For disjoint early and late groups `J1<=a` and `a<J1<=b`, both surviving
to b, their **next-step** survival contrast is exactly

    survival_early - survival_late
      = [E(nu_2 | late,J2>b)-E(nu_2 | early,J2>b)]/(N-b).   (2.4)

This is a predictive risk-set comparison, not an intervention on birth
time or geometry. Equality of these cohort averages can occur even when
nu_2 varies substantially across individual rank-one configurations.

When every relevant survival factor is positive, the finite-clock
multiplicative defect also has the exact sum representation

    log[H(a,c)H(b,b)/(H(a,b)H(b,c))]
      = sum_(k=b)^(c-1)
          log[(1-m_a(k)/(N-k))/(1-m_b(k)/(N-k))].     (2.5)

If a factor is zero, use (2.2) and the undivided determinant rather than
assigning a finite logarithm. No interpolation of (2.5) is being assumed.

## 3. Uniform-label completion intensity and the kernel derivative

Let G_t contain the occupied history and the labels already revealed by t,
but no future labels. Conditional on G_t, all empty sites have independent
Uniform(t,1) residual labels. Each empty site therefore arrives in
`(t,t+h]` with probability `h/(1-t)` and has instantaneous hazard
`1/(1-t)`. The independence here concerns the remaining site clocks,
**not** the two topological births.

For a rank-one configuration, completion in the next short interval occurs
to first order exactly when one of its nu_2 completion sites arrives.
Any discrepancy involves at least two arrivals and is
`O_N(h^2/(1-t)^2)`. Consequently, conditional on G_t and rank one,

    P(T2<=t+h | G_t)
      = nu_2(A_t) h/(1-t) + O_N(h^2/(1-t)^2).        (3.1)

Fix `s<t<1`. The old cohort `C_s(t)={T1<=s,T2>t}` is already in rank
one at t; it cannot gain new members as t increases. Conditioning (3.1)
and taking expectations gives

    partial_t H_label(s,t)
      = -E[1{T1<=s,T2>t} nu_2(A_t)]/(1-t).           (3.2)

This derivative holds throughout the interior s<t<1, not just formally:
for a finite model the two-time configuration probabilities are finite
sums of terms in s, t-s and 1-t, so H is polynomial on that triangle.
Alternatively the uniform short-interval bound gives local absolute
continuity and the same derivative. A diagonal component `T1=T2` in the
joint birth law does not invalidate this off-diagonal argument.

For positive H(s,t), the cohort's logarithmic completion hazard is

    -partial_t log H_label(s,t)
      = E[nu_2(A_t) | T1<=s,T2>t]/(1-t).             (3.3)

Equivalently this is the hazard of the conditional survival
`H(s,t)/H(s,s)`, when H(s,s)>0. The cutoff s is held fixed.

### Do not differentiate the moving diagonal as a fixed cohort

Along `H(t,t)=P1(t)`, new rank-one entrants arrive while old ones complete.
Let gamma_01(t), gamma_12(t), gamma_02(t) be the unconditional probability
fluxes for the indicated rank transitions. The formula above gives

    gamma_12(t)=E[1{r(A_t)=1}nu_2(A_t)]/(1-t),
    P1'(t)=gamma_01(t)-gamma_12(t),
    F_2'(t)=gamma_12(t)+gamma_02(t).                  (3.4)

The other fluxes similarly count the individual empty sites causing their
specified transition, multiplied by `1/(1-t)`. Thus neither
`-d log P1(t)/dt` nor, when direct jumps occur, `F_2'(t)/P1(t)` is the
rank-one completion hazard. The latter coincides with the mean rank-one
exit intensity only when the direct `0 -> 2` flux is absent.

Continuous independent site labels are distinct almost surely, but a
single site can still cause `T1=T2`. Such direct double jumps never enter
the rank-one cohorts in (2.2) or (3.2); they must not be silently included
in their completion count.

## 4. The Markov defect is an integrated cohort hazard difference

Take `0<=a<b<c<1` with all four H values below positive, and define

    L(a,b,c)=log[H(a,c)H(b,b)/(H(a,b)H(b,c))],
    m_s(t)=E[nu_2(A_t) | T1<=s,T2>t].

Both conditional survival ratios start at one at t=b. Integrating (3.3)
therefore proves the exact identity

    L(a,b,c) = -integral_b^c [m_a(t)-m_b(t)]/(1-t) dt.   (4.1)

In particular, higher early-cohort completion counts throughout a time
interval lower its relative survival. A single zero defect can also
result from cancellation of hazard differences over that interval; it
does not establish history independence at all intermediate times.

**The second cohort in (4.1) is birth-by-b, not all rank-one sites at t.**
Its event is `T1<=b,T2>t`, with b fixed. It excludes paths that first enter
rank one between b and t. Replacing it by `{T1<=t<T2}` changes the
comparison and breaks the identity.

To relate it to the early/late contrast, let

    alpha(t)=H(a,t)/H(b,t),
    m_late(t)=E[nu_2(A_t) | a<T1<=b,T2>t].

When the late cohort also has positive risk,

    m_b(t)=alpha(t)m_a(t)+(1-alpha(t))m_late(t),
    m_a(t)-m_b(t)=(1-alpha(t))[m_a(t)-m_late(t)].      (4.2)

The weights are survivor weights at t, not the original weights at b.
This identifies exactly how different geometric mixtures are transmitted
into the observable two-time kernel.

The earlier Markov theorem states that the undivided H identity must hold
for **all** deterministic triples to give ordinary rank-only Markovness.
Equation (4.1) is its finite-geometry mechanism interface, not a replacement
by one fitted hazard or one triple. In label time, positive cohort risk at
b remains positive for every c<1: the conditional event of no further site
arrivals has probability at least `((1-c)/(1-b))^N`. Thus interior logs do
not hide an exhausted positive cohort. At a zero initial risk or at the
endpoint 1 use the undivided equations or justified one-sided limits.

## 5. Persistent direction marks and symmetry

Let D(A) be a mark on rank-one configurations, taking values in a finite
set of sectors. Assume it stays unchanged along every occupied-set
inclusion path that remains rank one. This is a real hypothesis: a generic
geometric descriptor need not be persistent. For an ambient homology rank,
the unoriented one-dimensional homology image over Q is one natural
persistent direction, because those images are nested under inclusion.
An arbitrary chosen winding representative or other geometric statistic
must not be assumed constant for this reason.

Define the **unconditional sector mass**

    H_d(s,t)=P(T1<=s,T2>t,D(A_t)=d).

Its birth-by-s cohort remains in sector d until completion, so exactly the
same argument gives

    partial_t H_d(s,t)
      = -E[1{T1<=s,T2>t,D(A_t)=d}nu_2(A_t)]/(1-t),  (5.1)

and, for the count clock,

    H_d(a,k+1)=H_d(a,k)
      -E[1{J1<=a,J2>k,D(A_k)=d}nu_2(A_k)]/(N-k).   (5.2)

Conditional sector logarithms and the defect integral follow by replacing
H and m with H_d and the corresponding within-sector conditional means.
If D can change while rank remains one, inter-sector transfer terms are
needed; (5.1)--(5.2) cannot then be used without modification.

Suppose a group acts on the sites and directions, preserves r and the
sampling law of complete labelled/permuted histories, and satisfies
`D(gA)=gD(A)`. The same bijection gives `nu_2(gA)=nu_2(A)`. If the action
is transitive on m admissible directions, then, for either clock,

    H_d(s,t)=H(s,t)/m,
    E[1{cohort,D=d}nu_2]=E[1{cohort}nu_2]/m.         (5.3)

Consequently each sector has the same conditional completion mean as
the pooled cohort, the same early/late survival contrast, and the same
logarithmic Markov defect. The sector's **unnormalized** determinant is
the pooled determinant divided by m^2; its risk-normalized conditional
covariance is unchanged. Under these exact symmetry hypotheses, splitting
by pure direction cannot remove the entry-history defect.

When there are several symmetry orbits, equality is guaranteed only within
an orbit. Anisotropic sources or a geometry without the required symmetry
need not satisfy (5.3). Also H_d is not a probability conditional on ever
visiting sector d: that conditional kernel would divide by the visitation
probability. Direct double-jump paths visit no rank-one sector, so that
normalization must not count them as though they had a direction.

In particular, unmarked Markovness does **not** imply Markovness after
revealing a persistent direction. The pooled kernel is `H=sum_d H_d`,
but its determinant is not the sum of sector determinants: cross-sector
terms matter. Without transitivity on the whole direction set, averaging
can hide entry-history dependence that remains within a revealed sector.
Conversely, adding the mark supplies more information, not less; it is
the Markov sufficiency property, rather than prediction quality, that need
not be preserved by this change of observer.

## 6. What nu_2 does and does not summarize

At a fixed occupied set it is exactly sufficient to determine the
**immediate completion hazard** under these clocks. It is not thereby a
recursively sufficient state for future predictions. Two configurations
with the same nu_2 can have different counts of non-completing successor
configurations with each possible future nu_2, or different mark transitions.
Those successor laws matter for longer survival and are absent from (1.1).
A claimed `(k,nu_2,D)` Markov compression would require a separate closure
argument for these transition distributions.

If nu_2 is a function of k alone on every accessible rank-one configuration
in the count model, then rank-one exit is indeed independent of all past
history. Together with the all-zero past in state zero and absorption of
state two, this suffices for count-clock rank-only Markovness. The converse
does not require pointwise constancy: suitable cohort averages can agree
even when microscopic counts differ. In label time, even a function of
the hidden count |A_t| need not reduce to a function of t and current rank.

These identities therefore provide a genuine predictive transmission law,
but not a causal intervention, a proof that the births are independent,
or an explanation of the original normalized U observable with its own
source, normalization and moving-root contract. No large-L behaviour or
scaling exponent follows from them alone.

Uniformity of the sampling law is also substantive. A biased permutation
requires its actual conditional probabilities of choosing each empty
site. For independent nonuniform site-label laws with hazards h_v(t),
the instantaneous completion intensity is the sum of h_v(t) over the
completion sites, not nu_2/(1-t). Correlated or history-controlled sources
require the corresponding conditional intensities and may invalidate the
remaining-site independence used here.

## 7. Supplied L3 illustration: geometric mixtures, not an age mechanism

The parent reports the following exact results from a new **512-subset
geometric dynamic programme** for N=9. The designated companion is
[birth-completion-geometry-20260929.md](birth-completion-geometry-20260929.md),
where the parent supplies the finite-mechanism derivation, exact code and
count-clock realization. These numbers are quoted, not recomputed or
independently certified in this note.

- On the reported square rank-one layers `k=3,4,5,6`, every configuration
  has `nu_2=k-3`. Thus entry cohorts on each such layer have identical
  instantaneous count-clock exit probabilities. The general all-risk-layer
  sufficiency statement is the one in section 6.
- Triangular rank-one configurations have `nu_2 in {2,3}` at k=4 and
  `nu_2 in {3,4}` at k=5. Mere variation does not by itself establish a
  rank-only memory defect; the cohort mixture is the decisive quantity.
- The parent finds both previous triangular Markov-defect triples still
  present separately in each of the three reported direction sectors
  `(1,0),(0,1),(1,1)`. This is consistent with the symmetry mechanism in
  section 5; pure direction is not the missing predictive coordinate there.

For the triangular comparison at b=5, early entry `J1<=4` has reported
mean `E nu_2=13/4`, whereas late entry `J1=5` has mean `24/7`, conditional
on survival to b in each case. The weight of the faster `nu_2=4` class
is respectively `1/4` and `3/7`; the other class has nu_2=3. With four
empty sites left, (2.4) directly gives

    survival_early = 1-(13/4)/4 = 3/16,
    survival_late  = 1-(24/7)/4 = 1/7,
    survival_early-survival_late = 5/112.

This is the already reported one-step contrast explained through a
specified microscopic exit mechanism: the later-entry risk set contains
more configurations with four completion sites. It is not evidence that
changing birth age itself causes an exit change. It does not establish
that nu_2 alone closes subsequent dynamics; the parent's separate finite
state construction addresses that additional question.

The parent's subsequent exact marked-kernel calculation gives a particularly
useful **clock/observer boundary**. For square L3, all direction sectors
are Markov in insertion time. In iid-label time the unmarked rank kernel
is still Markov, but each axial sector has a strictly negative determinant
for `0<a<b<c<1`; each diagonal sector remains Markov. Axial and diagonal
directions are different symmetry orbits, so (5.3) never licensed dividing
the unmarked square kernel equally among all directions. For triangular
L3 the three marked kernels instead satisfy `H_d=H/3`, consistently
retaining the history defect under pure direction conditioning.

The companion note is the source of these closed forms and their exact
certificates; they are not rederived here. Their separation/coefficient
matrix ranks are not automatically dimensions of a positive hidden-state
model. In particular, neither count-clock sector Markovness nor unmarked
label-clock Markovness can substitute for checking the actual marked
label-clock contract. None of these finite L3 classifications is a
large-L conclusion.
