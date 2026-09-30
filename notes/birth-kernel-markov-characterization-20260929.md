# Two-birth kernels: an exact rank-only Markov characterization

2026-09-29. Bounded author derivation for the joint-birth programme and its
three-time rank-one memory statistic. No claim of novelty, general
percolation Markov/non-Markov behaviour, continuum identification, or unique
scaling limit is made. Section 7 links the companion exact finite archive
analysis; this theorem derivation is not an independent reproduction of
those calculations.

**Main statement.** For an ordered two-birth process, the identities

    H(a,c) H(b,b) = H(a,b) H(b,c),    all a <= b <= c,          (M)

characterize ordinary time-inhomogeneous Markov behaviour in its natural
rank filtration. Direct `0 -> 2` jumps and arbitrary birth-time atoms are
allowed. A single three-time identity is necessary, not sufficient.
The continuous marginal completion below additionally assumes absolutely
continuous birth marginals and no direct double birth; those assumptions
are not needed for the characterization itself.

## 1. Process, filtration and the precise Markov contract

Let `T1<=T2` be random variables in `R union {+infinity}` and, for every
deterministic real t, set

    R(t) = 1{T1<=t} + 1{T2<=t},
    F_t^R = sigma(R(r): r<=t),
    H(s,t) = P(T1<=s, T2>t),    s<=t.

One can equivalently work on a fixed observation interval, extending the
initial all-zero process to earlier times when appropriate. We use the
raw natural filtration, or its completion by null sets; no enlargement by
hidden geometry, future information, or an encoder/random-clock state is
implicit. The right-continuous path convention is fixed by `T_i<=t`.

Write `F_i(t)=P(T_i<=t)` and

    P0(t)=1-F_1(t),   P1(t)=F_1(t)-F_2(t)=H(t,t),
    P2(t)=F_2(t).

Ordinary Markov means that for each deterministic `b<=c` there is a
stochastic kernel `K_(b,c)` on `{0,1,2}` such that

    E[g(R(c)) | F_b^R] = K_(b,c)g(R(b))   almost surely       (1.1)

for every bounded g. Null sets may depend on the deterministic times.
This note does not assert a strong Markov property at arbitrary stopping
times, or silently substitute a right-continuous enlarged filtration.

The complete theorem is: (1.1) holds if and only if (M) holds for **every**
deterministic triple with `P1(b)>0`. When `P1(b)=0`, all terms of (M)
vanish and it supplies no restriction. Endpoint cases `a=b` or `b=c`
are identities, so the substantive triples have `a<b<c`.

## 2. Necessity and sufficiency, without density assumptions

### Necessity

Fix `a<=b<=c`. On `B_b={R(b)=1}={T1<=b<T2}`, rank at c is still one
exactly when `T2>c`. If `P1(b)>0`, Markovness therefore gives the survival
probability

    K_11(b,c) = H(b,c)/P1(b).

The early-entry event

    E_a = {T1<=a, T2>b} = {R(a)=1, R(b)=1}

belongs to `F_b^R`. Integrating (1.1) over `E_a` yields

    H(a,c) = H(a,b) H(b,c)/P1(b),

which is (M). This uses neither a density nor a claim that births are
independent. A path that jumped directly from zero to two is simply not
in `B_b`.

### Sufficiency on state one: CDF events generate the observed history

Assume (M) for all the stipulated triples and fix b,c with `P1(b)>0`.
On `B_b`, for every r<=b,

    R(r)=1{T1<=r},

because `T2>b`. Thus the trace of `F_b^R` on `B_b` is generated precisely
by the sets `B_b intersect {T1<=a}`, a<=b. This is the censored first-birth
CDF information, not an additional hidden-state assumption.

Put `s_bc=H(b,c)/P1(b)`. Identity (M) says

    P(B_b, T1<=a, T2>c) = s_bc P(B_b, T1<=a)            (2.1)

for every a<=b. The CDF sets form a pi-system containing `B_b` (take a=b).
The two finite measures in (2.1) therefore agree on their generated
sigma-field, by the pi-lambda theorem/monotone class argument. Consequently

    E[1{T2>c} | F_b^R] = s_bc   on B_b, a.s.           (2.2)

It follows that the entire conditional law of `R(c)` on state one is
`(0,s_bc,1-s_bc)`. This proof includes atomic and singular birth laws;
it does not condition on an event `T1=t` of probability zero.

On `B_b^0={R(b)=0}` the whole observed history is identically zero.
The trace of `F_b^R` on this event is trivial, so the conditional future
law depends on b and the current state zero only. On `{R(b)=2}` the future
is identically two, irrespective of the possibly informative past. These
facts together with (2.2) prove (1.1). The null-risk state-one case is
irrelevant to conditional expectations.

For clarity, the on-risk kernel is explicitly

    K_00(b,c) = P0(c)/P0(b),
    K_01(b,c) = [P1(c)-H(b,c)]/P0(b),
    K_02(b,c) = 1-K_00(b,c)-K_01(b,c),                  if P0(b)>0;

    K_11(b,c) = H(b,c)/P1(b),
    K_12(b,c) = 1-K_11(b,c),   K_10(b,c)=0,            if P1(b)>0;

    K_22(b,c)=1,   K_20(b,c)=K_21(b,c)=0.

Rows at zero-probability current states may be filled arbitrarily as
probability rows; they do not affect the law started from the specified
process. Chapman-Kolmogorov follows on reachable rows. Arbitrary off-risk
fillings are not claimed to produce a canonical kernel for every artificial
starting state.

This is why the all-triples condition suffices here, although a generic
process is not certified Markov by a few three-point checks: on rank one,
the entire relevant past is generated by just the first-birth CDF events.

## 3. A readable survival contrast and its four-cell determinant

For a fixed `a<b<c`, let

    p = H(b,b),   e = H(a,b),   x = H(a,c),   z = H(b,c).

Within the state-one risk set at b, distinguish

    early: T1<=a, T2>b;
    late:  a<T1<=b<T2.

Their unconditional masses are e and p-e. If both are positive, the
conditional survival difference is exactly

    Delta(a,b,c)
      = P(T2>c | T1<=a, T2>b)
        - P(T2>c | a<T1<=b<T2)
      = x/e - (z-x)/(p-e)
      = [x*p-e*z]/[e*(p-e)].                          (3.1)

Here the four-cell table contains **unconditional** masses:

| Entry group within R(b)=1 | Survives: T2>c | Completes: b<T2<=c |
|---|---:|---:|
| Early | x | e-x |
| Late | z-x | p-e-z+x |

Its determinant, with exactly this row/column ordering, is

    D = x*(p-e-z+x) - (e-x)*(z-x)
      = x*p-e*z
      = H(a,c)H(b,b)-H(a,b)H(b,c).                    (3.2)

After conditioning the whole table on `R(b)=1`, its determinant is
`D/p^2`, also the covariance of the early-entry and survival indicators
under that conditional law. On positive early and late risks, `D=0`,
`Delta=0`, and conditional independence of these two indicators are
equivalent. A positive Delta means earlier entrants among survivors at b
are more likely to remain rank one at c; a negative value means the reverse.

If p=0 there is no rank-one comparison. If e=0 or p-e=0, Delta is undefined,
the determinant is automatically zero, and the triple is uninformative
about entry-history dependence. Do not replace an undefined contrast by
a measured zero or rely on an odds ratio with empty cells.

Conditioning on `R(b)=1` is **selection of a risk set**, not a causal
intervention setting rank to one. The criterion tests sufficiency of the
current rank for the ordinary predictive law. It does not identify a
causal effect of early birth, lifetime, geometry, or a microscopic source.

## 4. Continuous marginals: the unique no-double-jump Markov completion

For this section only, suppose T1,T2 are finite almost surely,
`P(T1<T2)=1`, and both CDFs F_1,F_2 are absolutely continuous, with densities
f_1,f_2. They are assumed compatible with such an ordered strict pair.
Their joint law need not have a two-dimensional density.

The marginal occupancies are the P_i above. On positive risks the only
possible absolutely continuous transition rates of a Markov completion are

    q_01(t)=f_1(t)/P0(t),
    q_12(t)=f_2(t)/P1(t),
    q_02(t)=0.                                       (4.1)

With `q_00=-q_01`, `q_11=-q_12`, and state two absorbing, these rates
produce the marginally matching time-inhomogeneous pure-birth completion.
Its rank-one kernel is

    H_M(s,t) = P1(s) exp[-integral_s^t q_12(r) dr].    (4.2)

The formula is first read within a connected positive-P1 interval.
It is zero when P1(s)=0, or when a zero-P1 time is crossed. Extended
integrals and risk exhaustion give the same convention, as detailed below.

### Uniqueness uses only ordinary deterministic-time Markovness

Here is a direct justification, avoiding a strong Markov assumption at T1.
On a compact interval `[s,t]` with P1 bounded below by m>0, take partitions
`s=t_0<...<t_n=t` of vanishing mesh. Define

    epsilon_i = P(t_i<T1<=T2<=t_(i+1)),
    r_i = [F_2(t_(i+1))-F_2(t_i)-epsilon_i]/P1(t_i).

The second numerator is precisely the mass that leaves rank one between
the two deterministic times having already entered by t_i. Therefore
`K_11(t_i,t_(i+1))=1-r_i`. Markovness gives

    K_11(s,t)=product_i (1-r_i).

There is no diagonal birth mass, so

    sum_i epsilon_i <= P(0<T2-T1<=mesh) -> 0.

Also `max_i r_i->0` by continuity of F_2 and the risk lower bound, while
`sum_i r_i` is bounded. Thus the logarithm of the product differs from
`-sum_i r_i` by a term tending to zero. Continuity and positivity of P1
give the Riemann-Stieltjes limit

    sum_i r_i -> integral_s^t dF_2(r)/P1(r)
                 = integral_s^t f_2(r)/P1(r) dr.

This proves (4.2) for any ordinary Markov process with the stipulated
marginals and no double jumps. Similarly the state-zero survival is
`P0(t)/P0(s)=exp[-integral_s^t f_1/P0]` on positive state-zero risk.
The on-risk kernels in section 2 are consequently fixed. Thus the
Markov **path law**, not merely an arbitrarily chosen off-risk generator,
is unique with these marginals and the no-double-jump convention.

### Existence and zero-risk boundaries

Let `O={t:P1(t)>0}`, an open set since the CDFs are continuous. On every
component `(l,r)` of O, the rate `q=f_2/P1` is locally integrable and

    P1' = f_1-f_2 = f_1-q*P1   almost everywhere.

Variation of constants, with zero occupancy at the left boundary, yields

    P1(t) = integral_l^t f_1(x)
              exp[-integral_x^t q(y)dy] dx.           (4.3)

For l=-infinity, use `P1(l)=0` as the limiting value. More formally start
at l+delta and let delta decrease to zero; the initial survival term is
bounded by P1(l+delta) and vanishes.

An ordered strict pair puts its births inside these components almost
surely. To see the only boundary issue, choose a countable dense subset
of `O^c`. Each of its points has zero rank-one probability. Hence almost
surely the open random interval `(T1,T2)` meets none of those points and
therefore none of `O^c`; it lies inside one component of O. Its endpoints
can lie outside O only at endpoints of these countably many components.
Absolute continuity gives those endpoints probability zero. Thus f_1 and
f_2 vanish almost everywhere off O; no unmodelled birth mass is placed
at a zero-risk boundary.

Construct T1 with density f_1. Given its sampled value x in O, draw its
second birth using survival `exp[-integral_x^t q]` within that component.
This defines a new process; it is not an invocation of a stopping-time
property for the original one. Local integrability gives T2>T1 almost
surely. At a right endpoint r where P1 tends to zero, the inequality

    exp[-integral_x^t q] <= P1(t)/P1(x),    x<t<r,

follows from `P1' >= -q*P1`. Thus survival falls to zero there and the
integral diverges. This also applies at r=+infinity, because both births
are finite and `P1(t)->0`. No residual atom or direct double jump needs
to be inserted at the boundary.

Equation (4.3) shows that the constructed occupancy is exactly P1(t);
its completion density is `q(t)P1(t)=f_2(t)`. Conditional survival within
state one multiplies at deterministic intermediate times, so its H is
(4.2) and satisfies (M). Section 2 then proves ordinary Markovness.

Set ratios with zero numerator and zero risk to zero for notation, or
choose any off-risk row convention; their values do not affect the
specified law. Rates may diverge near risk exhaustion, and the extended
survival integral, not a finite rate assigned at the endpoint, is what
matters. If absolute continuity, strict birth order, or compatibility
fails, this construction is not asserted. In particular an unknown
direct-double-jump flow is not determined by the two marginals alone.

**Interpretation:** fitting all static rank laws by this Markov completion
does not show that the actual process is Markov. The observed joint H
may differ from H_M while both marginals, and hence every static rank
law, agree exactly.

## 5. Finite insertion counts: marginals plus direct-jump mass

For a finite count clock, take `1<=J1<=J2<=N`,
`R_k=1{J1<=k}+1{J2<=k}`, `k=0,...,N`. Define F_i(k), P_i(k) analogously.
At step k+1 let

    m_1=P(J1=k+1),   m_2=P(J2=k+1),
    d_(k+1)=P(J1=J2=k+1).

The unconditional one-step flow table is fixed exactly:

| From / to | 0 | 1 | 2 |
|---|---:|---:|---:|
| 0 | P0(k)-m_1 | m_1-d_(k+1) | d_(k+1) |
| 1 | 0 | P1(k)-m_2+d_(k+1) | m_2-d_(k+1) |
| 2 | 0 | 0 | P2(k) |

Divide each positive-risk row by P_i(k). The resulting time-inhomogeneous
Markov chain, started at R_0=0, is the unique Markov completion matching
these marginals **and** these per-step direct-jump masses, on reachable
states. Row/column sums reproduce P_i(k) and P_i(k+1), respectively, so
existence follows directly by forward multiplication. It also reproduces
the specified d_(k+1). Its survival kernel is

    H_M(k,l)=P1(k) product_(j=k)^(l-1)
                 [1-(P(J2=j+1)-d_(j+1))/P1(j)]       (5.1)

along positive risks; it is zero once the intervening rank-one risk is
exhausted. Arbitrary off-risk rows cannot change this path law.

For compatible data all flows are nonnegative. Explicitly the allowed
direct-jump mass lies in

    max(0,m_2-P1(k)) <= d_(k+1) <= min(m_1,m_2).

If P0(k)=0 then m_1=d_(k+1)=0. If P1(k)=0 then compatibility forces
`m_2=d_(k+1)`, so no state-one exit is divided by a fictitious risk.
Static marginals without d generally leave a family of completions.
Conversely, supplying both marginals and d fixes only adjacent-time
flows of the actual process; it still does not establish its Markovness.
The CDF proof in section 2 applies on this finite time grid as well.

## 6. Minimal exact four-cell example

Take T1 in `{1,2}`, T2 in `{3,4}`, and assign the following joint masses:

| | T2=4 | T2=3 |
|---|---:|---:|
| T1=1 | 3/8 | 1/8 |
| T1=2 | 1/8 | 3/8 |

For `(a,b,c)=(3/2,5/2,7/2)`, every path is rank one at b and

    p=1,  e=1/2,  z=1/2,  x=3/8;
    D=1/8,   Delta=3/4-1/4=1/2.

The past first-birth information changes the future survival law even
though the current rank is the same. This process is not Markov in the
specified filtration. Replacing all four masses by 1/4 preserves both
birth marginals, all static rank laws and the absence of direct jumps.
The resulting independent pair (whose supports already enforce order)
has `H(s,t)=F_1(s)[1-F_2(t)]`, satisfies every identity (M), and is Markov.
This is hand algebra, not an enumeration or a percolation example.

## 7. Exact L3 archive controls from the companion analysis

The [companion analysis](birth-history-memory-readout-20260929.md) reports
an integer evaluation of all 120 strict triples
`0<=a<b<c<=9` from the already existing full `9!` paired-birth tables:

- [Square L3 paired table](../analysis/birth-gap-20260929/data/exact-control-square-L3-b00.json.gz):
  zero violations of (M) on the complete insertion-time grid.
- [Triangular L3 paired table](../analysis/birth-gap-20260929/data/exact-control-triangular-L3-b00.json.gz):
  exactly two violations, at `(3,4,5)` and `(4,5,6)`; the conditional
  early-entry/survival covariance is `1/90` for each.

These are supplied exact-archive results, not an independent re-reading of
the files in this note. The [reproducible archive script](../analysis/birth-markov-kernel-20260929/analyze.py)
does not rerun the `9!` configuration/permutation enumeration.
Applying the sufficient direction of section 2, the square result is an
**exact finite insertion-clock rank-only Markov positive control** for this
L3 ensemble. The triangular result is an exact negative control for its
corresponding ensemble. This is not a presumption that every microscopic
projection must be non-Markov, nor an extrapolation to other sizes or clocks.

For the triangular triple `(4,5,6)`, the reported count table, in the same
early/late and survives/completes ordering as section 3, is

| Entry group within R(5)=1 | J2>6 | 5<J2<=6 |
|---|---:|---:|
| Early: J1<=4 | 23,328 | 101,088 |
| Late: 4<J1<=5 | 15,552 | 93,312 |

The risk count is 233,280 out of 362,880 total permutations. The reported
conditional survivals are `3/16` and `1/7`, giving `Delta=5/112`; the
conditional covariance is `1/90`. This last number is the determinant of
the **risk-normalized** table, not the unnormalized count determinant or
the determinant of counts divided by the full 362,880 population. Production
three-quantile estimates and their uncertainties belong in the separate
results note, not in this proof's exact-control claims.

The [exact clock script](../analysis/birth-markov-kernel-20260929/exact_clock.py)
supplies a separate **exact iid-label clock transformation** of these same
tables, using multinomial order-count
weights, not new configuration enumeration. For `0<=p<=q<=1` the reported
closed kernels are

    H_square(p,q) = 6 p^3 (1-q^2)^3,
    H_tri(p,q) = 9 p^3 (1-q)^3 [A(p)+B(q)],
    A(x)=2x^3-6x^2+3x,   B(x)=-2x^3+3x+1.

The square kernel is separated and satisfies (M) for **all continuous
label-time triples**, so it is also a label-clock Markov positive control.
This conclusion comes from the explicit transformed kernel, not from
assuming that count-clock Markovness survives a random time change.
For the triangular kernel the exact factorization is

    D(a,b,c) = 81 a^3 b^3 (1-b)^3 (1-c)^3
                 [A(b)-A(a)] [B(c)-B(b)],

which is not identically zero. At `(a,b,c)=(1/3,1/2,2/3)` its reported
survival contrast is `-1156/716639`, giving a continuous-label negative
control as well. These closed forms are recorded in the companion
[polynomial certificate](../analysis/birth-markov-kernel-20260929/exact-clock.json)
and are not independently recomputed here. Their polynomial coefficient
matrices have reported algebraic ranks one and two, respectively; those
ranks are not automatically dimensions of a hidden positive state model.
The finite L3 conclusions do not establish a general lattice or scaling
limit theorem.

## 8. What a three-quantile archive comparison can and cannot say

The birth-mixture CDF is `(F_1+F_2)/2`. Its **population** .25, .5 and .75
quantiles provide deterministic times a,b,c for one application of (3.1).
They must be distinct with positive early and late risk to give an
informative survival contrast. Empirical plug-in quantiles are estimated
times; sampling uncertainty and their dependence on the same archive do
not disappear by calling them fixed quantiles.

A nonzero population determinant excludes rank-only ordinary Markovness
for this time/observer/ensemble contract. A finite estimate requires its
statistical uncertainty to be assessed. A zero value at one triple cannot
establish the all-triples condition; it may miss dependence at other
times. A finite-N violation need not persist in a scaling limit, and
neither outcome proves uniqueness of a continuum birth copula. Selection
of rank one remains predictive conditioning, not an intervention.

**Birth age is a history coordinate, not a discovered mechanism.** On
rank one, the observed past is already encoded by T1. Thus augmenting
the state by the observed birth time, or by age `t-T1` while retaining
calendar time t, supplies a complete predictive state for this two-birth
rank process. On rank zero there is no extra observed history; on rank
two the future is absorbing. This gives a time-inhomogeneous ordinary
Markov description of the augmented process, regardless of whether rank
alone is Markov. It does not automatically give a time-homogeneous
semi-Markov law: the completion holding-time distribution may still
depend on the calendar entry time (with direct-jump routing handled
separately). Finding prediction improves after adding age is therefore
not itself a microscopic mechanism result or a reason to search more
age descriptors. A mechanism would explain how T1, through geometry,
direction or other specified marks, changes the completion law.

Finally, the count clock k is not the finite-N iid-label clock p. For iid
continuous labels, let `K(p)` be the number of sites with labels <=p;
the ordering permutation is independent of the order statistics. Then

    R_label(p)=R_count(K(p)),
    H_label(s,t)=P(J1<=K(s), J2>K(t)).

At one time, static laws are binomial mixtures over k; at several times,
the count increments have the corresponding joint multinomial law.
The random hidden count can affect transition probabilities even when
the count-clock rank has a Markov completion. Neither the measured
three-time contrast nor the Markov property transfers automatically by
substituting k/N for p. A deterministic invertible time reparameterization
is a different matter and preserves the ordinary Markov property on its
time image. Clock and observer definitions must accompany the result.
