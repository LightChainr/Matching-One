# Directional Markov completion and its exact count-clock criterion

2026-09-29. Bounded author derivation for the existing finite-site model.
This is a concrete model interface for the
[completion-pivotal hazard identity](completion-pivotal-hazard-kernel-20260929.md)
and the [two-birth Markov criterion](birth-kernel-markov-characterization-20260929.md),
not a novelty claim or an analysis of production data. The designated empirical
companion is [directional-hazard-contrast-20260929.md](directional-hazard-contrast-20260929.md).
No new source, observable, cutoff scan, simulation, or test is introduced here.

**Main statements.** Actual directional one-step probability flows always
define a time-inhomogeneous Markov completion preserving every marked
single-time law and every unconditional adjacent-time joint table. The actual
marked process has that Markov law if and only if, at every count and within
every direction, the risk-set mean completion count is independent of every
first-entry cutoff. A zero survival contrast over just one window does not
give this condition.

## 1. Observer, probability law, and persistence

Let V have N sites, and let `r:2^V -> {0,1,2}` be monotone, with
`r(empty)=0` and `r(V)=2`. A uniform random permutation gives the occupied
prefix A_k at count k, for `k=0,...,N`. Set

    J_i = min{k:r(A_k)>=i},     i=1,2.

Direct `0 -> 2` jumps, hence `J1=J2`, are allowed. On rank-one sets let
D(A) take values in a finite set of directions, and assume persistence:

    A subset B, r(A)=r(B)=1  =>  D(A)=D(B).

The primitive unoriented ambient homology direction has this property when
the rank and direction are defined by the nested homology images. A general
geometric descriptor need not have it.

The complete *current marked state* under consideration is

    X_k = 0                  if r(A_k)=0,
          (1,D(A_k))         if r(A_k)=1,
          2                 if r(A_k)=2.

Its natural observation filtration is `F_k^X=sigma(X_0,...,X_k)`. It does
not additionally reveal the occupied set, the insertion order, or nu_2.
The complete microscopic prefix filtration is denoted G_k. Write `d` as
shorthand for the state `(1,d)` when indexing probabilities. Once a path
leaves d it is absorbed in 2; different rank-one directions cannot succeed
one another on the same path. No direction is assigned to a direct double
jump. The current state 2 does not retain the old mark, although the
observation history can retain it; this is harmless because its future is
constant.

For rank-one A define the existing completion count

    nu_2(A) = #{v in V\A:r(A union {v})=2}.

Uniform permutation sampling gives, on `{X_k=d}` and for `k<N`,

    P(X_(k+1)=2 | G_k) = nu_2(A_k)/(N-k).             (1.1)

All claims below concern this preparation and deterministic **count** times.
They are not strong lumpability statements for arbitrary microscopic initial
preparations.

## 2. The flow-defined directional Markov completion

Let `p_x(k)=P(X_k=x)`. Define the actual, unconditional one-step flows

    alpha_d(k) = P(X_k=0, X_(k+1)=d),
    beta_d(k)  = P(X_k=d, X_(k+1)=2),
    gamma(k)   = P(X_k=0, X_(k+1)=2).

In particular,

    beta_d(k) = E[1{X_k=d} nu_2(A_k)]/(N-k).          (2.1)

The entry and direct-jump flows also have microscopic counting forms:
count empty sites taking a rank-zero A respectively to `(1,d)` or to 2,
multiply by `1{r(A_k)=0}`, take expectations, and divide by N-k.
There is no assumed independence between entry and completion.

Persistence and monotonicity leave exactly the following nonzero entries
of the actual adjacent joint table `F_k(x,y)=P(X_k=x,X_(k+1)=y)`:

    F_k(0,0) = p_0(k) - sum_d alpha_d(k) - gamma(k),
    F_k(0,d) = alpha_d(k),      F_k(0,2) = gamma(k),
    F_k(d,d) = p_d(k) - beta_d(k),
    F_k(d,2) = beta_d(k),       F_k(2,2) = p_2(k).    (2.2)

These are actual probability masses, so they are nonnegative; their row
sums are p_x(k) and their column sums are p_y(k+1). On a positive-risk row
define

    Q_k(x,y) = F_k(x,y)/p_x(k).                       (2.3)

On zero-risk rows choose, for example, a self-loop; set state 2 absorbing
in all cases. Start the new time-inhomogeneous Markov chain Y at `Y_0=0`
and use Q_k for its step k to k+1.

**Exact preservation and uniqueness.** Induction gives

    P(Y_k=x)=p_x(k),
    P(Y_k=x,Y_(k+1)=y)=p_x(k)Q_k(x,y)=F_k(x,y).

Indeed the induction step is the column-sum identity for (2.2). If
p_x(k)=0 the whole corresponding flow row is zero, so its arbitrary
transition convention contributes no mass. Any Markov chain with the
same starting law and these adjacent joint tables must use (2.3) on all
positive-risk rows. Thus its path law from this preparation is unique;
off-risk rows are not uniquely specified.

For example, the completion's directional birth kernel is

    H_d^M(a,k) = p_d(a) product_(j=a)^(k-1)
                            [1-beta_d(j)/p_d(j)],    (2.4)

when the displayed risks are positive. In all cases it is defined by
the chain itself: it is zero if p_d(a)=0 or if the surviving cohort is
exhausted. A factor at an unreachable zero-risk row is not evaluated as
0/0. The event in this kernel is `Y_a=d,Y_k=d`; because d is persistent,
it is precisely entry by a followed by survival in d to k.

This completion is an abstract marked-state process. It is not asserted
to be another physical uniform-site model with the same geometry.
Matching even all adjacent tables does **not** show that X and Y have
the same nonadjacent or full-path law. Static marked occupancies alone
do not even supply the separate entry, exit, and direct-jump flows above.

## 3. Exact characterization by completion-count cohort means

For `0<=a<=k<N` define the unconditional sector kernel and its cohort mean

    H_d(a,k) = P(J1<=a, X_k=d),
    m_d(k)   = E[nu_2(A_k) | X_k=d],
    m_(d,a)(k) = E[nu_2(A_k) | J1<=a, X_k=d].

Means are defined only on positive risks; `H_d(k,k)=p_d(k)`.

**Theorem.** The actual process X is ordinary time-inhomogeneous Markov
in its natural marked filtration if and only if, for every count `k<N`,
every direction d with `p_d(k)>0`, and every cutoff `a<=k` with
`H_d(a,k)>0`,

    m_(d,a)(k) = m_d(k).                             (C)

Every count is required, including future counts beyond a selected
measurement window. The condition is on conditional **means**, not on
pointwise constancy of nu_2, its entire conditional distribution, or the
whole microscopic geometry. Zero-risk cohorts impose no conditional test.

### Necessity

By the tower property applied to (1.1),

    P(X_(k+1)=2 | J1<=a,X_k=d) = m_(d,a)(k)/(N-k),
    P(X_(k+1)=2 | X_k=d)      = m_d(k)/(N-k).         (3.1)

The first-entry cutoff is observable from `F_k^X`. If X is Markov, its
one-step law conditioned on this additional past event cannot change,
which proves (C). The denominator N-k is common and strictly positive.

### Sufficiency

On `{X_k=d}`, the entire observed path has the form

    0,...,0, d,...,d,

with its sole entry at J1. Thus the trace of `F_k^X` on this risk set is
generated by the finite collection of events `{J1=j}`, or equivalently
the entry-CDF events `{J1<=a}`. Condition (C) is equivalently the
undivided relation

    E[1{J1<=a,X_k=d}nu_2(A_k)]
      = H_d(a,k)m_d(k).                              (3.2)

It also holds when H_d(a,k)=0, since then both sides vanish. Subtract
(3.2) at a=j and a=j-1. Every positive-mass entry atom consequently has
mean m_d(k). It follows that

    E[nu_2(A_k) | F_k^X] = m_d(k)  on {X_k=d}.

Equation (1.1), another application of the tower property, and persistence
now give the complete one-step observed law: exit to 2 has probability
`m_d(k)/(N-k)`, otherwise the next state is d. Both depend only on k,d.

On `{X_k=0}` the observed past is all zero, so its trace sigma-algebra
is trivial. Its one-step law is already the row Q_k(0,.), without any
extra criterion on the hidden geometry. On `{X_k=2}` the future is
identically 2. Therefore the one-step conditional law everywhere is
`Q_k(X_k,.)`. Repeated conditioning on successive deterministic counts
gives the usual multi-step Markov property and the path law of section 2.
This proves sufficiency, including direct double jumps and null risks.

### Kernel form and its all-count requirement

The sector hazard recursion from the companion identity is

    H_d(a,k+1)=H_d(a,k)[1-m_(d,a)(k)/(N-k)].          (3.3)

At positive risks, the adjacent triple determinant is exactly

    H_d(a,k+1)p_d(k)-H_d(a,k)H_d(k,k+1)
      = -H_d(a,k)p_d(k)[m_(d,a)(k)-m_d(k)]/(N-k).   (3.4)

Its undivided expectation form remains valid at zero cohort risk.
Thus checking (C) at all counts and entry cutoffs is equivalently
checking all these adjacent triple identities. Iteration then gives,
for **all** `a<=b<=c`,

    H_d(a,c)p_d(b)=H_d(a,b)H_d(b,c).                 (3.5)

The theorem is also the directional version of the earlier all-triples
criterion. It is not a proposal that one selected triple, one aggregated
directional contrast, or one finite-sample null result proves closure.

## 4. The immediate necessary condition versus a finite-window zero

Fix `a<b<N` and direction d. Within `X_b=d` let

    E = {J1<=a,X_b=d},       L = {a<J1<=b,X_b=d},

with both masses positive. Write m_E(b),m_L(b) for their mean nu_2.
Condition (C) implies their equality: m_d(b) is their probability-weighted
average, and the early cutoff mean must equal this average. Conversely
equality for this split supplies only this cutoff's necessary condition.
Exactly,

    P(J2>b+1 | E)-P(J2>b+1 | L)
      = [m_L(b)-m_E(b)]/(N-b).                      (4.1)

A nonzero population within-direction completion-count contrast therefore
rules out the current marked-state Markov model for this finite count
contract. No extrapolation to another size, clock, or scaling limit is
included. A sample estimate must retain its uncertainty; a zero estimate
is not an all-count identity.

For a longer window, retain the same initial E or L cohort and condition
on its survival to each j. Where those survivor risks are positive, let

    h_G(j)=E[nu_2(A_j) | G,J2>j]/(N-j),   G=E,L.

Then, without assuming Markovness,

    P(J2>c | G)=product_(j=b)^(c-1)[1-h_G(j)].        (4.2)

If a cohort is exhausted its survival is zero, with no conditional mean
assigned afterwards. The products can agree even though their first
factors, or several intermediate factors, differ. These are survivors of
the fixed birth-by-b groups, not new rank-one entrants after b. Averaging
across directions can introduce further cancellation; equality of an
overlap-weighted contrast is weaker than equality in every direction.

The parent's designated saved-geometry comparison concerns only (4.1)
at its already chosen cutoffs. It is a **post-hoc** use of the existing
140k block, not new independent evidence, a new coordinate/cutoff scan,
or a test of every condition in the theorem. This note neither reads that
block nor anticipates its result; interpretation belongs in the
[directional hazard companion](directional-hazard-contrast-20260929.md).

## 5. Exact nonphysical toy: window cancellation

Consider an abstract two-birth process with one persistent direction and
the following probability masses; all unlisted pairs have zero mass:

| (J1,J2) | Probability |
|---|---:|
| (1,3) | 1/8 |
| (1,4) | 1/8 |
| (1,5) | 1/4 |
| (2,3) | 1/4 |
| (2,5) | 1/4 |

This is hand algebra for a probability toy, **not a percolation model**
or a claimed realization of a site-completion count. At b=2 all paths
are rank one. Early `J1=1` and late `J1=2` each have probability 1/2.

At c=3 their survivals are respectively

    (1/8+1/4)/(1/2)=3/4,       (1/4)/(1/2)=1/2.

At c=4 they are both `(1/4)/(1/2)=1/2`. The step 2 to 3 completion
probabilities are 1/4 early and 1/2 late. Among survivors at count 3,
the next completion probabilities are 1/3 early and zero late. Thus

    (3/4)(2/3)=1/2=(1/2)(1),

which explicitly exhibits cancellation of the first-step difference.

Here `H(1,2)=1/2`, `H(2,2)=1`. At c=3,
`H(1,3)=3/8`, `H(2,3)=5/8`, so the determinant is
`3/8-(1/2)(5/8)=1/16`. At c=4, `H(1,4)=1/4` and
`H(2,4)=1/2`, so it is zero. The process is not Markov despite the
second triple's exact zero.

Its flow-defined Markov completion has exit probabilities 3/8 at count
2 and 1/5 at count 3. Both entry groups in that completion consequently
have survival 5/8 to count 3 and 1/2 to count 4. It matches the toy's
single-time laws and adjacent joint tables but not its entry-conditioned
survival to count 3. This is the precise distinction between constructing
a compatible Markov model and identifying the actual path law as Markov.

## 6. Clock and interpretation boundaries

- **Count versus label.** All transitions and the criterion (C) above
  use the uniform-permutation count clock. For iid Uniform site labels,
  the instantaneous completion intensity is `nu_2(A_t)/(1-t)`, not
  `nu_2/(N-k)`. The label observer hides the random occupied count;
  count-clock closure does not automatically survive this random time
  change. It requires its own sector kernels or intensity analysis.
- **Observer.** A persistent direction makes the rank-one past exactly
  `(J1,d)`. If the mark can change, its previous transitions remain in
  the history and the entry-cutoff proof is no longer sufficient. Neither
  revealing nor forgetting directions generally preserves Markovness.
- **Immediate versus recursive.** Pointwise nu_2 determines immediate
  completion hazard, not the successor distribution of nu_2. No Markov
  claim for `(D,nu_2)` follows here. Conversely, marked-state Markovness
  requires equality of the cohort means, not identical microscopic
  configurations or identical nu_2 distributions.
- **Predictive, not causal.** Selecting rank-one survivors and early/late
  entries is observational conditioning, not an intervention on birth or
  geometry. These formulas assert neither independent births nor an
  explanation of the original normalized U observable.
- **What is proved.** The completion and iff criterion are finite-model
  probability deductions. They contain no large-L or continuum claim,
  no certification of the empirical companion's outcome, and no general
  novelty assessment.
