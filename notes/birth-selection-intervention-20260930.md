# A clock-preserving geometric source: from history correlation to intervention

2026-09-30. Author finite-model derivation in Draft #838. The main research
advance is a **source with zero immediate rank response but a specified
delayed response**, not a more precise estimate of the old contrast.
One separate worker handles the old target's independent-archive replication.

**Outcome.** We can change which safe site arrives while preserving both the
entire occupied-count clock and the current completion hazard. This changes
the two-step rank law with an exact, nonpositive susceptibility. Birth-history
association, meanwhile, factors into a backward homology-retention functional
and a forward completion functional. They are not the same geometric object.
An exact path-weight decomposition separates safe geometric transport from
survivor selection. These give a concrete next mechanism experiment.

All statements below are for the specified finite occupied-site torus and
uniform source at theta=0. They do not identify original U, a conformal
operator, a continuum memory coefficient or the sign of a long-window effect.

## 1. The intervention, not another passive descriptor

Fix a rank-one occupied set A with k sites, m=N-k vacancies, c completing
vacancies and s=m-c safe vacancies. In its existing synergy graph, let
d_v be the degree of safe vacancy v and e the edge count. Then

    mean_safe(d) = 2e/s,
    c(A+v) = c+d_v.

For s>0, define a new insertion source:

    P_theta(v | A) = 1/m                         if v completes,
                    (s/m) exp(theta*d_v)/Z_A(theta)  if v is safe,
    Z_A(theta) = sum_safe exp(theta*d_v).                    (1)

Use the original uniform source outside rank one or when s=0. At theta=0
this is exactly the original source, not a sampled approximation. For every
theta the total probability of completing on the next insertion is c/m.
The rank-one direction also stays fixed on every safe insertion. Thus this
source changes continuation geometry while leaving the immediate coarse
rank response **identically zero**, configuration by configuration.

Let S2(theta|A) be rank-one survival after two insertions, with the kick (1)
on the first insertion and the ordinary uniform source on the second. Then

    S2(theta|A) = (s/m) [1 - (c+E_theta[d])/(m-1)],
    d S2/d theta |0 = -s Var_safe(d) / [m(m-1)].             (2)

Proof: after the first safe choice v, exactly c+d_v vacancies complete;
the second insertion exits with probability (c+d_v)/(m-1). Differentiate
the normalized exponential mean. Its derivative is Var_safe(d).
The result also holds when (1) is active on both insertions: the second
step's total exit mass is still unchanged. Assume m>=2; s=0 has zero
survival and zero source response, with no need to define a degree variance.

Consequently the two-step response is strictly negative whenever the safe
degrees vary. It is zero for a regular synergy graph, even if e is large.
The number e controls the baseline two-step law; degree variance controls
this *specified perturbation's* gain. This is a consumer for graph structure,
not a proposal to rank arbitrary extra features by explained variance.

### Continuous clock: the occupation process itself is unchanged

In the iid-label baseline at time t, every remaining site has arrival rate
1/(1-t). Replace safe-site rates by

    r_theta(v|A,t) = [s/(1-t)] exp(theta*d_v)/Z_A(theta),

and leave completing-site rates equal to 1/(1-t). Total arrival rate stays
m/(1-t) for every A. The entire occupied-count process K(t) therefore has
the same law for all theta, and can be coupled with identical arrival times.
The current rank-exit intensity c/(1-t) is also unchanged. At nonzero theta
these are **modified geometry-dependent arrivals, not iid site labels**.

The second survival derivative at a fixed starting geometry is

    S_theta''(0|A) = [c(c-1)-s E_theta[d]]/(1-t)^2,
    partial_theta S_theta''(0|A)|0
                  = -s Var_safe(d)/(1-t)^2.                 (3)

The killed-generator calculation replaces the old 2e creation term by
s E_theta[d]; its loss term is still c^2, and differentiating 1/(1-t)
supplies the -c term. Hence the geometric source first reaches rank at
order h^2, while changing neither K(t) nor the first-order exit rate.
This distinguishes it from a mere change of the arrival clock. It does not
by itself identify a CFT normal direction or rule out thermal projections
in a different readout contract. Remainders here are finite-N, not uniform
near-critical bounds.

## 2. Birth memory is a backward/forward coupling, not an ageing force

For a rank-one k-site set A define its backward retention probability

    R_a(A) = #{B subset A: |B|=a, rank(B)=1} / binom(k,a).

Under the original uniform permutation, conditional on A_k=A the first a
sites form a uniform a-subset. Monotonicity ensures no rank-two subset is
possible. Therefore

    P(J1<=a | A_k=A) = R_a(A).                              (4)

The remaining permutation is independent of the internal ordering of A.
Let C be any positive-probability current-state cell, for example exact
(rank=1,D,c), and p_C=E[R_a|C] in (0,1). For any future readout Z, write
z(A)=E[Z|A_k=A]. Conditional past/future independence at full A gives

    E[Z|J1<=a,C] - E[Z|J1>a,C]
        = Cov_C(R_a,z) / [p_C(1-p_C)].                      (5)

For two-step survival, all c-dependent terms are constant within C, so

    Delta S2 = -2 Cov_C(R_a,e) / [m(m-1)p_C(1-p_C)],
    Delta E[Y] = 2 Cov_C(R_a,e) / [s p_C(1-p_C)],
    Y(A)=2e/s.                                             (6)

The current L512 matched negative Y contrast thus estimates a negative
overlap-weighted sum of these normalized covariances, not a universal
negative sign in every cell. More first-direction robustness under thinning
is associated, in that weighted readout, with less second-direction creation.
This interpretation is exact for the *population target*; finite estimates
and replication retain their own uncertainty. Merely finding e variation,
or saying that older paths are different, does not establish this coupling.

In label time, conditional on A_t=A, retained earlier sites are independent
with probability u=s0/t. Set R_A(u)=P(rank(thin_u A)=1). Equation (4) becomes
P(T1<=s0|A_t=A)=R_A(s0/t), and the analogous covariance identity holds under
the correct label-time law of A_t. This is not permission to substitute
t=k/N into a fixed-count archive.

### Three specified physical controls

On square L4, v=x+4y, take

    A={0,1,4,5,8,12}, B={0,2,4,5,8,12}, C={0,1,5,8,9,12}.

All have k=6, D=(0,1), c=0. A/B are the existing physical pair; C is a
specified six-site bent winding loop. Only these three preparations and
their subsets/extensions are evaluated, not a new whole-lattice census.

| Preparation | e | R_A(u) | Var_safe(d) | S2(0) | S2'(0) | S2(log 2) |
|---|---:|---:|---:|---:|---:|---:|
| A | 2 | u^4 | 6/25 | 43/45 | -2/75 | 59/63 |
| B | 3 | u^4 | 16/25 | 14/15 | -16/225 | 71/81 |
| C | 2 | u^6 | 6/25 | 43/45 | -2/75 | 59/63 |

A and B have the **same entire conditional first-birth law** but different
future and intervention laws. A and C have the same c,e and this two-step
response but different birth laws. Thus backward robustness and forward
completion are distinct axes; neither can generally be reconstructed from
the other. A 50:50 A/B preparation has zero birth-history association at any
nondegenerate cutoff despite heterogeneous futures, because its two R_a
values coincide. These are specified-preparation controls, not population
estimates for natural L512 prefixes.

The longer-window response has also been **calculated**, not left as a
request to extend the local formula. Exact continuation over subsets
reachable from these starts gives single-kick survival derivatives:

| Further insertions | A | B | C |
|---|---:|---:|---:|
| 1 | 0 | 0 | 0 |
| 2 | -2/75 | -16/225 | -2/75 |
| 3 | -7/180 | -173/1800 | -17/450 |
| 4 | -17/525 | -307/4200 | -6/175 |
| 5 | -2/105 | -47/1260 | -2/105 |
| 6 | -1/150 | -67/6300 | -2/315 |
| 7 through 10 | 0 | 0 | 0 |

Mean remaining steps-to-completion derivatives are -779/6300, -454/1575
and -391/3150 respectively. A/C first diverge at lag three despite equal
initial c,e and degree variance. The response peaks in magnitude at lag
three in these controls and vanishes when all continuations have completed;
it is not a constant curvature accumulated over time. All nonzero responses
here are negative; these three controls do not prove that sign for arbitrary
geometries or windows. No whole-lattice census was rerun.

## 3. Separating geometric transport from survivor selection

The previous covariance identifies what is coupled now, not how that
coupling developed. There is an exact decomposition with a stated reference
process; it is not a unique causal attribution from endpoint observations.

Start with a rank-one entrance distribution nu_j at count j. Let P^S_i
choose uniformly among safe vacancies at count i, keeping one trajectory
per starting prefix instead of discarding those with large c. Along this
reference safe path put

    W_{j:k} = product_{i=j}^{k-1} (1-c(A_i)/(N-i)).

The natural killed transition is Q_i(A,A+v)=1/(N-i) for safe v, hence
Q_i=(1-c/(N-i))P^S_i. Multiplying along a path proves

    E_natural[1{rank_k=1} f(A_k)] = E_safe[W_{j:k} f(A_k)],
    E_natural[f | rank_k=1]
        = E_safe[f] + Cov_safe(f,W)/E_safe[W].               (7)

If a reference path has no safe vacancy before k, send it to a cemetery
and set its subsequent W to zero. The same identity holds. To condition
on an endpoint cell C, first condition the reference law on A_k in C and
use its conditional covariance and weight mean; exclude the cemetery.

For actual birth cohorts let alpha_j(A) be the unnormalized original
probability of first entering rank one at count j in configuration A.
Direct rank0-to-rank2 jumps are not in this flux. Then, for example,

    H(a,k) = sum_{j<=a} sum_A alpha_j(A) E_safe,A[W_{j:k}].  (8)

Late cohorts use a<j<=k. Thus the first-birth ensemble, safe transport and
selection weights enter separately. Apply (7) to each normalized cohort's
mixture over j to decompose its endpoint f, and subtract the two identities.
Different cohort weights and endpoint-cell normalizers must be kept; a
single uniform average of safe paths is not the original survivor law.

This is the precise reason that a large current e-history association
cannot alone decide “imported at birth” versus “selected during survival.”
Post-hoc relabelling ages on fixed geometries is not a physical intervention.
Nor is this a request for a larger descriptor catalogue: the reference
process in (7) directly switches off survival reweighting while retaining
explicit within-rank-one geometric evolution.

## 4. The next decision, with distinct predictions

The main lane is now **mechanism intervention and propagation**, not further
precision rounds. Independent archival replication runs on a support lane.

1. **Geometric continuation kick.** Use (1) on the same starting prefixes,
   sharing the arrival clock with the unperturbed branch. One-step rank
   response must be zero; two-step susceptibility is (2), and for the three
   controls the finite responses are already calculated above. For a
   specified longer window, predict the response by the actual continuation
   law rather than extending the two-step sign. A model that represents the
   perturbation only as a clock change misses this zero-first/nonzero-delayed
   response. The useful large-size question is whether this channel has a
   nonvanishing *scaled finite-window gain*, not another confirmation of (2).
2. **Selection versus transport.** For specified first-birth ensembles and
   a common destination count, compare the unweighted safe reference to
   its W-weighted counterpart in (7). “Selection-only for the chosen
   contrast” predicts zero unweighted contrast and a nonzero weight
   covariance contribution. A nonzero unweighted contrast excludes that
   explanation: entrance/transport already carries it. Both nonzero means
   a mixed mechanism; cancellation is also allowed. This decomposition
   does not by itself split entrance geometry from subsequent safe evolution.
3. **Source at entrance, if that split is needed.** At fixed j,D,c, tilt an
   explicitly specified entrance ensemble by exp(theta*g(A)), then use the
   original continuation. Its two-step derivative is
   -2 Cov_nu(g,e)/[m(m-1)]. Choosing g=e gives a definite negative result
   when e varies (equal A/B: -1/180). This is a separate *preparation*
   source, unlike (1)'s *transition* source. In a longer-window comparison,
   keep j and source normalization fixed, not an age swap across sizes.

For an instantaneous transition kick at i, the finite-window derivative is

    nu Q_k ... Q_{i-1} B_i Q_{i+1} ... Q_{l-1} 1,
    B_i(A,A+v) = [d_v-2e/s]/(N-i)  for safe v.              (9)

Other entries of B_i are zero; its row sum is zero. Summing (9) gives the
response to a field active over several counts. This is the explicit
source -> geometric transition -> future-rank map to evaluate, not a
generic request for a new certificate. A long-window sign and a continuum
coefficient are unresolved predictions, not consequences of local negativity.

No new large production is launched by this note. First fix which of these
mechanisms the longer-window calculation distinguishes; medium/large runs
go to Huawei. These are attention priorities, not locks on other research.

## Artifacts

[Exact controls](../analysis/birth-selection-mechanism-20260930/RESULT.md),
[fractions and subset counts](../analysis/birth-selection-mechanism-20260930/result.json),
[short calculation](../analysis/birth-selection-mechanism-20260930/derive.py).
The script reuses only the old lifted topology routine; no old census,
Monte Carlo, full test suite or numerical differentiation is repeated.
