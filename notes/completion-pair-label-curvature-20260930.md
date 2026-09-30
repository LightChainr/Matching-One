# Completion pairs determine label-time survival curvature

2026-09-30. Finite-model author derivation, connecting the
[completion-pair identities](completion-pair-synergy-20260929.md) to the
[label-time kernel](completion-pivotal-hazard-kernel-20260929.md).
No new production, fitted exponent, literature-novelty claim or continuum
limit assertion. This is a consumer for the same geometric quantity e,
not an additional independent experimental outcome.

## 1. Conditional short-time survival

At label time t<1, fix the entire revealed past with current rank-one set A.
The m remaining site labels are independent Uniform(t,1). Let c be the number
of individually completing vacancies and e the number of safe synergy pairs.
During (t,t+h], each vacancy arrives independently with probability
q=h/(1-t). For this fixed finite graph,

    S_A(h) = P(T2>t+h | revealed past)
           = 1 - c*q + [c*(c-1)/2 - e]*q^2 + O_N(q^3).       (1)

Proof: survival requires that no completing singleton arrives, contributing
(1-q)^c. Among safe sites, each completing unordered pair contributes q^2
to failure at second order. Intersections of distinct pairs involve at least
three sites, as do completing sets with no completing singleton or pair.
Their contribution starts at degree three. Thus the remaining factor is
1-e*q^2+O_N(q^3), proving (1) by multiplication. This is a finite polynomial
expansion; its remainder bound need not be uniform as N grows.

In particular,

    S_A'(0)  = -c/(1-t),
    S_A''(0) = [c*(c-1)-2e]/(1-t)^2.                         (2)

An explicit hand control is one completing singleton plus two safe sites
whose pair completes: survival is exactly (1-q)*(1-q^2). Its second-degree
coefficient is -1, as (1) gives for c=e=1. This is a monotone finite control,
not a claimed torus witness or a new census.

For two revealed-past cohorts with the same t and exact c (and any matched
persistent direction D), their initial hazards coincide. Their curvature
difference is nevertheless

    Delta S''(0) = -2*Delta E[e]/(1-t)^2.                    (3)

Condition on fixed initial cohorts, then let them evolve; do not continually
reselect a fixed-c cohort as time advances. A difference in (3) is a local
history-dependent survival response despite identical initial hazard.

## 2. The fixed-entry kernel, without conditioning on c

Fix s<t and write I_t=1{T1<=s,T2>t}, H(s,t)=E[I_t]. No new member enters
this cohort as t advances. On positive risk, let

    mu_s(t)=E[c(A_t)|I_t=1],   eps_s(t)=E[e(A_t)|I_t=1],
    v_s(t)=Var(c(A_t)|I_t=1).

Set the killed count to zero after completion. Summing over all m potential
insertions, each with rate 1/(1-t), gives

    d/dt E[I_t*c(A_t)] = E[I_t*(2e-c^2)]/(1-t).              (4)

Indeed a completing site loses c, giving -c^2 in total; safe sites add their
synergy degrees, summing to 2e. This is the existing killed count drift with
the actual label clock, not with k replaced by Nt. Since
partial_t H=-E[I_t*c]/(1-t), differentiation yields

    partial_t^2 H(s,t) = E[I_t*(c*(c-1)-2e)]/(1-t)^2,       (5)
    mu_s'(t) = [2*eps_s(t)-v_s(t)]/(1-t).                   (6)

The hazard lambda_s=-partial_t log H=mu_s/(1-t) therefore obeys

    lambda_s'(t) = [mu_s + 2*eps_s - v_s]/(1-t)^2.          (7)

The variance term is selection among surviving geometries: high-c
configurations leave sooner. It must not be omitted. Safe completion
creation is nonnegative at each fixed geometry, while the risk-set mean c
can decrease because v_s can exceed 2*eps_s. Neither (6) nor (7) requires
rank-only Markovness or a recursively closed pair graph. Fixed disjoint
early/late entry cohorts satisfy the same formulas after their common
starting time. Direct rank0-to-rank2 births do not enter these rank-one cohorts.

There is also an exact count-clock counterpart. For a fixed rank-one cohort
at count k, put m=N-k, mu=E[c], eps=E[e], v=Var(c). If m-mu>0, its mean
completion count after the next insertion, conditional on surviving, obeys

    mu_next - mu = (2*eps-v)/(m-mu).                          (7a)

The surviving numerator is E[(m-c)*c+2e], and its normalizer is m-mu;
subtracting mu gives (7a). The negative variance term is again selection.
Forcing one safe insertion in every original prefix would instead average
c+2e/(m-c) without survival reweighting. These agree on an exact-c cell,
not generally on a mixed-c cohort. This distinction is why the current
exact-(D,c) target is well-defined but is not a whole-cohort hazard derivative.

## 3. What this changes, and what it does not

This provides the next theoretical map in explicit coefficients:

    completion count -> instantaneous exit;
    pair creation minus survival selection -> hazard change;
    their cohort contrasts -> history-dependent kernel curvature.

At a deterministic near-critical rescaling t=t0+a_L*tau, a candidate limiting
kernel must account for the scaled quantity

    d lambda_scaled/d tau
      = a_L^2*[mu_s+2*eps_s-v_s]/(1-t)^2.                   (8)

Equation (8) is just a finite-chain-rule identity, with the entry cutoff held
fixed when differentiating. Passing it to a limit needs moment control and
justified differentiation/convergence. A small fixed-N Taylor remainder
does not control a near-critical interval containing many microscopic arrivals.
The present notes supply no limit of eps_s or v_s, no cross-lattice exponent,
and no nonzero continuum memory coefficient.

The current L512 replay fixes an **insertion count**, not a label time. It
exactly scores its original two-insertion/safe-step target. It cannot be
substituted into (5)--(8) by taking t=b/N. A label-time numerical consumer
needs the correct count/history mixture under iid labels, or an explicitly
matched label-time acquisition. Nor may a two-step contrast be extrapolated
to the old 735-insertion endpoint by assuming constant curvature.

The useful next theory question is thus quantitative: do competing continuum
descriptions predict the same or different scaled pair/selection combination
in (8), under the same observer and entry-history contract? Merely appending
another geometric feature does not answer this question. This is attention
guidance, not a restriction on parallel research.
