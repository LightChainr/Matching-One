# Safe transport versus survivor selection: one mechanism comparison

2026-09-30. Contract and implementation committed before new continuations
at `5054c02383242473834045fcbddc7379d83acfb9`, in Draft #838.
The result section is pending the single bounded cloud run; this is not a
claim that the experiment is complete.

## What changes the explanation

The previous natural-history contrast has two possible origins. Birth
cohorts can enter/evolve through different safe geometries, or natural
exits can select different survivors out of otherwise comparable cohorts.
Matching the endpoint direction and completion count does not match the
**path of past exit hazards**. Endpoint matching alone cannot separate them.

The [exact rejection/flag coupling](safe-path-selection-coupling-20260930.md)
now supplies an operational distinction. At count i, with m_i vacancies
and c_i completing sites, the killed natural transition on rank-one states is

    K_i(A,A+v) = (1-c_i/m_i) P_safe,i(A,A+v),
    P_safe,i(A,A+v) = 1/(m_i-c_i), for safe v.

Thus a natural safe path has relative weight

    W = product_i (1-c_i/m_i)

under the always-safe reference source. This is a pathwise killing weight,
not a proposed additional endpoint descriptor. The natural-alive flag I
has E[I | safe path]=W. For each cohort g and endpoint cell C,

    mu_natural(g,C)-mu_safe(g,C)
       = Cov_safe(Y,I | g,C)/E_safe[I | g,C].

Subtract early and late means using the same overlap weights throughout:

    Delta_natural = Delta_safe + Delta_selection.

**Selection-only for this specified comparison requires Delta_safe=0.** A
resolved safe-reference contrast requires information already in entrance
and/or safe transport. If both terms are appreciable with opposite signs,
survival selection masks part of that contrast rather than creates it.
These predictions are fixed before reading the new continuations. No fourth
descriptor or adaptive time window will be added to explain the outcome.

### What the selection term means physically

Completing vacancies are persistent along a safe rank-one path: once adding
v would create rank two, enlarging the occupied set without adding v cannot
make it non-completing. Safe insertion never occupies one of these sites.
Consequently c_(i+1)=c_i+d_i, with d_i>=0 the chosen safe site's synergy
degree. The killing weight therefore remembers **when completion sites
were created**, not just how many exist at the final endpoint.

For two paths with the same entry and terminal counts and equal terminal c,
if c_i on one path is at least as large at every intermediate count, its W
is no larger; the inequality is strict when an interior factor differs and
both weights are positive. Natural survival prefers later completion
creation under this pointwise comparison. This is an exact ordering, not a
claim that an arbitrary endpoint Y contrast must have a particular sign.
Endpoint c matching cannot remove it. No new feature is needed to test the
aggregate consequence: the coupled flag already samples the full weight.

| Fixed result pattern | Mechanism consequence for this comparison |
|---|---|
| Safe-reference contrast unresolved, selection resolved | Selection can account for the observed contrast; not proof of all-history closure |
| Safe-reference contrast resolved, selection unresolved | Entrance/safe transport is needed; quantify selection as unresolved, not exactly zero |
| Both resolved, same sign | Entrance/transport and selection reinforce |
| Both resolved, opposite signs | Selection attenuates or reverses the entrance/transport contrast |

These are interpretation branches, not significance-driven instructions to
extend the sample. Endpoint conditioning is retained in both laws; "without
selection" here means without the **natural survival reweighting**, not
without all conditioning or with an unmodified iid source.

## One endpoint, unchanged history groups

- Source: all 140000 original square-L512 permutations, grouped in their
  original 14 batches; replay to the **first** birth, not to the old risk set.
  A prefix that originally completed before b is still an eligible entrance.
- Entrance: rank-one first births J1<=b; direct rank-zero-to-two transitions
  and no-birth prefixes are counted separately. Early J1<=154646, late
  154646<J1<=155385.
- One distinct seeded uniform-safe continuation per entrance to b=155385.
  No-safe intermediate states enter a cemetery. The first proposal's
  completion status updates the natural flag; extra rejection queries do
  not advance count time.
- Endpoint C=(primitive direction, integer c), Y=2e/(262144-b-c).
  States with no safe endpoint have an undefined Y and are counted, not
  replaced by zero. Only cells containing both natural-survivor cohorts
  contribute; weight n_E*n_L/(n_E+n_L) is reused for reference means too.
- One three-component result and its singular full covariance, recomputing
  every weight and support on 14 aligned whole-batch deletions. New futures
  are conditional branches of archived entrances, **not new independent
  preparations or untouched holdout**.

The reference is explicitly modified dynamics. Its surviving branch has the
natural killed-path law, but the reference contrast is not a unique causal
split between imported birth geometry and subsequent safe growth. This
finite count-clock calculation does not identify original U or continuum
fields. The two earlier original-block contrasts are context, not numbers
to subtract from this new realization or pool as additional votes.

## Reproduction

See [contract](../analysis/safe-transport-selection-20260930/contract.json),
[runner](../analysis/safe-transport-selection-20260930/run.py) and
[scorer](../analysis/safe-transport-selection-20260930/score.py).
Medium computation is assigned to Huawei TgFr7R in a new task directory;
one support worker reviews the coupling and a tiny physical query control.
No precision production, new size or repeated all-repository test suite.
