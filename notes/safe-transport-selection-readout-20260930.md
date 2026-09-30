# Safe transport versus survivor selection: one mechanism comparison

2026-09-30. Contract and implementation committed before new continuations
at `5054c02383242473834045fcbddc7379d83acfb9`, in Draft #838.
**Completed.** The negative early-minus-late completion-creation contrast
persists after removing natural survival reweighting. For this specified
reference comparison, selection-only is disfavoured; entrance/safe transport
is needed. The selection contribution itself remains unresolved, not zero.

| Same1198 endpoint cells, same weights | Estimate | Aligned14-batch SE |
|---|---:|---:|
| Natural-alive contrast | -0.002205065286 | 0.000787091700 |
| Uniform-safe reference contrast | **-0.002537058794** | **0.000636492674** |
| Selection = natural minus reference | +0.000331993508 | 0.000287375294 |

The reference contrast is about3.99 batch SE from zero. Selection's point
estimate opposes it, but its sign is not resolved (about1.16SE); do not call
that a proved cancellation or a measured mediation percentage. These are
correlated components of one experiment, not three independent confirmations.

All140000 original J1 values match. There are97666 rank-one entrances
(56286early,41380late),42329 no-entrance prefixes and5 direct rank-two
entrances. Every reference entrance reaches b; none has an undefined safe
endpoint. The natural flag survives on54850 endpoints (24216early,30634late).
Common support contains54295natural endpoints and94123reference endpoints.
The strong cohort difference in overall survival is retained; an unresolved
selection term for this matched Y contrast does **not** make selection
irrelevant to survival generally.

Weighted natural early/late means are .09610311808/.09830818336; reference
means are .09540936877/.09794642757. The exact decomposition and singular
3x3 covariance are saved with all14 deletions and cell means in
[result.json](../analysis/safe-transport-selection-20260930/results/result.json).
This is new conditional continuation of old entrances, not a new independent
prefix block. The outcome supports moving the main question to entrance
versus safe growth under explicit sources, rather than more precision or
more descriptions of the same endpoint.

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

## If safe transport persists: an identifiable next comparison

Do not try to "hold J1 fixed and swap early/late labels": the cohorts are
defined by J1, and the birth configuration's cardinality is J1. That is not
an intervention identifying age separately from birth geometry.

Instead specify two actual source columns: an entrance-law perturbation
delta alpha_j with zero total mass at each j (preserving the first-birth
count law), and a row-sum-zero safe-transition perturbation delta P_i
(preserving safety and the insertion clock). Let P_(j:b) be the reference
safe propagator, Z=sum_j alpha_j P_(j:b)1_C, and mu its conditional endpoint
mean. Their distinct response formulas are

    entrance: Z^(-1) sum_j delta alpha_j P_(j:b)[(f-mu)1_C],
    growth:   Z^(-1) sum_(j<=i<b) alpha_j P_(j:i)
                        delta P_i P_(i+1:b)[(f-mu)1_C].

The centered terminal function includes the changing endpoint-cell
normalizer. The sources must be realizable, not arbitrary fitted columns;
these formulas define the comparison rather than supplying such a source
automatically. At least two nonredundant readouts are needed to distinguish
the two response columns; a single scalar contrast cannot identify them.
These are within-cell formulas. An overlap-weighted multi-cell response
must either fix its external cell weights in advance or include their
derivatives; it must not silently hold data-dependent weights constant.

In particular, **Y and two-step survival are a bad two-readout pair** here.
Within fixed (k,D,c), with m=N-k,s=m-c,

    survival_2 = s(s-1)/(m(m-1)) - s*Y/(m(m-1)).

Their source columns are exactly proportional after the within-cell
normalizer. More precision cannot restore the missing rank. Pooling
different c cells with different coefficients can change the apparent
weighted relation but would not create a new within-cell physical readout.

A same-contract long-window survival readout is a genuine candidate for the
second row: the already calculated L4 preparations A and C have equal c,e,Y
but different lag-three source responses. Use a real propagator rather than
replace it with two-step curvature. Whether entrance and growth columns
are independent for a specified source pair is the next mechanism question,
not answered by the current decomposition and not a request for an automatic
size/precision ladder. The earlier735-step susceptibility cannot be pasted
in as a column under a different preparation or normalization.

### Static M/E_top also cannot separate two first-birth-preserving sources

There is a second structural rank bound, complementary to the previous
[one-sided-source positive control](one-sided-birth-source-normal-response-20260930.md).
Suppose BOTH proposed sources preserve the complete first-birth marginal
and hence P0(p), and are evaluated at the same p. For column r let
g_r=partial_r P2. Probability conservation gives

    partial_r (P0,P1,P2) = g_r (0,-1,1),
    partial_r (M,E_top) = g_r (1,1).

All such source columns are proportional in these static readouts. At the
moving root their (root position, root E_top) columns remain proportional:

    partial_r (p_star,E_top(p_star))
       = g_r (-1/M_p, -2 P0_p/M_p).

Thus the genuine rank-two **thermal p plus one geometric source** Jacobian
already proved is not rank-two identification of **two geometric sources**.
More precise M/E_top data at one p cannot distinguish the latter. Multiple
p values or a joint-time/geometry readout may supply a nonproportional row,
but that must be demonstrated under the same source definition. This bound
does not apply indiscriminately to the missing original norm-4/U columns,
whose sources and observables have a different contract.

## Reproduction

See [contract](../analysis/safe-transport-selection-20260930/contract.json),
[runner](../analysis/safe-transport-selection-20260930/run.py) and
[scorer](../analysis/safe-transport-selection-20260930/score.py).
Huawei TgFr7R completed14workers in664.094seconds, ARM64/Python3.9.9/GCC10.3.1;
actual14.5CPU/25GiB limit. All14 engine exits were zero. The new task directory
reused only this task's persistent compiler in a distinct compilation subtree.
The32-old-prefix cost probe predicted430.389seconds and was excluded; the
full task took longer, with no sample extension. Output hashes and commands
are in [execution.json](../analysis/safe-transport-selection-20260930/execution.json).

One support worker reviewed source/coupling/covariance and ran one fixedL4
vacancy-query control plus one exact analytic coupling example; see
[verification.md](../analysis/safe-transport-selection-20260930/verification.md).
Main ran the fixed scorer once in research-py311 after retrieving the files.
No extra precision production, new size or repeated all-repository test suite.
