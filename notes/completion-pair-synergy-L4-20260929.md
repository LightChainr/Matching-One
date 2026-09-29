# Completion count does not close the L4 process: the two-site transmission

2026-09-29. Exact finite calculation under the uniform-permutation **count
clock**, on square NN and triangular (+diagonal) occupied-site tori.

The later [independent L512 safe-insertion block](safe-insertion-independent-block-20260929.md)
has now measured the proposed mean successor contrast without resolving a
history difference. That empirical outcome does not alter the exact L4
counterexamples below or promote them to large-size conclusions.

**Result.** Current rank, persistent direction and completion count do not
give a Markov state even from the uniform empty preparation at L4. There
are positive-probability early/late first-birth groups with identical
current `(k,D,nu_2)` and identical next-step exit probability, but different
next completion-count distributions and two-step rank survival. This is
stronger than finding two specially prepared configurations with different
futures. It is not a large-size or continuum non-Markov theorem.

The microscopic transmission is explicit: two empty sites that are each
safe alone can jointly complete the second homology direction. Their
pair count, rather than another fitted descriptor, gives the exact next
term after the already known one-site completion hazard.

## 1. A fixed two-step law

For rank-one A with k occupied sites, put `m=N-k`, `c=nu_2(A)`.
Let C be its completing vacancies and S its safe vacancies. On S draw an
edge between v and w exactly when `r(A union {v,w})=2`; let e be the number
of these unordered edges. The [general derivation](completion-pair-synergy-20260929.md)
gives, for m>=2 in the two-step formula,

    nu_2(A union {v}) = c + degree(v),  v in S,

    P(rank remains 1 after two insertions | A)
      = [(m-c)(m-c-1)-2e] / [m(m-1)].

Existing completion sites cannot cease completing after a safe insertion,
by monotonicity. New completing sites are exactly the neighbours of the
inserted safe site in this graph. The graph degree distribution therefore
gives the full one-step successor law of the completion count.

At fixed k,D,c, averaging the last formula over early and late histories
produces the **exact** identity

    Delta(survival for two insertions)
      = -2 Delta(mean e) / [m(m-1)].

This is an actual source-to-transition map, not a fitted association. It
does not say that adding e closes all subsequent steps. Even the whole
current degree distribution need not specify how the graph evolves.

## 2. What was computed

One local run enumerated the `2^16=65,536` occupied subsets of each L4
lattice. The existing lifted-graph rank/direction routine was reused with
L=4; no old L3 computation was rerun. For each rank-one subset, direct
one- and two-site extensions give C, the synergy graph and successor
counts. Integer dynamic programming counts every possible J1 for its
occupied prefixes; it does **not** enumerate 16! complete permutations.

At a fixed k every k-site set is equally likely, and every ordering of
its sites has the same probability. Consequently prefix-order counts are
the exact appropriate weights for the conditional early/late comparison;
the common remaining `(N-k)!` factor cancels. All fractions below are
exact, not Monte Carlo estimates.

| Lattice | Rank-one subsets | `(k,D,c)` cells | Cells with differing e / successor laws |
|---|---:|---:|---:|
| square | 19,932 | 72 | 34 / 34 |
| triangular | 23,502 | 90 | 48 / 48 |

The saved result includes all per-cell integer birth masses, edge sums,
two-step-survival numerators and successor-count totals. Witnesses are
selected lexicographically to demonstrate existence, not ranked by effect
size or significance. This is an exact finite structural search, not a
prespecified stochastic experiment.

## 3. Failure under the original uniform preparation, not only strong lumping

Within the following current states compare `early={J1<=a}` with
`late={a<J1<=k}`. Both groups have positive probability.

| Quantity | Square | Triangular |
|---|---:|---:|
| k, a, D, c | 6, 4, (0,1), 0 | 5, 4, (0,1), 0 |
| Early / late prefix mass | 10,368 / 156,672 | 1,152 / 8,448 |
| Early mean e | 70/27 | 4 |
| Late mean e | 130/51 | 54/11 |
| Both one-step rank-one survival | 1 | 1 |
| Early mean next completion count | 14/27 | 8/11 |
| Late mean next completion count | 26/51 | 108/121 |
| Early two-step survival | 229/243 | 51/55 |
| Late two-step survival | 433/459 | 551/605 |
| Early minus late two-step survival | **−4/4131** | **2/121** |

Here c=0, so the next step always survives; conditioning on its survival
does not alter these particular successor means. Their differences are
already a **one-step** obstruction to Markovness of the augmented
`(rank,D,nu_2)` process. They become a **two-step** obstruction if the
future readout is rank alone.

J1 is observable in the augmented process history. If its future law were
determined by the current augmented state, conditioning further on that
birth event could not change the future law. Thus these are weak-Markov
counterexamples from the uniform empty start, not merely failures of
strong lumpability under arbitrary microscopic preparations.

The square and triangular signs need not agree. These finite geometries
are not a matched-modulus or universality experiment.

## 4. A small physical witness that can be checked directly

For square L4 index sites by `v=x+4y`, with x,y=0,...,3 and periodic edges.
Two six-site preparations, shown in increasing y, are

```text
A = {0,1,4,5,8,12}    B = {0,2,4,5,8,12}
##..                  #.#.
##..                  ##..
#...                  #...
#...                  #...
```

Both have vertical direction `(0,1)`, no horizontal winding and zero
one-site completions. For A the only safe pairs completing another
direction are `{2,3}` and `{6,7}`. For B they are `{1,3}`, `{3,6}` and
`{6,7}`. In B the middle pair completes the bent winding route

    0 -> 4 -> 5 -> 6 -> 2 -> 3 -> 0 across the x boundary.

The other pairs close the two partially occupied horizontal rows. A
two-site route cannot close through the two other rows, which each lack
three sites. This supplies a direct geometric check of one census
witness, in addition to the general pair-count identity.

There are 45 unordered vacancy pairs. Exactly 2 or 3 cause completion,
giving survival `43/45` versus `42/45=14/15`. The respective next-nu counts
are `{0:6,1:4}` and `{0:6,1:2,2:2}`, over 10 equally likely insertions.
This particular pair proves strong-preparation failure. The distinct
cohort calculation in section 3 is what establishes uniform-start failure.

## 5. Consequence for the research queue

The earlier L3 four-state completion-phase realization remains correct
at L3. The L4 result prevents treating it as a universal closure principle.
The [L512 conditional finite-window result](completion-hazard-independent-block-20260929.md)
also remains unchanged: it did not resolve a same-(D,nu) residual.
Neither a small-size exact counterexample nor a large-size unresolved
contrast settles the strength of this channel in a scaling limit.

The next useful large-size question is now **how much birth-history
dependence enters through completion creation**, not whether identical
current completion counts force identical immediate exit hazards (they
do by definition). One fixed prospective target is the conditional safe-step
increment of nu, which directly estimates `2 E[e]/(m-c)` at matched k,D,c.
Its corresponding two-step rank effect has the displayed known factor.

A single uniformly chosen safe insertion, followed by a new completion
count, is an unbiased probe of this drift; it does not require enumerating
all O(m^2) pairs. Independent prefixes remain the sampling units. Multiple
continuations of one prefix share its dependency group. At large N the
two-step rank difference can be tiny, so a meaningful acquisition should
compare this mechanistic drift or a specified finite-window consequence,
not blindly count rare exits or request the next full subset census.

Existing L512 files save nu at b and an independent random later count,
not a one-safe-step successor. They cannot be relabelled as this new
measurement. No such large-size production was launched in this calculation.

## 6. Reproduction

```bash
python3 analysis/completion-pair-synergy-20260929/analyze.py
```

[Result JSON](../analysis/completion-pair-synergy-20260929/result.json)
and [readable output](../analysis/completion-pair-synergy-20260929/RESULT.md).
The one standard-library Python 3.11 run took 0.540 s square and 0.693 s
triangular. The new enumeration checks the graph successor formula and
direct two-insertion count on its rank-one rows; the displayed square pair
is a separate hand-geometric check. No full suite, 9!/16! rerun, package
installation, sampling or cloud operation.
