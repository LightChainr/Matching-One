# Direction-only closure meets a sharper next-insertion prediction

2026-09-29. Same completed 140k square-L512 block; **one retrospective
necessity check**, not a new experiment. The
[original prospective readout](completion-hazard-independent-block-20260929.md)
and every original result, cutoff and contract remain unchanged.

**Independent follow-up:** the [new safe-insertion block](safe-insertion-independent-block-20260929.md)
predeclared and repeated the same-D completion-count contrast, obtaining
`-0.03987 +/- .56988`, not reproducing the `+1.62917 +/- .39098` below.
The original post-hoc result is preserved, but must not be read as a stable
empirical exclusion. L512 direction-only closure remains unresolved.

**Outcome.** Within the same primitive direction, early histories have
`1.6292 +/- .3910` more completion sites than late histories, using the
existing exact-D overlap weights and aligned batch SE. Direction-only
Markovness requires that population contrast to vanish. This is a concrete
challenge to that candidate, despite its unresolved longer-lag survival
contrast. It is not an independently validated rejection or a continuum
non-Markov theorem.

## 1. The necessary prediction comes from the sampling law

At fixed insertion count b, current rank one, and persistent direction d,
the uniform-source identity gives

    P(J2=b+1 | entry cohort g, D=d, J2>b)
      = E[nu_2(A_b) | g,D=d,J2>b] / (N-b).

If current `(rank,D)` is a time-inhomogeneous Markov state in its observed
history, the right side cannot depend on the first-birth cohort g. Thus

    E[nu_b | early,d] = E[nu_b | late,d]

is an exact necessary condition. It is more informative than simply
requesting another fitted transition curve. Averaging zero differences
with positive common-support weights must still give zero.

The [companion marked-model note](directional-markov-completion-20260929.md)
supplies the full-cutoff criterion, the flow-matching Markov surrogate,
and an exact example showing why a zero long-lag contrast cannot replace
this next-step requirement. That surrogate preserves single-time marked
laws and adjacent transition tables but need not preserve true history.

No claim here depends on choosing a new direction subset, cutoff, size,
binning scheme, regression feature or fitted exponent. The question was
selected **after** the original block had been read, including its
direction-resolved means, so its uncertainty must not be presented as a
prespecified independent confirmation.

## 2. Compute one new contrast, keep the original comparison fixed

Use the same early/late cohorts, all exact primitive unoriented directions,
and weights `w_d=e_d*l_d/(e_d+l_d)` as the original exact-D survival
contrast. Replace only the outcome by the already saved nu_b. Recompute
those weights after deleting each whole batch. All 140,000 records remain
in the same source block; rank-one cohort sizes are 24,218 and 30,505.
Shared directions retain all early histories and 30,490 late histories.

| Quantity, early minus late | Estimate | Batch SE |
|---|---:|---:|
| Pooled mean completion count | -2.18473 sites | .47950 sites |
| Within-D weighted mean completion count | +1.62917 sites | .39098 sites |
| Within-D **next-insertion exit probability** | +0.0000152602 | .0000036622 |
| Same exit contrast multiplied by c-b=735 | +0.0112163 | .0026917 |
| Original within-D survival through c | -0.0041500 | .0045119 |

The fourth row is a scaled intensity contrast, **not** a 735-insertion
exit probability. It is included solely to put the microscopic one-step
quantity on the existing window's scale. The old survival estimate and
SE are reproduced by the same weights, without changing its definition.

The direction-adjusted completion contrast is about 4.17 batch SE from
zero. Fourteen-batch uncertainty is approximate, this comparison is
retrospective, and no multiplicity-adjusted certification is claimed.
Its positive sign means greater immediate early-cohort exit propensity.
The longer-lag survival estimate has the consistent negative sign but
substantially less separation from zero. The data do **not** demonstrate
temporal cancellation; lower precision or changing hazard differences
can both leave the finite-lag comparison unresolved.

The pooled completion difference has the opposite sign. This is consistent
with the already observed direction-composition difference between cohorts.
It does not license a causal or mediated fraction: the pooled and
within-direction comparisons target differently weighted populations.
The narrow statement is that direction composition and within-direction
completion geometry need not point in the same direction.

## 3. What to retain, revise, and ask next

- Retain the original result: the prespecified same-D and same-(D,nu_b)
  finite-lag survival contrasts were unresolved. Do not rewrite them as
  significant, or discard the experiment because a primary outcome is weak.
- Revise the working model: a persistent-direction **Markov** mixture is
  now a challenged candidate, not a preferred closure justified by one
  unresolved window. Direction remains scientifically useful; only its
  purported sufficiency is under question.
- Keep `(D,nu_b)` predictive sufficiency unresolved. Its immediate exit
  probability is fixed by construction, so equality of nu_b inside an
  exact-(D,nu_b) cell is a tautology, not a successful test of recursive
  dynamics. Its longer-lag original comparison remains the relevant datum.

The next discrimination should concern **evolution of completion geometry**,
not another list of static descriptors. A concrete future experiment could
compare conditional continuations from matched `(D,nu_b)` preparations,
with first-birth history as the contrasting label, and estimate a fixed
future survival or completion-count transition. Multiple continuations of
one prefix remain one prefix-level dependency group; they are not extra
independent configurations. Define that prediction and its cost before any
new acquisition. No such continuation campaign is started by this note.

An independent confirmation of the present same-D immediate-hazard
contrast would likewise have one fixed direction weighting, the existing
cutoffs, and a declared sample budget, not a search for a more favourable
time or sign. The present analysis alone does not determine a scaling
limit, an original-U operator, or a minimum hidden-state dimension.

## 4. Reproduction and evidence scope

```bash
python3 analysis/directional-hazard-contrast-20260929/analyze.py
```

The new small scorer reuses the unchanged production reader and batch
sufficient-statistic code. It reads each compressed batch once, verifies
the original hashes, and saves the three-dimensional joint covariance and
all aligned delete-one vectors in
[result.json](../analysis/directional-hazard-contrast-20260929/result.json).
[RESULT.md](../analysis/directional-hazard-contrast-20260929/RESULT.md)
provides the numeric view. Both input manifest and reader hashes are saved.

The source is the `87f6644d` block, not a second random block. One local
standard-library Python 3.11 analysis took under one second; no packages,
cloud access, simulation, source regeneration, cutoff scan or full test
suite was needed. A narrow arithmetic identity confirms that the retained
old within-D survival point/SE are unchanged; this is not an independent
validation of the entire acquisition pipeline.
