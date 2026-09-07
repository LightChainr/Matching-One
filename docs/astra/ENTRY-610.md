# Entry note for #610 — orientation, not instructions

**This file changes the mode of the pack.** Every other file in `docs/astra/` is written for
a model that **cannot** be shown the repository: self-contained, one file, no context.
[#610](https://github.com/LightChainr/Matching-One/issues/610) assumes you **can** read the
tree, and that you should read whatever you judge to be worth reading.

Everything below is *information*, not a route. There is no required reading order, no
budget, no file you are supposed to stay out of, and no format your answer has to take. The
one thing this file is for is keeping you from wasting effort on things that are already
settled or already known to be stale. If any of it turns out to be wrong, say so — that is
a result too.

---

## The one navigational fact

**Read `claude/matching-one-workspace-pwr5pv`.** `main` is several days behind and none of
the four gates #610 is about have landed there. Two other branches carry work that has not
been merged yet:

- `claude/p581-empirical` — Gate 4, the typed `Q`-score control on square bond.
- `docs/literature-officer-20260906-issue601` — the literature answers to #601.

Anything else you find on `main` is real but old. Open pull requests are the best index of
what is in flight.

---

## What the vocabulary means here

Worth two minutes because it changes how the evidence reads.

`GOVERNANCE.md` §2 is the whole rule set while exploring — five items. Two of them shape
almost every artifact you will open: *don't fool yourself about a number, check once by
independent means*, and *count one random block once*. When something in this repository
says a check "fired," it means an independent check returned a contradiction and it was
reported rather than absorbed.

`docs/STATUS.md` carries claim levels. `C3` is a scored prospective test, `C2` is
reanalysis of committed productions. A row can read **negative** and still be the strongest
thing on its line — a large fraction of them do. That is the house style, not modesty, and
a negative answer from you is treated the same way.

Scripts in this repository are written to be read: the module docstring usually states what
the object is and what the wrong answer would look like, and test docstrings name the
specific wrong number the test exists to stop us believing. If you want to know what
someone thought they were measuring, the docstring is normally more honest than the note.

---

## Where things are

Not a reading list — a map, so you can go straight to whatever you actually want.

**The five #582 transitions and everything built on them** (this is what Q1 is about):

```text
scripts/threshold_quantile_lineage.py       quantile reconstruction, jackknife,
                                            spin-0 orientation weighting
scripts/score_wasserstein_shape_flow.py     the affine projection and the shape flow
scripts/score_type582_residual.py           Gate 3: the label screen, exact permutation null
scripts/p582_amplitude_law.py               the exponent, the curvature check, the four-cell
                                            weighting control
notes/p582-amplitude-law-20260906.md        the result in prose, including a correction to
                                            my own first reading of it
notes/type582-residual-20260906.md          Gate 3 in prose
results/wasserstein-shape-flow/latest.json  the raw material: per-size quantiles, standard
                                            errors, covariance inputs, orientation data
results/p582-amplitude-law/latest.json      every number quoted in #610's Facts 3 and 4
results/type582-residual/latest.json        the label screen's full output
```

`results/wasserstein-shape-flow/latest.json` is worth knowing about specifically: nine
quantiles and their standard errors at eight sizes, plus the log geometry, is the **entire**
input to the exponent fit and the curvature check. If you would rather not trust our
pipeline at all, you can rebuild all of Q1 from that one file.

**The exactly-solvable side** (Facts 5, and Q3):

```text
scripts/p398_reflection_parity.py           Gate 1: the reflection, the parity split, and
                                            selection_rule(), where the Duhamel integrand
                                            is evaluated pointwise in tau
scripts/p398_intervention_transport.py      #580, the frozen-realization transport
scripts/noncrossing_connectivity_codec.py   the state space itself
notes/p398-reflection-parity-20260906.md    Gate 1 and Phase C in prose
```

**The typed decomposition** (Fact 6, Q2), on `claude/p581-empirical`:

```text
notes/qtangent-empirical-20260907.md
scripts/score_qtangent_empirical.py
results/qtangent-empirical/latest.json
scripts/qtangent_scale_decomposition.py     (on the main workspace branch — the exact gate)
```

**Context you may or may not want:** `docs/ROADMAP.md` ranks work by information gained per
unit effort and is the closest thing to a narrative; `docs/STATUS.md` is the claim ledger;
`docs/astra/Q1`–`Q4` are four older questions in the self-contained mode, of which
`Q4-why-square-site-resists.md` is aimed at the same target as #610's Q4 and is arguably
better posed. `docs/PUBLICATION-CHECKLIST.md`, the cut-network no-go results (#435, #491,
#549, #550) and the exhaustive censuses are all there and none of the four questions touch
them.

---

## What is already settled, so you need not spend effort on it

Offered so you can skip it, not to fence it off. If you think any of it is wrong, that
outranks the questions.

- **The noncrossing reflection-fixed count is classical.** Callan–Smiley, arXiv:math/0510447,
  Theorem 1, plus Burnside. Ding 2016 supplies the Kreweras anti-isomorphism showing both
  even-`n` reflection classes give the same count. Cyclic sieving (Reiner–Stanton–White) is
  the *rotation* action and is not the reflection theorem.
- **"Coarsest lumping = automorphism-orbit partition" is not a theorem** without extra
  hypotheses. Kemeny–Snell Thm 6.3.2 is the criterion; D'Angeli–Donno 2013 Prop. 13 is a
  counterexample and Thm 12 gives conditions. Our `w = 4..8` agreement is an instance.
- **No named selection rule for equivariant Markov generators was found** in a deliberate
  search. The abstract vanishing is Schur / Wigner–Eckart; Hänggi 1978 covers symmetries of
  the master equation but not this. That gap is why #610's Q3 exists.
- **Fieller sample-size for a ratio of means is standard**; choosing among *observables* by
  a Fieller criterion was not found as a named method.

---

## What you can have

- **The whole repository**, including the branches above and anything on `main`.
- **All the underlying data.** The eight square-site productions are committed as raw batch
  histograms under `results/server-*/`, with full delete-one jackknife covariance
  reconstructible from them. Nothing in Q1 rests on a number you cannot recompute.
- **Compute, if you want a specific thing measured.** We have the exact P398 generator, the
  exact lumping refinement by partition refinement, a matrix-free memory kernel, exact
  enumeration on square bond to `L = 4`, Monte Carlo beyond that, and budget for roughly one
  new production run of a new site count. If your answer needs a number we do not have,
  name it precisely and we will get it.
- **Permission to reject the questions.** The four in #610 are what we currently believe is
  blocking. If the material suggests a better question, answer that one and say why. If a
  question is malformed — wrong object, false premise, a distinction we collapsed — saying
  so is more valuable than answering it as asked. If you want to go somewhere the four
  questions do not point, go.

---

## Three things we would rather you heard from us than discovered late

Observations, not instructions.

1. **If Q1's answer is "reconstruction artefact," it lives in
   `scripts/threshold_quantile_lineage.py`** — most likely `quantile` or `profile_cdf`. We
   argued the textbook quantile-estimator bias is about six orders of magnitude too small
   (`O(1/M)` at `M = 10^8`, against the `~7e-4` per-size bias the discrepancy needs). That
   arithmetic is the load-bearing step of our dismissal and it has not been independently
   checked.

2. **We named a channel `B_even`, meaning "duality-even."** Whether the construction earns
   that name is not something we have verified as carefully as the name implies. Q2 leans on
   its sign behaviour, so if the naming is doing work it has not earned, Q2's support is one
   item shorter.

3. **`notes/p582-amplitude-law-20260906.md` contains a correction we made to ourselves** —
   the section "The two lineages are not as independent as they look." We first read a 1.3%
   agreement as cross-lineage replication; it was shared structure at the third rung. We
   found that one. We do not assume it was the only one of its kind.

---

## What happens to your answer on our side

So you know the shape of the container, not to constrain what goes in it.

Your answer is recorded as a **theory input** in `docs/astra/ANSWERS.md`. It does not enter
`docs/STATUS.md`, does not score a frozen block, and does not by itself close a ticket —
the same standing as a literature packet. That is about our evidence accounting, not about
how seriously it is taken; several of the project's turns have come from exactly this
channel.

Two practical notes. If you propose a computation with a stated expected value and a stated
falsifier, we will run it and report the result whichever way it goes — that is the highest-
leverage form an answer can take here, and it is why the compute offer above is real. And
if an answer contains a claimed *exact threshold value*, we put it through
`scripts/threshold_claim_intake.py` before anything else is written down; the filter never
confirms, it only reports whether the claim is already dead against the committed
certificates and censuses.

It also helps to know roughly what you looked at — not as an audit, but because "I answered
this from the note without opening the reconstruction code" and "I read the code and the
artefact branch is still open" are different answers, and we would act on them differently.
