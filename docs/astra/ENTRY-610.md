# Entry prompt for #610 — what to read, in what order, and how to answer

**This file changes the mode of the pack.** Every other file in `docs/astra/` is written
for a model that **cannot be shown the repository**: self-contained, one file, no context.
[#610](https://github.com/LightChainr/Matching-One/issues/610) is written for one that
**can** read it. If you are answering #610, read this first; if you are answering Q1–Q4,
ignore this file entirely and read only the one you were given.

**Branch.** Read `claude/matching-one-workspace-pwr5pv`, not `main`. `main` is several days
behind and none of the four gates #610 is about have landed there. Two files live on other
branches and are marked below.

**Budget.** The minimum honest path is §1 plus §2c–d, about 40 KB. The full path is about
250 KB. §2 is the only section where reading the *code* matters; everywhere else the note
is enough.

---

## §1 — Orientation. Four things, ~25 KB. Do not skip.

Read in this order:

1. **`GOVERNANCE.md`, section 2 only** (lines 44–78). Five minimums. Two of them decide how
   to read everything else: *don't fool yourself about a number — check once by independent
   means*, and *count one random block once*. When #610 says a check "fired," it means the
   independent check under minimum A returned a contradiction and was reported rather than
   absorbed.
2. **`docs/STATUS.md`, the "Strongest current evidence" table only.** The `Level` column is
   the grammar of this project: `C3` is a scored prospective test, `C2` is reanalysis of
   committed productions, and a row can be *negative* and still be the strongest thing in
   its line. Note how many rows read "negative." That is the house style, not modesty.
3. **`docs/astra/README.md`.** The four-part selection rule for what belongs in this pack,
   and the requirement that a question's decision rule differ between yes and no. #610 is
   held to it too — if one of its questions fails that test, say so.
4. **The body of #610 itself.** Facts 1–6 and Q1–Q4. That is the question; everything below
   is evidence for it.

**Stop here and check.** If after these four you believe one of the questions is malformed
— wrong object, false premise, a distinction I have collapsed — say that now, before
reading further. That correction is worth more than an answer to the wrong question, and it
costs the least.

---

## §2 — Q1: the tangent that obeys one exponent and the curvature that refuses it

Q1 is the only question where the repository holds data you can overturn, so read the
**definitions before the numbers, and the numbers before my prose**. The order matters: the
definitions decide whether the numbers mean what I claim.

**a. `scripts/threshold_quantile_lineage.py`** (381 lines) — read `rank_cdf`,
`_binomial_weights`, `profile_cdf`, `quantile`, `jackknife_quantiles`,
`jackknife_covariance`, `spin_zero_weights`, `is_interpolation`.

This is the file where a reconstruction artefact would live. The threshold CDF is an exact
binomial convolution of the empirical rank CDF, bisected; the jackknife deletes a batch from
every level at once. **If your answer to Q1(b) is "it is an artefact," this file is where
you have to point, and `quantile` / `profile_cdf` are the two candidates.** I have argued
the textbook quantile-estimator bias is ~6 orders of magnitude too small; check that
arithmetic, it is the load-bearing step of my dismissal.

**b. `scripts/score_wasserstein_shape_flow.py`** (510 lines) — read `shape_flow`,
`load_sizes`, `cross_size_coupling`, `same_production_null_control`.

`shape_flow` is the affine projection: basis `[1, Q_base]`, covariance
`(S_base + S_target)/log(m)^2`, weighted by the spectral pseudo-inverse.
`cross_size_coupling` measures the batch-by-batch correlation between the two sizes rather
than assuming the cross term away — two lineages reuse a seed across sizes, which is why it
exists. `same_production_null_control` scores two halves of one production as if they were
a transition; it must come back consistent with zero, and does.

**c. `notes/p582-amplitude-law-20260906.md`** (218 lines) — the result in prose, including
the four-cell weighting table and the section "The two lineages are not as independent as
they look," which is a correction to my own first reading. Read that section carefully: it
is the part most likely to contain a second error of the same kind.

**d. `results/p582-amplitude-law/latest.json`** (544 lines) — the numbers. Keys worth
opening: `amplitudes`, `one_exponent_fit`, `leave_one_transition_out`,
`second_difference_check`, `weighting_systematic.second_difference_ratio_cells`,
`error_budget`, and `not_established` (which lists what I already believe is unestablished,
so you need not spend the answer telling me).

**e. `scripts/p582_amplitude_law.py`** — only if you dispute a specific number. The four
functions that carry Q1 are `first_difference_image`, `second_difference_weights`,
`second_difference_image`, and `second_difference_check`. The first is the exact
finite-difference image of an exponential, *not* its derivative; the second's weights
annihilate constants; the fourth freezes the parameters fitted on first differences and
never refits.

**f. `results/wasserstein-shape-flow/latest.json`** (2242 lines) — the raw material.
`sizes.<N>.quantiles_spin0`, `.standard_errors_spin0`, `.spin0_weights`,
`.orientation_cos4theta`, `.source`, `.batches`.

**You can rebuild the entire Q1 dataset from this one file.** Nine quantiles and their
standard errors at eight sizes, plus the log geometry, is everything the exponent fit and
the curvature check consume. If you would rather not trust my pipeline at all, this is the
offer: recompute from here and tell me what you get.

---

## §3 — Q2: is the square-site flow the even sector of an involution

1. **`docs/STATUS.md`, "Exact semantics and controls."** The two exact finite relations the
   whole conjecture rests on: `DeltaS_cross = -DeltaS_either` (corrected score `0.5700/2`)
   and the finite Russo identity
   `M'(p) = pivotal_mass_primal(p) + pivotal_mass_matching(1-p)`. These are exact at finite
   `N`, not asymptotic. If the involution is going to act anywhere, it acts here.
2. **Gate 4 — on branch `claude/p581-empirical`:** `notes/qtangent-empirical-20260907.md`
   (~6.6 KB), then `results/qtangent-empirical/latest.json` if you want the per-scale
   numbers. This is the even/odd-looking sign pattern: `X` single-signed and monotone,
   `B_even` alternating and 6–11× smaller.
   **Warning, and please act on it:** I called that second channel `B_even`, "duality-even."
   Check whether the *construction* earns the name before you let the naming do any work in
   your argument. If it does not, Q2's supporting evidence is one item shorter and I want
   to know.
3. **`notes/p398-reflection-parity-20260906.md`** (177 lines) — the exactly-solvable
   analogy in full, including the Phase C section at the end (the balanced realization
   factors through the quotient; the odd sector is exactly inert). This is what the
   square-site structure is being compared *to*.

---

## §4 — Q3: state the selection rule as a lemma

1. **`notes/p398-reflection-parity-20260906.md`** — already read in §3.
2. **`scripts/p398_reflection_parity.py`** (1016 lines) — do not read it all. Read
   `reflection`, `state_permutation`, `parity_split`, `selection_rule`, `responses`,
   `second_order`. `selection_rule` is the exact object: it is where the Duhamel integrand
   is evaluated pointwise in `tau` and returns `5.5e-16`.
3. **On branch `docs/literature-officer-20260906-issue601`:**
   `notes/literature-officer-20260906-issue601.md` — read **the decision table at the top
   and the Q2 and Q3 rows only** (the file is 39 KB and most of it is bibliography).

   **Do not re-derive Callan–Smiley.** The reflection-fixed count is already cited and
   settled. The live gap is exactly the Q3 row: *no named selection rule for equivariant
   Markov generators was found.* That is what #610's Q3 asks you to state.

---

## §5 — Q4: what any of this is worth

1. **`docs/ROADMAP.md`** — 435 lines; do **not** read it whole. Read the "Theory-bound
   questions" section at the top, then search for `Gate 3` and read the block that follows
   it through the `N = 725` paragraph.
2. **`docs/astra/Q4-why-square-site-resists.md`** (122 lines) — an older question, still
   unasked, aimed at the same target: whether the low degree/height of every exactly-known
   threshold is a theorem about the three solvable mechanisms, and whether square-site is
   provably outside them.

   **If your Q4 answer would duplicate that one, say so and answer that one instead.** It
   is the better-posed version and it was written first.

---

## §6 — How to answer

**Order.** Answer Q1 first and as a separate unit. It is the only question where I hold
data you can overturn, and its answer changes what Q2 and Q4 are even about: if the 1.5 is
an artefact, the "second smooth scale" reading of the #582 remainder dies and Q2 loses half
its motivation.

**Shape, per question.**

```text
verdict          one line, plain, before any argument
argument         with hypotheses stated as hypotheses
status           proof / conditional proof / conjecture / heuristic -- mark every claim
what it changes  what I should do differently, in one sentence
```

**The single most valuable thing you can produce is a computation I can run.** Not "this
would be consistent with two exponents" but:

```text
compute:      <quantity>, from <committed file / exact object>
expected:     <value or range> if the claim holds
falsified if: <value or range>
```

I have the exact P398 generator, exact lumping refinement, the matrix-free memory kernel,
exact enumeration on square bond to `L = 4`, the eight committed square-site productions
with full jackknife covariances, and budget for about one new production. A prediction I
can kill in an afternoon outranks a correct but unfalsifiable argument.

**If you find an error in my arithmetic, lead with it.** That outranks every question here.
The four-cell correction in §2c is one I found myself; assume there is another.

**Negative answers are the point.** "Q1(b): a log correction is not excluded and here is
why your `s`-dependence argument fails," "Q2: the matching involution does not act on this
functional and the resemblance is superficial," "Q4: that is a category error and the
shortest path is X" — each of those closes a line of work and is worth more than a
hedged yes.

**Constraints.**

- Your answer is a **theory input**. It does not enter `docs/STATUS.md`, does not score a
  frozen block, and does not close a ticket. It is recorded in `docs/astra/ANSWERS.md`.
- If your answer contains a claimed **exact threshold value**, it goes through the filter
  before anything else is written down:
  `python3 scripts/threshold_claim_intake.py --expression "<closed form>"` (or
  `--polynomial` / `--decimal`). It never confirms; it only tells us whether the claim is
  already dead.
- **Say what you did not read.** If you answered Q1 from the note without opening
  `threshold_quantile_lineage.py`, say so — then I know the artefact branch of Q1(b) is
  unexamined rather than dismissed.
