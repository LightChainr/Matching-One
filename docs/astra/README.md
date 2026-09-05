# Astra query pack

A small set of questions for an expensive external mathematical model (GPT‑6 Astra).
The model cannot be shown the repository. It can be shown one of these files.

## The selection rule

The cost of a query is high enough that "it would be interesting to know" is not a
reason to ask. A question belongs here only if all four hold:

1. **Theory-bound.** No amount of compute or Monte Carlo we can run would answer it.
   If more samples would settle it, run the samples instead.
2. **Blocking.** Something in the project is stopped, and the answer restarts it —
   or tells us to abandon it, which is equally useful.
3. **Self-contained.** A mathematician who has never seen this repository can answer
   from the file alone. If the question needs our data, it is the wrong question.
4. **Decision-changing in both directions.** We can write down, before asking, what
   we do if the answer is yes and what we do if it is no, and those differ.

Each file states its decision rule explicitly. If a question's decision rule reads
"we would carry on as before" under either branch, delete the question.

## The questions

| | Question | Blocks | Status |
|---|---|---|---|
| [Q4](Q4-why-square-site-resists.md) | Is "degree ≤ 6, height ≤ 3" a theorem about the solvable mechanisms, and is square-site provably outside them? | the interpretation of our main exclusion result, and the original target | unasked |
| [Q1](Q1-descendant-log-coupling.md) | Does the `h=5/8` logarithmic coupling descend to the level‑4 `x=21/4, s=4` state? | issue #275, the project's only P0 | unasked |
| [Q2](Q2-additive-shape-ambiguity.md) | What admissible additive shape carries the measured ratio from `11/4` to `1.88 ± 0.18`? | `ROADMAP` item 2 — **the block has been run and came back negative** | unasked, and sharper than when written |
| [Q3](Q3-unit-automorphism-escape.md) | What is the smallest escape from the Gaussian-unit no-go? | the Gaussian-cover production line, declared dead by a theorem | unasked |

**Q4 first.** It is the one aimed at the problem the project was started for, it is
the only one whose answer changes a manuscript already drafted, and it is the only one
where a lucky answer could be very large. Q1 is the project's P0 and comes next. Q2
and Q3 are independent of both and of each other.

Q2 has *improved* since it was written: the production block it was asking us to
approve has since been frozen, run and reported negative, so it now carries a
measured number to explain rather than a design to bless.

## If an answer contains a claimed exact threshold

Run it through the filter before doing anything else:

```bash
python3 scripts/threshold_claim_intake.py --polynomial <ascending integer coefficients>
python3 scripts/threshold_claim_intake.py --decimal <as many digits as are claimed>
python3 scripts/threshold_claim_intake.py --expression "<mpmath closed form>"
```

It places the claim against the four published intervals and both exhaustive censuses
and returns one of four verdicts: refuted by a committed certificate, already
catalogued as a width artifact, contradicts every published interval, or survives our
checks. **It never confirms** — the last verdict means only that the claim is not
already dead, and the report lists what would still have to be done.

The filter takes seconds and costs nothing, so it goes first, before any assessment
is written into `ANSWERS.md`.

## What an answer is, and is not

An answer from Astra is a **theory input**, at the same epistemic level as a
derivation we did ourselves or a claim in a paper we have not checked. It is not
evidence about the lattice.

Concretely, an answer may **not** be cited as support for any statement in
`docs/STATUS.md`. What it may do is supply the missing column, contract or
prediction vector that an existing frozen dataset can then score. The order is:

```text
answer -> a candidate prediction vector with units and a normalizer
       -> profile-rank / covariance score against an existing frozen block
       -> only then a claim, at whatever level that score earns
```

This ordering is what keeps a confident paragraph from an expensive model out of
the claim ledger. It applies even when the answer is a proof: a proof about the
continuum is still a hypothesis about our lattice observable until the overlap is
established, which is the substance of Q1 in the first place.

If an answer contains a derivation we cannot follow, that is a result too — record
it as `NOT_VERIFIED_BY_US` and do not build on it.

## Recording answers

Append to [`ANSWERS.md`](ANSWERS.md), one section per query, using the template
there. Record the answer verbatim before any commentary, including the parts we
think are wrong. Then, separately, our own assessment.

Do not edit a question file after it has been asked. If the answer shows the
question was malformed, add a follow-up file (`Q1b-...`) explaining what was wrong
with the original — the failure mode is itself worth keeping.

## Asking economically

- Send **one** question file as the whole prompt. Do not concatenate.
- Do not paste the repository, the issue threads, or the data.
- Every question file ends with a **Do not spend output on** list. It is there to
  stop the model re-deriving things we already have, which is where the budget
  actually goes.
- If a follow-up is needed, quote only the specific sentence being followed up on.
