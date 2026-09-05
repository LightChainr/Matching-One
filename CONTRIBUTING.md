# Contributing

Read `GOVERNANCE.md` §0 and §2 first. They are short, and they are the whole rule set
for exploratory work.

Everything in this file that sounds like a requirement applies **at publication
time**, and lives in `docs/PUBLICATION-CHECKLIST.md`. While exploring, the list is:

1. don't fool yourself about a number;
2. don't destroy data;
3. don't misdate a freeze;
4. say which observable;
5. count one random block once.

## While exploring

**Issues.** Open one when you want the discussion or want to hand work off. Not
required before starting. An issue that exists only to record permission is noise.

**Branches.** One focused branch, any name that says what it is. Don't stack.

**Committing.** Say what you did and why in the commit message. That is the
documentation. If the work is interesting, a note in `notes/` is worth more than a
manifest.

**Pull requests.** A sentence on what and why. The template is a *publication-time*
form; for exploratory work, delete the sections that don't apply, or don't use it.
Nobody is waiting to review — merge your own work when it is useful.

**Checks before pushing.** Whatever convinces *you* the number is right. Usually:

```bash
python3 -m compileall -q scripts tests
python3 -m unittest tests.test_<the_thing_you_changed>
```

Run the whole suite when you have reason to think you broke something far away —
not as a ritual before every push, and not to feel finished.

## Writing tests

The bar: **name the wrong number this test would stop us believing.** One sentence.
If you can't, don't write it.

Worth writing:

- an exact vector or closed form the implementation must reproduce;
- an independent method agreeing on one input — union-find against BFS, a tiny case
  done by hand;
- an invariance the mathematics requires — basis, geometry, batch partition;
- a leakage check where held-out integrity is the point.

Not worth writing, ever — see `GOVERNANCE.md` §4:

- tamper tests and other validator error paths;
- "fails closed when a frozen constant drifts";
- assertions that a document contains a sentence;
- digest and checksum re-verification;
- anything auditing the repository's own structure.

A test that repeats the implementation in different syntax protects nothing.

## Code

**Python.** Standard library for reference implementations. Type hints where they
help a reader. Exact arithmetic wherever exactness is part of the claim — never
binary floating point inside an exactness claim.

Validate an input when a bad value would silently produce a plausible wrong answer.
Do not validate inputs to satisfy a policy: a `ValueError` for a negative CLI
argument protects nothing.

**C++.** RNG and reductions independent of thread scheduling where promised.
Fixed-width integers for counters, with the overflow bound written down. No
`fast-math` where it can move a Bernoulli decision or a rank ordering. Keep a slow
reference path — that is the check that earns its place.

## Data and results

Text formats: CSV, JSON, YAML, Markdown, exact integer coefficients.

An imported dataset needs a citation precise enough to find the table again, and the
decimals as printed. That is it, while exploring. Checksums, row-count tests and
sentinel tests are publication-time.

One thing that is *not* optional, because it cannot be recovered later: if you
transcribe numbers from a source, check the digits against the source **once**, when
you transcribe them, and say in the commit that you did. A hash cannot do this — it
proves a file did not change, never that it was right to begin with.

Commit raw sufficient statistics, not only fitted coefficients. Don't overwrite a
committed result; add beside it and link the superseded one.

## Experiments

If a run is meant as confirmatory evidence, freeze the hypothesis, sizes, orientation
order, splits, seeds, estimand and stopping rule before you look. That is what earns
C3, and it is the reason to do it — not compliance.

If it is exploratory, run it. Report signed effects and uncertainties. Keep the
negative results.

## Scientific writing

Precise language, per `GOVERNANCE.md` §9. "Compatible with", not "equal to".
"Candidate", not "identified". "Exact" only when it is.

Do not remove failed attempts to improve the narrative.

## Review

There is usually no reviewer. When there is, look in this order: definitions and sign
conventions; whether the number is right; whether the claim language matches the
evidence. Then stop.

Correctness matters. Auditability is a property of published work, not of a research
line in progress.
