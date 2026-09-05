# Governance

Matching One is an exploratory computational-mathematics project. Its output is
mathematics. Everything else — checks, manifests, digests, templates, ledgers — is
overhead, justified only when it buys more discovery than it costs.

Time spent on assurance is time not spent exploring. That trade is real and it has
gone the wrong way here before.

## 0. Two speeds

This is the rule the rest of the document serves.

| | **Exploring** (default) | **Publishing** (a specific paper, a specific claim) |
|---|---|---|
| Verification | minimum necessary — §2 | full — `docs/PUBLICATION-CHECKLIST.md` |
| Freezing, preregistration | only when the run is meant as confirmatory | as required by the claim |
| Digests, manifests, provenance chains | no | yes, for the artifacts the paper cites |
| PR ceremony | a sentence saying what and why | the checklist |
| Reviewer | none needed | whatever the venue needs |

**Everything that reads like a requirement in this repository applies at publication
time unless it appears in §2.** If you find a rule elsewhere that seems to gate
exploratory work, it is stale; delete it.

The order matters: explore first, verify what survived. Verification is *cheaper*
after the fact, not more expensive, because by then you know which numbers mattered.
Verifying during exploration means verifying everything, including the nine ideas
out of ten that were about to be abandoned.

## 1. Default mode: run and integrate

`main` is the shared research line, not a publication-only branch.

Useful scripts, exact calculations, source-data reanalyses, pilots, theory notes,
frozen predictions, negative results, and result archives enter `main` as soon as
they are understandable. External approval is not required.

A registry or documentation conflict must never block a scientifically useful asset.
Integrate the asset; repair navigation later, or not at all.

Branches and PRs are coordination tools, not permission gates.

## 2. The minimum, while exploring

These five are the whole list. They are here because each is either free, or
impossible to reconstruct later, or the kind of mistake that yields a confidently
wrong answer instead of a missing check.

**A. Don't fool yourself about a number.** If a result is claimed exact, compute it
exactly — exact arithmetic, no binary floating point in the claim. If a computation
*is* the result, check it once, by the cheapest independent means available: a closed
form, a tiny case done by hand, a different method on one input. Once. Not a suite.

**B. Don't destroy data.** Commit raw sufficient statistics, not just fitted
coefficients. Don't overwrite a committed result; add alongside it. Reruns are
expensive; disk is not.

**C. Don't misdate.** If something was frozen before a target was seen, say so. If it
was not, don't imply it was. This is a sentence, not a process — but it is the one
fact that cannot be recovered afterwards, which is why it survives the cut.

**D. Say which observable.** When comparing two numbers, name the convention each
uses. If they differ, say what the map is. A wrong comparison is worse than no
comparison and much harder to notice later.

**E. Count one random block once.** Roots, slopes, derivatives, quantiles and score
modes taken from the same histograms are all useful, and they are all the same block.
Adding them as independent evidence inflates the conclusion. This one stays in the
minimum not because it is ceremony-free — it is — but because violating it produces a
confidently wrong answer rather than a missing check.

That is all. Not: digests, tamper tests, manifests, registries, adoption audits,
templates, power analyses, or independent reimplementations. Those are §3.

## 3. What moves to publication time

Everything else. Concretely, and non-exhaustively: SHA-256 chains and provenance
manifests; independent second implementations; preregistration and held-out design
for claims that need them; full covariance treatment; power and sensitivity; digest
re-verification; artifact registries; reviewer checklists; the PR template's long
form; and anything phrased as "every X must have Y".

`docs/PUBLICATION-CHECKLIST.md` holds the full form. Reach for it when there is a
specific paper and a specific claim. Not before.

## 4. What not to build

The failure mode here is additive: every individual check looks reasonable, and the
sum is a compliance system with a research project attached. Some checks are always
a bad trade, so they are named:

- **Tamper tests.** Verifying that a validator raises when you deliberately corrupt
  its input. Nobody is corrupting anything. This tests an error path, not the
  mathematics.
- **Error-path tests for frozen constants** — "fails closed when a manifest field
  drifts". The constant is right there in version control.
- **Tests that assert a document contains a sentence.** Prose is not an invariant.
- **Digest re-verification.** A hash proves a file did not change. It never proved
  the file was right, and nothing here changes files behind our back. Git already
  does this job.
- **Meta-tooling that audits the repository itself** — inventories of which files
  import which type, adoption reports, registry-consistency checkers.
- **Wrappers whose only function is to be audited.**

A test earns its place if you can name the wrong number it would prevent us from
believing. If you cannot name that number in one sentence, do not write the test.

Remember that a check written today runs forever. Its cost is not the hour you spend
writing it; it is the permanent tax on every run afterwards, and the pressure it puts
on everyone later to keep the pattern going. That is why "it's cheap, might as well"
is wrong.

## 5. Scientific claim levels

| Level | Meaning |
|---|---|
| C0 | hypothesis, conjecture, design, or theory candidate |
| C1 | method/control validated by exact identity, oracle, or deterministic regression |
| C2 | exploratory numerical signal; analysis may be adaptive |
| C3 | reproduced/frozen finite-size numerical result — independent seed, prospective or held-out |
| C4 | asymptotic/mechanistic interpretation supported by multiple discriminating tests |
| C5 | rigorous result or independently checkable certificate/proof |

A result can be on `main` at any level. Merging is not a claim upgrade, and lack of
preregistration is not a reason to discard useful C2 evidence.

Claim levels are labels on conclusions, **not** gates on work. Nothing needs to reach
a level before it can be committed, discussed or built on. The only thing a level
governs is the language used in `docs/STATUS.md`.

## 6. Research execution

Existing-data analysis, exact calculations, controls, pilots, and new production:
**run them**. There is no production gate. Choose by expected information gain.

For an expensive confirmatory question, freezing the target, model, sign and score
before reading the target is worth it, because it is what earns C3. If that was not
done, the run is still useful C2 evidence rather than something to discard.

Large campaigns — GPU, Pell, N=1105, norm-4, norm-5, modulus scans — are **priority
decisions, not permission classes**.

## 7. High-risk machinery

Topology, homology, RNG, threshold-rank reconstruction, covariance propagation and
exact polynomial machinery contaminate everything downstream if wrong. These deserve
a real check — one, at the point of reuse, proportional to consequence.

This is the one place where more than the §2 minimum is worth it while exploring, and
it is worth it because the blast radius is large, not because checking is virtuous.

## 8. Results and corrections

Negative, null, failed, underpowered and contradictory results are first-class assets.
Keep them.

On finding an error: preserve the old artifact, add the correction, say what changed.
Update `docs/STATUS.md` only if the claim boundary moved. Do not build a correction
workflow more elaborate than the risk.

## 9. Scientific language

Language stays conservative even though execution is permissive:

- "observed at these sizes", not "asymptotic", without an asymptotic test;
- "compatible with", not "equal to", for a numerical candidate;
- "candidate operator", not "identified operator", until competitors are excluded;
- "exact" only for a proved identity, exact arithmetic, or a certified computation.

A result may be valuable because it kills a promising route. Do not tidy failures out
of the narrative.

## 10. Operating principle

**Explore first. Verify what survives.**

- integrate useful analysis immediately;
- use existing sufficient statistics harder before assuming new data are needed;
- run cheap exact controls early — to kill ideas, not to certify them;
- choose expensive work by information gain, never by ceremony;
- let a failed test redirect the programme the same day;
- when a rule and a discovery conflict, the rule is what gives way.
