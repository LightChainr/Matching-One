# Astra answers

One section per query. Nothing here is evidence about the lattice — see the
epistemic rule in [`README.md`](README.md). Record the answer before assessing it.

## Template

```markdown
## Qn — <short title>

- **Asked:** YYYY-MM-DD
- **Model / version:**
- **Prompt sent:** `docs/astra/Qn-....md` at commit `<sha>`, unmodified / with these changes: …
- **Cost:**

### Answer (verbatim)

> …

### Our assessment

- **Verification status:** VERIFIED_BY_US / PARTIALLY_VERIFIED / NOT_VERIFIED_BY_US
  — with what we checked and how.
- **Which branch of the decision rule this triggers:**
- **Disagreements:** parts we think are wrong, and why. Keep these even if we later
  turn out to be the ones who were wrong.

### What was done about it

- Issue / PR / commit, or "nothing yet".
```

## Log

### #610 — 2026-09-07, Codex / GPT-6

[Answer and local calculations](ANSWER-610-20260907.md), based on workspace
`8b5f9d1a`. The 55% curvature comparison changes affine chart; transporting the
chart explains most of its excess. Includes exact square-bond covariance
identities, the transfer-matrix selection lemma, and a conditional quantile
convergence proof for #276. Local verification and unverified scope are stated
in the answer. Independent team assessment is pending; no C3/C4 promotion.
