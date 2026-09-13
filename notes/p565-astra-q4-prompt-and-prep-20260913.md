# #565 — Astra Q4 prompt freeze and preparation note

**Branch:** `retrieval/p565-astra-q4-20260913`
**Ticket:** #565  ·  **Source:** `docs/astra/Q4-why-square-site-resists.md` (read in full via `gh api`; exists on `main`)

This ticket asks me to (1) **freeze the prompt** that the repo intends to send to an external
mathematical model, (2) verify the "do-not-spend" list against the ticket's update evidence and
independently confirm the Ziff-2006 A-lattice quintic is irreducible of height 4, and (3) record that
the literal answer is **absent** (I have no external-model channel) and that this note is **not** an
answer to Q4.

---

## 1. Prompt-freeze record

| Field | Value |
|-------|-------|
| File | `docs/astra/Q4-why-square-site-resists.md` (on `main`) |
| SHA-256 | `785d9f42594070045807de2be851c59a8841722e632cc2cdbb3800680142777f` |
| Characters | 6972 |
| Lines | 122 |
| Self-contained? | **Yes** — the file is a complete, standalone mathematical question (table of known thresholds; the regularity question; the "why the answer decides something" table; the do-not-spend list; provenance pointers). |
| Violates "do not concatenate other question files"? | **No.** Sibling files `Q1/Q2/Q3/Q5/Q6` exist in `docs/astra/`, but Q4 contains **no instruction** to fetch or concatenate them. Its "Provenance of the framing above" cites *repo artifacts* (`manuscript.md`, `results/*.json`, `analysis/pslq_search_contract.json`) as claims the model "does not need re-derived" — these are data references, not other question prompts, and are not required to be read for the prompt to stand alone. |
| References `docs/astra/ANSWERS.md`? | No (that filename appears only in the ticket's workflow instruction, not in the prompt file). |

## 2. do-not-spend list — consistency check

The file's "Do not spend output on" list reads:
1. Re-deriving the known thresholds in the table.
2. A survey of percolation / exact solvability in statistical mechanics.
3. Numerical estimation of the threshold (already at 13 digits).
4. Arguing the threshold is "probably transcendental" without a mechanism.

**Consistency with the ticket's update evidence:** the ticket states the update evidence is that
"degree ≤ 6 **and** height ≤ 3" is false, because Ziff 2006's A-lattice threshold is an irreducible
quintic of **height 4**, `p⁵ − 4p⁴ + 3p³ + 2p² − 1`. The file's threshold table indeed lists the
A-lattice row with degree 5 and height **4**, and the prose states the "degree ≤ 6 and height ≤ 3"
census class was broken by this row. → **The do-not-spend list and the update evidence are mutually
consistent** (the file does not ask the model to re-derive the A-lattice row; it uses it as the
counterexample). No inconsistency found.

## 3. Independent exact verification — Ziff 2006 A-lattice quintic

Polynomial (as printed in the file): `f(p) = p⁵ − 4p⁴ + 3p³ + 2p² − 1`.

**Claim to verify:** irreducible over ℚ, degree 5, coefficient height 4.

### 3.1 Degree and height (trivial)
- `deg f = 5`.
- coefficients (descending): `1, −4, 3, 2, 0, −1` → max absolute value = **4**. ⇒ height = 4. ✓

### 3.2 Irreducibility — rational-root test (no linear factor)
Rational-root candidates are divisors of the constant term over divisors of the leading coefficient:
±1.
- `f(1)  = 1 − 4 + 3 + 2 + 0 − 1 = 1  ≠ 0`
- `f(−1) = −1 − 4 − 3 + 2 + 0 − 1 = −7 ≠ 0`
No rational root ⇒ no linear factor over ℚ.

### 3.3 Irreducibility — exclude quadratic × cubic (exact, by Gauss's lemma)
Since `f` is monic, any factorization over ℚ into degree-2 × degree-3 is, up to sign, into **monic
primitive integer** polynomials. Write
`f(p) = (p² + p·x + q)(p³ + r·p² + s·p + t)` with `q,t ∈ ℤ`, `q·t = −1` (constant term) ⇒ `q = ±1`,
`t = −1/q`. Matching coefficients:
- `x + r = −4`
- `q + x·r + s = 3`
- `q·r + x·s + t = 2`
- `q·s + x·t = 0`
- `q·t = −1`

From `q·s + x·t = 0` with `t = −1/q`: `q·s = x/q` ⇒ `s = x` (since `q² = 1`).
From `x + r = −4`: `r = −4 − x`.
Plug into `q + x·r + s = 3` ⇒ `q + x(−4 − x) + x = 3` ⇒ `q − 3x − x² = 3` ⇒ `x² + 3x + (3 − q) = 0`.
- **Case q = 1:** `x² + 3x + 2 = 0` ⇒ `x = −1` or `x = −2`.
  - `x = −1`: `s = −1`, `r = −3`. Check `q·r + x·s + t = (1)(−3) + (−1)(−1) + (−1) = −3 + 1 − 1 = −3 ≠ 2`. ✗
  - `x = −2`: `s = −2`, `r = −2`. Check `= (1)(−2) + (−2)(−2) + (−1) = −2 + 4 − 1 = 1 ≠ 2`. ✗
- **Case q = −1:** `x² + 3x + 4 = 0` ⇒ discriminant `9 − 16 = −7 < 0` ⇒ no integer `x`. ✗

No monic quadratic × cubic factorization exists. Combined with §3.2, **`f` is irreducible over ℚ**. ✓

### 3.4 Machine re-check
A brute-force enumeration over monic quadratics `x² + p x + q` (`q = ±1`, `p ∈ [−10,10]`) with monic
cubic quotient found **no** factorization (Python, exact integer arithmetic). Consistent with §3.3.

**Conclusion:** the A-lattice threshold polynomial is an **irreducible quintic (degree 5) of height 4**,
exactly as the file and the ticket claim. This **confirms** that "degree ≤ 6 **and** height ≤ 3" is
false (degree 5 ≤ 6 holds, but height 4 > 3), validating the counterexample the ticket relies on.

> Provenance note: the file states its A-lattice row is "corroborated from two independent indexes
> rather than read from the primary text." The irreducibility/height proved here is a property of the
> printed polynomial itself, independent of that provenance.

---

## 4. Literal answer — ABSENT (required explicit statement)

- **The literal answer to Q4 is missing from this note.** I have **no channel** to the external
  mathematical model the ticket describes; I cannot execute, invoke, or receive output from it.
- **This note does NOT constitute an answer to Q4.** It is a prompt-freeze + preparation record only.
- **No external-model output has been fabricated.** Every statement above is either a file-property
  measurement (SHA-256, counts), a consistency observation, or an exact arithmetic verification I
  performed locally. Where the file makes claims (e.g. the regularity heuristic, the census results),
  they are reported as the file's claims, not as established fact.

## 5. Mandatory disclaimer
This note performs no novelty judgement and answers no mathematical question about square-site
percolation. The irreducibility proof in §3 is exact and self-contained. The Q4 answer remains the
external model's job, which I could not perform.

*Full Matching-One repository CI has not been run for this commit.*
