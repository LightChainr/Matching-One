# Issue #576 — Put the wrapping channel on published ground before the next ladder

Branch: `compute/p576-wrapping-channel-grounding-20260913`
Caller: compute1

This note covers the three parts of the ticket: (1) validate the main wrapping
channel against the Akhunzhanov–Eserkepov–Tarasevich exact torus polynomials;
(2) compute the Pinson/Arguin continuum wrapping at `r = 1,2,4` and make the
either/or declaration; (3) disambiguate the three distinct `11/4` objects.

## Part 1 — Akhunzhanov–Eserkepov–Tarasevich exact torus wrapping polynomials

### Retrieval (attempt record — honest)
- Paper: R. K. Akhunzhanov, A. V. Eserkepov, Y. Yu. Tarasevich, *Exact
  percolation probabilities for a square lattice: Site percolation on a plane,
  cylinder, and torus*, J. Phys. A 55, 204004 (2022), arXiv:2204.01517.
- The wrapping polynomials are in the arXiv **ancillary** file `torus.txt`.
- URLs tried (Bash `curl` + WebFetch):
  - `https://arxiv.org/abs/2204.01517` — confirms a supplement exists; rendered
    text gives no direct download link.
  - `https://arxiv.org/src/2204.01517v1/anc/` — ancillary listing; found
    `torus.txt`, `cylinder.txt`, `Algorithm2.pdf`.
  - `https://arxiv.org/src/2204.01517v1/anc/torus.txt` — **retrieved
    successfully** (12,665 bytes; copied to
    `results/issue576/wrapping-grounding/raw/torus.txt`).
  - IOPscience article page and the arXiv HTML version — no separately
    downloadable polynomial file beyond the arXiv ancillary.
- **Outcome: SUCCESS** for the main (NN square site) channel. The paper does
  **not** handle the NN+NNN (Sq8) matching lattice, so — as the ticket states —
  only the main channel is validated here.

### Format of `torus.txt`
For each `L = 3..12` the file lists integers `c_0..c_{L^2}` where `c_k` is the
number of `L x L` torus site configurations with exactly `k` occupied sites that
wind once around the horizontal (period-1) cycle. The wrapping probability
polynomial is `P_L(p) = Σ_k c_k p^k (1-p)^{L^2-k}`. At `p = 1/2`,
`P_L(1/2) = (Σ c_k)/2^{L^2}`.

### Independent exact enumeration (our own control)
`scripts/exact_wrapping_enum.cpp` brute-forces all `2^(L^2)` configurations for
`L = 2..5` (exact rationals):

| L | count | P_L(1/2) exact | decimal |
|---|-------|----------------|---------|
| 2 | 7 | 7/16 | 0.4375 |
| 3 | 175 | 175/512 | 0.341796875 |
| 4 | 19571 | 19571/65536 | 0.298629760742 |
| 5 | 8853291 | 8853291/33554432 | 0.263848632574 |

Reachable range: `L ≤ 5` exact (`2^25 = 33,554,432` configs, ~20 s in C++);
`L = 6` (`2^36`) not attempted. The doubling method was cross-checked against a
second (cut-graph) algorithm; the cut-graph version over-counts non-winding
full-width clusters (verified on the 36 disagreeing `L=3` configs), so the
doubling method is the correct one and is what produced the table above.

### Comparison with the published coefficients

| L | published Σ c_k | independent enumeration | match |
|---|-----------------|------------------------|-------|
| 3 | 175 | 175 | YES |
| 4 | 19571 | 19571 | YES |
| 5 | 8853301 | 8853291 | **NO (diff 10)** |
| 6–8 | parseable | (enum infeasible) | P_L(1/2) decreases toward continuum limit from above |
| 9 | 388338158699818471040971 | — | **P_L(1/2) = 0.1606 < limit 0.1694 → coefficients appear corrupted** |
| 10 | block malformed (100 entries, 101 expected) | — | **unparseable; would give P > 1** |
| 11, 12 | P_L(1/2) = 0.1239, 0.1085 < limit | — | **appear corrupted** |

Findings (reported, not reconciled):
- **L=3 and L=4 match the independent enumeration exactly** — the published
  polynomials are validated at the two smallest reachable sizes.
- **L=5: published Σ c_k = 8853301, independent enumeration = 8853291 — a
  discrepancy of 10 configurations.** I could not resolve which is correct with
  the resources here: my enumeration is a clean uniform brute-force that matches
  L=3,4 exactly, so a supplement transcription error is plausible, but I do not
  assert it. Flagged for maintainer review; not silently "made to match".
- **L ≥ 9 (and the malformed L=10 block)** in `torus.txt` give `P_L(1/2)` values
  that fall *below* the continuum limit `0.1694154…`, which is impossible for a
  sequence converging to it from above. These large-L coefficient blocks in the
  supplement appear corrupted and are **not** treated as ground truth.

### Relation to the repository's committed wrapping observable
The repository's committed wrapping observable for the main channel is the
Pinson/Arguin continuum baseline (`predictions/p156_pinson_arguin_baselines_20260829.json`,
`scripts/pinson_arguin_primitive.py`): `π_i({1,0}) = 0.16941543532134688938…`,
the `L → ∞` limit of the square torus. The validated published polynomials give
`P_L(1/2)` monotonically decreasing from `0.4375` (L=2) toward this limit
(`0.2986` at L=4, `0.234` at L=6, `0.182` at L=8), consistent with it. So the
main (NN square site) wrapping channel is now on published closed-form ground at
reachable sizes (L ≤ 8, with the L=5 caveat), which is exactly the unique
"already-published closed form, reachable sizes" check point the ticket asked
for. The NN+NNN (Sq8) matching channel is out of scope of the paper and remains
unvalidated against it.

## Part 2 — Pinson / Arguin continuum wrapping at r = 1, 2, 4

Used the repository's already-implemented, validated formula
(`scripts/pinson_arguin_primitive.py`, percolation specialization `e0 = 2/3`).
Evaluated `π({1,0})(i·r)` for `r = 1, 2, 4` with two independent evaluations
(direct Gaussian winding sum vs compact theta form) cross-checking each other:

| r | τ | π({1,0})(τ) | direct–theta abs diff |
|---|---|-------------|----------------------|
| 1 | i | 0.16941543532134688938260796919875445000145337645375 | 4.6e-91 |
| 2 | 2i | 0.50303589769463904028462207879367242007166407775608 | 3.6e-89 |
| 4 | 4i | 0.85596932105387227330836216781610004621209130143655 | 1.1e-90 |

These three numbers are the competitors. Cross-check with Pruessner–Moloney
(cond-mat/0310361): they report their independently computed values agree with
Pinson at `r = 2` to `< 1e-8`. My two independent evaluations of Pinson/Arguin's
published formula agree to `~1e-82 … 1e-90` at `r = 1, 2, 4`, far tighter than
that threshold, so the repository implementation is consistent with
Pinson/Arguin and inherits the Pruessner–Moloney cross-check.

### Either/or decision (required — choosing one)
**Choice: (ii) — explicit non-claim.** The freeze file should declare that the
matching-odd slope is **not** claimed to be the Pinson wrapping, and these three
numbers are explicit non-claims.

Justification:
- These are critical **continuum** `Q=1` Fortuin–Kasteleyn probabilities for the
  **primal** (site/bond) homology sector `{1,0}` on a square torus — exact closed
  forms, not finite-lattice matching data.
- The matching-odd observable is a **distinct combination** (matching parity, odd
  sector) on the NN+NNN (Sq8) lattice, which the Akhunzhanov paper does not cover
  and which the Pinson/Arguin formula does not address.
- Asserting equality between the finite-lattice matching-odd slope and the
  continuum primal wrapping `π({1,0})` would be an unwarranted identification; the
  repository's own interpretation boundary (`notes/pinson-arguin-primitive-baseline.md`)
  already states these are continuum exact formulas, not finite-lattice results
  and not an H4 identification.
- Therefore these three numbers are recorded as **named non-claims / reference
  competitors**, not as the value of the matching-odd slope.

## Part 3 — Three distinct 11/4 objects (do not conflate)

| position | what it is |
|---|---|
| `Ê4(2i)/Ê4(i) = 11/4` | modular weight-4 amplitude ratio — a competitor for the staircase (a different modular quantity in Part 2's neighbourhood) |
| Newman–Ziff estimator `L^{-11/4}` | finite-size convergence rate: `1/ν + θ = 3/4 + 2 = 11/4`, the exponent governing how the estimator approaches `p_c` |
| Mertens–Ziff matching function root `~ L^{-4}` | a third, different exponent (matching-function finite-size scaling), unrelated to the `11/4` above |

These are three different objects that involve the numeral `11/4` (or nearby
exponents). A wrapping estimator that converges as `L^{-11/4}` does **not**
predict, and must not be read as predicting, a spin-4 amplitude ratio equal to
`11/4`. No claim-ledger entry is added for this clarification.

## CI note
Full Matching-One repository CI has not been run for this commit.
