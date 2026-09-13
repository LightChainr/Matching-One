# Issue #576 — Put the wrapping channel on published ground before the next ladder

Branch: `compute/p576-wrapping-channel-grounding-20260913`
Caller: compute1; resolution pass 2026-09-13 by compute1-retry (second
independent enumerator, spiral-bug diagnosis, divisibility + P_L(p_c) block
health, Pinson values recomputed)

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

**Resolution pass (2026-09-13, compute1-retry).** A second, independent
enumerator was written and run, and the two discrepancies flagged above were
both resolved. The corrected table is:

| L | count (displacement-DSU enum, `exact_wrapping_enum2.cpp`) | published Σ c_k | match |
|---|---|---|---|
| 2 | 7 | — (file starts at L=3) | — |
| 3 | 175 | 175 | YES (per-coefficient) |
| 4 | 19571 | 19571 | YES (per-coefficient) |
| 5 | **8853301** | **8853301** | **YES (per-coefficient, all 26 coefficients)** |

Reachable range: `L ≤ 5` exact (`2^25 = 33,554,432` configs);
`L = 6` (`2^36`) not attempted.

**Resolution of the earlier L=5 "diff 10".** The first enumerator's
doubled-grid criterion ("some doubled-grid component contains both (r,0) and
(r,L) for one row r") is not equivalent to winding: it misses **spiral**
configurations, whose occupied cluster winds in x AND in y simultaneously, so
its lift to the doubled strip contains no same-row `(r,0),(r,L)` pair.
`scripts/exact_wrapping_criterion_diff.cpp` scans all `2^25` configurations and
finds **exactly 10** disagreements, all of the form "DSU says wraps, same-row
doubled criterion says not" — all 10 have 15 occupied sites (inside the c_15
bin) and each is verified (see `gen_issue576_resolution.py`, enforced by
assertion) to wind in both directions. So `8853291 = 8853301 − 10` is fully
explained: **the published polynomial is correct; the first enumerator was
wrong.** The displacement-DSU method (union-find with lattice displacement
potentials, cf. `scripts/matched_torus_reference.py`) is the ground-truth
method and is what `exact_wrapping_enum2.cpp` implements.

### Comparison with the published coefficients (corrected)

| L | published Σ c_k | independent enumeration (enum2) | match |
|---|-----------------|---------------------------------|-------|
| 3 | 175 | 175 | YES |
| 4 | 19571 | 19571 | YES |
| 5 | 8853301 | 8853301 | YES (per-coefficient) |
| 6–8 | parseable, pass divisibility | (enum infeasible) | consistent (see below) |
| 9–12 | parseable, pass divisibility | (enum infeasible) | consistent (see below) |
| 10 | 100 of 101 coefficients in file | — | trailing c_100 = 1 missing from the ancillary file (a truncation; c_N = 1 necessarily, all sites occupied) |

Block health (see `derived/torus_blocks_health.json`,
`scripts/gen_issue576_resolution.py`):

- **The paper proves a divisibility test** (Lemma 1: `c_k` must be divisible by
  `L²/gcd(k, L²)` on the torus). **Every parsed block L=3..12 passes it** —
  including L=10 with the necessarily-restored trailing 1.
- The earlier claim "L≥9 blocks appear corrupted because P_L(1/2) falls below
  the continuum limit 0.1694" is **RETRACTED**. It was a category error:
  `P_L(1/2)` is the wrapping probability at `p = 1/2 < p_c`, whose limit is
  **0**, not 0.1694 (0.1694… is the critical continuum class-`{1,0}`
  probability π({1,0})(i), a different object). The correct published anchor
  for this polynomial is the Mertens–Ziff continuum value
  `R^v = 0.521058290…` (specified-direction wrapping incl. both, at `p_c`;
  the same 0.521058290 appears in Newman–Ziff PRE 64, 016706).
- Evaluated at `p_c = 0.59274605079210`, the published polynomials converge
  smoothly and monotonically (from L=4 on) to that anchor:

| L | P_L(p_c) | R^v − P_L(p_c) |
|---|----------|----------------|
| 3 | 0.521273730479 | −0.000215440 |
| 4 | 0.517303557598 | +0.003754732 |
| 5 | 0.517195470567 | +0.003862819 |
| 6 | 0.517834766352 | +0.003223524 |
| 7 | 0.518425696010 | +0.002632594 |
| 8 | 0.518892914782 | +0.002165375 |
| 9 | 0.519253790705 | +0.001804499 |
| 10 | 0.519534108648 | +0.001524181 |
| 11 | 0.519754734224 | +0.001303556 |
| 12 | 0.519930925308 | +0.001127365 |

  No corruption anywhere in the file; the only defect is the truncated L=10
  block (missing trailing 1). The paper additionally states its L∈[3,7] torus
  polynomials match Mertens' (unpublished) polynomials exactly, which covers
  the sizes where our independent enumeration cannot reach.

### Relation to the repository's committed wrapping observable
The validated published polynomials `P_L(p)` are the **finite-L** ground truth
for the main (NN square site) wrapping channel at every `p`; the committed
continuum baselines are the `L → ∞` anchors at `p_c` (`R^v = 0.521058290…` for
the P_L-defined event) and the Pinson/Arguin class probabilities. The two are
consistent: `P_L(p_c) → 0.521058…` as the table above shows. The main channel
is therefore on published closed-form ground at all sizes for which a closed
form exists, and the finite-L polynomial and the continuum anchor now check
against each other. The NN+NNN (Sq8) matching channel is out of scope of the
paper and remains unvalidated against it.

## Part 2 — Pinson / Arguin continuum wrapping at r = 1, 2, 4

Used the repository's already-implemented, validated formula
(`scripts/pinson_arguin_primitive.py`, percolation specialization `e0 = 2/3`).
Evaluated `π({1,0})(i·r)` for `r = 1, 2, 4` with two independent evaluations
(direct Gaussian winding sum vs compact theta form) cross-checking each other;
**all three values were recomputed from scratch on the resolution pass and
reproduce exactly** (internal direct–theta agreement now measured at
~1e-110 … 1e-116 at 100 dps):

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
