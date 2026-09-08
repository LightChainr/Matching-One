# Issue 112 retrieval slice: Riordan–Walters / BBW finite criterion vs our exact wrapping/M tables

Status: Phase 0 feasibility, retrieval + primary-source quotes only. No new Monte Carlo, no production
implementation, no new bound claimed. Companion to `notes/rigorous-pc-confidence-feasibility.md`
(PR #278), `notes/pc-block-event-exact.md` (#279), `notes/pc-block-event-linear-evaluator.md` (#280),
`notes/pc-confidence-preregistration.md` (#281). This slice answers the five questions in the issue
body against the repo's actual data holdings, with verbatim primary-source quotes.

Verdict: **methodologically interesting but weak** (unchanged from the #278 classification).
The published criterion applies directly, the certificate machinery already exists on main, but
the repo's exact wrapping/M tables cannot feed it, and the realistic output is a rigorous
confidence interval, not a deterministic decimal theorem.

## Primary sources (all quotes verbatim)

1. O. Riordan, M. Walters, "Rigorous confidence intervals for critical probabilities",
   Phys. Rev. E 76, 011110 (2007), arXiv:math/0702232. [Checked 2026-09-08, arXiv PDF.]
2. P. Balister, B. Bollobás, M. Walters, "Critical Probabilities of 1-Independent Percolation
   Models", Combin. Probab. Comput. 21(1-2), 11–22 (2012), doi:10.1017/S0963548311000538.
   [Bibliographic data checked 2026-09-08 via Cambridge Core; abstract quoted.]
3. P. Balister, T. Johnston, M. Savery, A. Scott, "Improved bounds for 1-independent percolation
   on Z^n", arXiv:2206.12335; published Electron. J. Probab. (2025), doi:10.1214/25-EJP1341.
   [Checked 2026-09-08, arXiv HTML v2.]

## Q1. Which published finite-size criterion applies directly to square-site percolation on Z^2?

**Riordan–Walters 2007, which is the BBW 1-independent renormalization reduction.** It applies to
site percolation on the square lattice without modification. The deterministic core is their
Lemma 1 (RW 2007, p. 3):

> "Let P~ be a 1-independent bond percolation measure on Z2 in which each bond is open with
> probability at least p0 = 0.8639. Then the probability that the origin lies in an infinite
> open cluster is positive."

RW attribute this lemma to Balister–Bollobás–Walters [ref 2 above] and note a weaker form follows
from first principles or from Liggett–Schonmann–Stacey. The renormalization statement (RW 2007,
pp. 3–4) fixes the route from site percolation to a 1-independent bond model:

> "given a 'scale parameter' s > 0, partition R2 into disjoint s by s squares Sv, v ∈ Z2. For each
> bond e of Z2 let Re be the corresponding rectangle ... Let Ee be some event that depends only on
> the states of the sites or bonds of Λ that lie within Re, and take the bond e of Z2 to be open
> with respect to P~ if and only if Ee holds. Since the rectangles corresponding to vertex-disjoint
> bonds of Z2 are disjoint, this defines a 1-independent measure P~."

The modernized deterministic constant replaces 0.8639 with the BJSS 2022 Theorem 1.3:

> "Theorem 1.3. p_max(Z^2) ≤ 0.8457."

and their usable certificate form (BJSS, discussion of the renormalized block event E_uv):

> "If P(E_uv) > 0.8457 for all edges uv, then by Corollary 5.1 below there are almost surely cycles
> enclosing any bounded region in the [original] model."

So: criterion = BBW 1-independent lemma + RW two-cell block event; strict marginal threshold
0.8457 (was 0.8639 in 2007).

## Q2. What finite rectangle/annulus/circuit events and numerical thresholds are required?

The block event on the open-boundary `2s x s` rectangle `R_e = S_u ∪ S_v` (RW 2007, p. 4):

> "Let Ee be the event that each of Λu and Λv contains a unique largest open cluster with these
> clusters part of the same open cluster in Λe. Here 'largest' simply means containing the most
> sites."

Requirements, all already frozen on main by PRs #279/#280:

- open-boundary (not periodic) `s x s` halves and `2s x s` union; frozen ordering:
  empty-half rejection, unique-largest-tie rejection, half-local selection, union connectivity;
- two graphs: square NN (upper side) and square site-matching NN+diagonals (lower side, via
  `p_c(square) >= 1 - q`);
- threshold: strict `P_p(E_e) > 0.8457` for every block edge, certified by an exact binomial tail
  at preallocated per-run alpha `10^-6/6 = 1/6000000` (family-wise `10^-6`, two sides, ≤3 attempts
  each — the 2007 allocation, RW 2007 p. 4: "we can perform up to three different runs with
  different parameters s and p ... so the probability that our final bound is incorrect is at most
  10−6/2");
- no annulus/circuit event is required for the basic certificate (BJSS's Corollary 5.1 provides
  the cycle/uniqueness consequences; they are outputs, not inputs).

Exact reproduction of the 2007 statistical step (RW 2007, p. 4):

> "If Pp(Ee) = π, then m has a binomial Bi(N, π) distribution with parameters N and π. In
> particular, if π < 0.8639, then P(m ≥ 378) ≤ P(Bi(400, 0.8639) ≥ 378) = 1.1489··· × 10−7 <
> 10−6/6."

Our `scripts/rigorous_pc_confidence_gate.py` reproduces this tail exactly (Fraction arithmetic):
legacy cutoff 378, tail `1.14899035289400922e-7`; modern cutoff at `p0 = 8457/10000` is 373,
tail `9.51451513592604452e-8`; 372 fails (`2.37086590688800581e-7`).

## Q3. Can our existing exact/M engines estimate those events with a mathematically valid error bound?

**The wrapping/M tables themselves: no — wrong observable, wrong boundary.** Inventory:

- `results/exact_small_matching_polynomials.md` + `scripts/exact_matching_polynomial.py`:
  exact Mertens–Ziff matching function M(p) as an integer polynomial, axis torus L=2..4
  (N = 4, 9, 16 sites) and diamond torus L=2,3. These are *periodic torus wrapping* observables
  (black NN wraps minus white matching wraps), not open-boundary unique-largest-cluster events.
  A torus wrapping probability is a global topology statement; `P_p(E_e)` is local with free
  boundaries. No algebraic identity connects the two families, so the tables cannot substitute
  for block-event probabilities at any s.
- `data/jacobsen_2015_square_site_cylinder.csv` (n=1..21) and
  `data/mertens_2022_square_site_estimators.csv` (p_med n=1..24, p_cell n=2..24): literature
  finite-size *estimator* sequences for threshold extrapolation. They are not event probabilities
  and cannot certify `P_p(E_e) > 0.8457`; at best they inform (heuristically) where to place p
  and s — the same role RW 2007 gave their own heuristic extrapolation.
- What *can* estimate the required event with a valid error bound already exists on main:
  - exact enumeration of all `2^(2s^2)` configurations with `Fraction`-exact probabilities
    (`scripts/pc_block_event_exact.py`, s=1,2 exhausted; coefficient vectors frozen, e.g.
    matching s=2 `[0,0,4,20,41,44,26,8,1]`, `P_{1/2}(E) = 9/16`);
  - an O(s^2) deterministic evaluator (`scripts/pc_block_event_linear.py`), differentially
    validated on all 520 s∈{1,2} graph/configuration pairs;
  - an exact-`Fraction` binomial-tail oracle with predeclared family-wise error and a
    fail-closed pre-registration ledger (`scripts/rigorous_pc_confidence_gate.py`,
    `scripts/pc_confidence_preregistration.py`, plan SHA-256 `8beeb397...`).

Caveat that carries over from #279: the square-NN block event is **not increasing** (28
success-to-failure one-site openings at s=2; frozen witness mask 6 → 22 `left_tie`). This does
not break the fixed-parameter 1-independence + binomial certificate, but it forbids using the
event as an increasing observable in coupling arguments, and it means polynomial-root/heuristic
placement of (s, p) has one more failure mode than for monotone crossing observables.

## Q4. What lattice size and sample count for a nontrivial improvement over textbook bounds?

Two quantitative ceilings, one statistical and one geometric:

- Statistical power is not the binding constraint. At N=400, cutoff 373, power is ≈0.015 / 0.206
  / 0.774 / 0.952 at true event probability 0.90 / 0.92 / 0.94 / 0.95. Larger N is cheap; the
  design requirement is that the block geometry already push `P_p(E_e)` into the mid-0.9 range.
- Geometry is binding. RW 2007 published square site
  `[0.5925, 0.5930]` (Table I: "Square [0.5925,0.5930] 5 × 10−4"), and the abstract claims
  "we obtain intervals of width at most 0.0005 in all cases". Under their (explicitly heuristic,
  conformal-invariance-motivated) planning law `δp ~ s^(-3/4)`, halving width to 2.5e-4 costs
  ≈2.5× linear scale / ≈6.4× area; width 1e-4 costs ≈8.6× linear / ≈73× area (table in
  `notes/rigorous-pc-confidence-feasibility.md`). Per-trial cost is O(s^2) with the #280
  evaluator, so the compute is modest in absolute terms — but any produced interval, however
  narrow, is dominated by ~9-10-digit published estimates (Jacobsen 2015, Mertens 2022, already
  canonical in `data/`). The theorem-level statement's value is auditability, not precision.
  Against *textbook* bounds (e.g. the Ermakov–van den Berg style lower bounds), even a modest
  rerun improves on nothing that is not already superseded by RW 2007's own published interval.

Conclusion: a rerun with s of a few hundred and N ≈ 400–2000 per attempt is computationally easy
(seconds-to-minutes per trial at O(s^2) per trial on a laptop for s ≤ 512) and would merely
reproduce-or-slightly-tighten a 2007 result. Nontrivial *scientific* improvement is not on offer.

## Q5. Is a deterministic exact enumeration/transfer-matrix variant feasible?

**No. The result is necessarily a rigorous confidence interval, not a deterministic theorem.**

- Exact enumeration sums over `2^(2s^2)` configurations: s=2 is 65,536 (done), s=3 is 2^18
  ≈ 2.6e5 × per-config union-find — already past this project's Mac compute budget per the
  #653 experience, and the polynomial coefficients blow up far faster than the count. RW 2007
  themselves dismiss this route in one sentence: "For sufficiently small scale parameters s, it
  is possible to find a p for which (2) holds by enumerating all possibilities for which
  sites/bonds in Re are open, and so writing Pp(Re) as a polynomial in p. Needless to say, this
  is impractical and gives poor results in practice."
- A transfer-matrix frontier method is exponential in the width s (the #279/#280 notes and the
  RW strip algorithm both reduce *memory* to O(s), not the enumeration count). The scales where
  the event probability clears 0.8457 with margin (s in the hundreds, since the event must
  succeed at p ≈ 0.5927 ≪ 0.8457 only because s is huge) are unreachable exactly.
- Even the statistical certificate has an irreducible non-mathematical assumption. RW 2007, p. 5:
  "...our procedure has probability at least 99.9999% of producing an interval containing the
  true value (assuming the random number generator we used is well behaved)." The output is a
  high-confidence interval conditional on fresh genuine IID Bernoulli trials — exactly the
  boundary stated in the #278/#281 notes.

## Outcome classification

**Methodologically interesting but weak.** Specifics:

- Go-components: the criterion applies verbatim; the certificate oracle, tiny-s exact semantics,
  linear evaluator, and pre-registration ledger are merged and validated; the modern constant
  0.8457 lowers the N=400 cutoff from 378 to 373 and is already wired in.
- Weak-components: the repo's exact wrapping/M tables (torus, periodic, L≤4; literature
  estimator tables n≤24) are the wrong observable family and contribute nothing to the
  certificate; the block event is nonmonotone in the square graph; exact enumeration dies at
  s=3; and any achievable interval is scientifically dominated by existing 9-10-digit
  estimates. The only unique deliverable — a fully auditable finite-computation statement —
  was already delivered at width 5e-4 by RW 2007, and our rerun would add no new mathematical
  content.

Recommendation: keep #112 open as parked (per the 2026-08-31 weekly decision; do not prepend
handoff blocks until `docs/STATUS.md` reprioritizes). If the project ever wants the
theorem-level side result, the remaining unmerged work is exactly the #278 production gate:
exploratory (s, p) placement using the heuristic planning law, then one fresh pre-registered
certification run. No server time is justified now.

## Gap found and fixed in this slice

`references.bib` had no entries for the three primary sources (they existed only as inline
arXiv IDs in notes). Added: `RiordanWalters2007RigorousCI`,
`BalisterBollobasWalters2012OneIndependent`, `BalisterJohnstonSaveryScott2022Improved`.
