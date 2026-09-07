# Probe baseline snapshot v2 — repository front moved (C1)

**Date:** 2026-09-06 (afternoon).  **Why:** the repository front moved
~1.5 h after round 1/2: three sibling agents produced #600/#601/#602/#603
that overlap round 1/2.  Before continuing the probe, every claim of
rounds 1–2 was checked against the new front.

## New activity map (2026-09-06, 12:07–12:58 UTC)

| item | actor | content | relation to this probe |
|---|---|---|---|
| #600 | DeepSeek ticket | Phase C compute job for #598; frozen-convention balanced build; names Wonham/Kalman as the expected theorem | overlaps round 1 |
| #603 (PR) | DeepSeek result | Phase C numbers: protected D0 A≡B bit-for-bit (1.5e-16..2.7e-15), odd sector inert (<6e-17), odd-width `halves_linked` teeth (1.8e-2, 8.7e-3), symmetrization restores 1e-16; continuity gate vs committed #588 numbers | **independent confirmation of round-1 Phase C** |
| #601 | Grok literature pass | Q1 Callan–Smiley/Ding: fixed-NC counts and orbit formula are known (OEIS A007123); Q2 Gate-1 "coarsest = orbits" is a conditional theorem, P398 is a 5-width instance (D'Angeli–Donno counterexample); **Q3: the Markov-generator pointwise selection-rule statement is a REAL literature gap**; Q4/Q5 method notes | sharpens the positioning of round 2 |
| #602 (PR) | Grok notes | literature officer notes (2 rounds) | see #601 |

## Conflict / overlap check against rounds 1–2

1. **Round-1 Phase C numbers.**  #603 reproduces them under the repository's
   own frozen convention (its own implementation), with the same structure:
   full-space ≡ quotient to machine precision on the even dictionary; odd
   sector exactly inert; odd-width `halves_linked` breaks and symmetrization
   restores.  **No conflict — mutual independent confirmation.**  Round-1's
   value is now that of a second independent implementation under a different
   (declared) convention; #603's value is the frozen-convention continuity
   gate.  Both agree; the factorization claim no longer rests on one script.
2. **Round-2 representation theorem.**  #601 Q3 searched for precisely this
   statement ("invariant observable, Markov generator, symmetry-breaking
   perturbation ⇒ first-order response identically zero, pointwise in the
   integrand") and reports it is **not found in the literature** — Diaconis
   1988 covers group (Cayley) walks only; Wigner–Eckart is Hamiltonian;
   Markov linear-response optimization (Antown et al.) runs the opposite way.
   The gap is declared real.  Round-2 therefore occupies exactly the slot #601
   Q3 marks empty, with two caveats to fix in that note:
     * the *representation-theoretic criterion itself* is Schur's lemma —
       round 2 must not claim that as new; the contribution is the explicit
       Markov-generator/Duhamel pointwise form, the full-order version
       (`rho_H^{x ell}`), and the proof chain;
     * Theorem 2(b) (quotient balanced data) is, for `C2`-even tasks, the
       Kalman decomposition of a `C2`-equivariant realization (Wonham) as
       #600 states — round 2 should cite that and keep only the genuinely
       additive parts (isotypic-vs-orbit correction, fibre reading of
       nontrivial sectors, inner-product diagnosis).
3. **Orbit-count arithmetic.**  #601 Q1 (Callan–Smiley 2005; Ding 2016)
   settles that the fixed-state count `C(w, floor(w/2))` and the orbit formula
   are known combinatorics (OEIS A007123); round-1's `verify_gate1` is a
   numerical check of it, and the partition-equality (Gate-1 branch A) stays a
   P398 5-width *instance*, not a theorem (Q2, D'Angeli–Donno counterexample).
   Round-1 already carried that caveat; no change needed beyond citing.
4. **experiment-language complexity theory** (the probe's north-star):
   confirmed **completely absent** from the whole front — no issue, note, or
   PR touches it.  The probe's atlas/ladder/rank-vs-classes programme remains
   the only occupant of that territory.

## Plan corrections (applied to PLAN-12-turns.md)

* C5 is re-scoped: do **not** repeat frozen-convention Phase C (done by
  #603).  New C5 = markov-generator selection-rule formal statement aligned
  with #601 Q3 (citation-disciplined), + isotypic-fibre numerics that #603 did
  not do (non-`C2` verification of Theorem 2(c)).
* C4 gains priority: no competitor occupies the strict-separation /
  predictive-class-vs-rank territory; it is the probe's most differentiated
  open workstream.
* C1's own product is this note + the PROGRESS update; no further C1 action.
