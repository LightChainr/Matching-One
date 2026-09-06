# C5: Markov-generator selection rule (formal statement) and the two numerics #603 did not run

**C5 note.**  Two parts.  (1) The round-2 finite-group selection statement,
rewritten as the formal gap-filling statement #601 Q3 asked for, with citation
discipline: nothing here claims Schur's lemma or Kalman decomposition as new;
the only claimed content is the *Markov-generator / Duhamel-pointwise* form,
its order-`ell` tensor-power extension, and the C2-accident observation.
(2) Two numerical checks that the sibling Phase-C ticket #603 did not run:
the P398 factorization in the **stationary `L2(pi)` inner product**, and the
**isotypic-fibre realisation** on the non-`C2` (C3) object (Theorem 2(c) of
the round-2 note).

Scripts: `scripts/probe/c5_l2pi_fibre.py`; results:
`results/probe-c5-l2pi-fibre/latest.json`.

---

## 1. Formal statement (what the Markov-generator gap actually is)

**Setup.**  `K` a finite group acting on a finite `X` by permutations; `V=R^X`
with the counting product (K-orthogonal); `G` a Markov generator with
`G k = k G`; sources `mu in V_{rho_B}`, readouts `f in V_{rho_C}`; a
perturbation operator `H` transforming in `rho_H` under conjugation.

**Statement (first order).**  If the trivial representation is absent from
`rho_B (x) rho_H (x) rho_C`, then for every `0 <= s <= t` the Duhamel
integrand

```
I(t,s) = < e^{(t-s)G^T} mu , H e^{sG} f >
```

vanishes, hence the first-order linear response of the Markov system is
identically zero.  **Statement (order ell).**  With `ell` insertions of
perturbations transforming in `rho_{H_1}, ..., rho_{H_ell}`, the `ell`-th
order response vanishes when `triv` is absent from
`rho_B (x) rho_{H_1} (x) ... (x) rho_{H_ell} (x) rho_C`.

**Citation discipline.**  The representation-theoretic content is Schur's
lemma (as #601 Q3 concludes); for `C2`-even tasks the quotient balanced-data
statement is the Kalman decomposition of a `C2`-equivariant realisation
(Wonham 1979; #600).  The literature pass of #601 Q3 found **no** statement of
the above *for Markov generators, pointwise in the integrand* (Diaconis 1988
covers Cayley walks only; Wigner–Eckart is Hamiltonian; Antown–Dragičević–
Froyland optimise the reverse direction).  This note therefore claims only:
the explicit Markov-generator form, the pointwise vanishing (stronger than
integral vanishing), the order-`ell` tensor-power criterion, and the proof
chain in the round-2 note.  C2 examples (P398 #598) and a non-C2 example (C3
synthetic) verify it numerically at 1e-17 / 1e-2 scale separation.

**C2-accident reminder.**  `V_triv = Fun(X/K)`, so "project to the invariant
sector" = "take the orbit quotient" only for `K = C2` (or `K` with a single
nontrivial real irrep).  For general `K` the quotient and the isotypic
decomposition are different objects, and exposing a nontrivial fibre requires
a character-marked channel.  (Round-2 note, Section 3.)

## 2. Numerics #603 did not run: L2(pi) factorization on P398

#603 verified Phase C under the repository's frozen counting-measure
convention.  The question left open by #598's decision table was whether the
factorization survives a different **inner product**.  The stationary measure
`pi` of `G0` is reflection-invariant (hence orbit-constant), so the `L2(pi)`
product is another `K`-compatible convention.  Full-space vs orbit-quotient
task spectra in `L2(pi)`, widths 4–8, **protected (all R-even) readouts**:

| w | L2(pi) rel diff full vs quotient | counting-L2 rel diff (round 1) |
|---|-----------------------------------|--------------------------------|
| 4 | 2.7e-16 | 5.2e-16 |
| 5 | 1.5e-16 | 2.6e-16 |
| 6 | 2.4e-16 | 2.4e-16 |
| 7 | 1.8e-16 | 2.3e-16 |
| 8 | 2.5e-16 | 2.3e-16 |

**Convention dependence of the numbers.**  The spectra themselves differ
between the two products (w=8 leading singular value 0.1197 in `L2(pi)` vs
0.6406 in counting), yet the *factorization* holds in both to ~1e-16.  So:
the balanced data are convention-dependent; the quotient identity is
convention-robust.  This is the cleanest form of the #598 "inner-product
diagnosis": a wrong convention changes the numbers (and breaks a mismatched
quotient implementation by tens of percent, round 1), but the theorem's
content does not depend on which compatible convention is chosen.

**Teeth.**  Including the R-odd readout `halves_linked` (odd widths) breaks
the comparison (w5: 0.12, w7: 0.03 in `L2(pi)`, exactly as in counting),
because the quotient cannot represent its odd part; removing exactly the odd
part restores agreement.  This reproduces #603's control under the second
inner product.

## 3. Numerics #603 did not run: isotypic-fibre realisation (C3)

On the six-state C3 synthetic object (`V = 2 triv + 2 W`), the W-sector task
(source and readout in W) was computed two ways:

1. full 6-state space;
2. the **isotypic fibre system**: `G` restricted to the 4-dimensional
   `W`-isotypic subspace, `A_W = W^T G W`.

Task singular spectra agree to **1.4e-15** (fibre [0.0446, 0.0114, 0.0009, …]
= full).  The joint two-sector task (triv + W channels) has a spectrum equal
to the union of the per-sector spectra (block-diagonality of `e^{tG}` on
isotypic sectors).  This is Theorem 2(c) of the round-2 note made numerical:
for a nontrivial sector `rho != triv` the correct "reduced system" is the
fibre system, not a Markov chain on `X/K`, and constructing it requires
choosing a marked basis in `rho`.

## 4. Boundary and position

* Both numerics are deterministic, no sampling, no new width, no fit.
* Statements are exact finite statements on P398 / the declared C3 object.
* Relative to #603: this note adds the second-inner-product robustness check
  and the non-C2 fibre check; it does not rerun the frozen-convention Phase C.
* The order-`ell` selection rule remains the single most reusable exact
  output of the probe for discriminating "genuine symmetry-breaking response"
  from "broken equivariance / misdeclared character" in future repository
  perturbations.

Related: #598, #600, #601, #603; round-2 note
`probe-finite-group-response-isotypic-factorization-20260906.md`.
