# C8: catalogue of what survives observer changes and what does not

**Status:** systematic catalogue over the probe's verified material
(Program G / #599 direction 7).  Each row gives the object, its behaviour
under the named change, the evidence, and where it lives.  "Invariant" here
always means *provably unchanged for the declared family*; several entries
are explicitly finite-object facts, not conjectures.

Legend: **=** invariant (verified exact/machine); **~** changes quantitatively
but some structure is preserved; **!** can change arbitrarily; **?** not yet
determined.

## 1. Under change of inner product (counting vs L2(pi))

| object | behaviour | evidence |
|---|---|---|
| quotient = full task spectrum identity (even task, P398) | **=** (both products, ~1e-16, w4–8) | round-1 + C5 |
| values of the task singular spectra | **!** (w8 leading sv: 0.6406 counting vs 0.1197 L2(pi)) | round-1 / C5 |
| balanced order at a frozen energy tolerance (rough order) | **~** | C5 heads still decay fast in both products |
| odd-sector inertia (zero response of even tasks to odd readouts) | **=** | round-2 / C5 teeth (restriction breaks, symmetrisation restores) |
| "the quotient inherits orbit weights" | **=** structure; using plain Euclidean on the quotient **!** breaks the identity (37–61%) | round-1 |

## 2. Under change of declared readout dictionary (P398)

| object | behaviour | evidence |
|---|---|---|
| D0/D1/D2 static signature class counts | **!** grows 32 → 156 → 209 (w8); 10 = orbits only at w4 | C4 note, results JSON |
| C2-orbit partition as the coarsest admissible strong lumping for J, D (D0-admissible) | **=** for D0; adding a non-R-invariant readout **!** destroys it (identity at w5+ with all 8 readouts forced? verified for non-even demand) | round-1 gate1; #598 comment |
| parity classification of readouts (halves_linked odd only at odd w) | **=** structural | verify_gate1 |
| balanced task spectra of an even task | **~** change when dictionary changes; full = quotient in each declared dictionary | phase_c runs |

## 3. Under symmetry sector surgery (C7 two-copy family)

| object | behaviour | evidence |
|---|---|---|
| every even-task datum under change of the odd-sector block (g in [0,1]) | **=** to 1e-16 | C7 |
| full-chain spectral gap / ergodicity | **!** (g·m vs m; ergodicity loss at g=0) | C7 |
| any claim that task order constrains the odd sector | **false** (no-go) | C7 |

## 4. Under symmetry language (K-invariance / characters)

| object | behaviour | evidence |
|---|---|---|
| state separation by K-invariant tests | **!**: can never separate states in one K-orbit | round-2 Thm 1/2; F1 of C4 |
| first-order response of a task (mu in rho_B, f in rho_C) to a rho_H perturbation | **= 0** iff triv absent from rho_B⊗rho_H⊗rho_C (any K; pointwise) | round-2 Thm 1a/1b; C5 formal statement; C3 demo; P398 #598 |
| C2 even = orbit-constant identification | holds for C2 only; for general K trivial sector = Fun(X/K), other isotypic fibres need marked channels | round-2 §3 (C2 accident) |

## 5. Under language depth/roots (#549 fork family)

| object | behaviour | evidence |
|---|---|---|
| single-root fork success probability as function of class coordinate a | affine (degree 1); rank <= 2 for any finite single-root language | C2 proposition + numerics |
| d-root (depth-d) response rank | degree d; tests d=0..k give exact rank k+1 | C2 |
| class count kappa vs rank r | **!**: kappa = k+1 with r = 2 (single-root), kappa = r = k+1 (depth k) | C2/C4 |
| nonnegative rank of the affine response matrix | = ordinary rank = 2 (columns are the nonnegative factor) | C4 fact 3 |

## 6. Objects that are invariant in the strongest sense seen so far

1. **Exact selection zeros**: the parity/tensor-product vanishing of Theorem
   1a/1b (holds for every K-compatible inner product, every width).
2. **Quotient-factorization identity** (full task spectrum = quotient task
   spectrum) for any K-compatible inner product (verified counting and
   L2(pi)).
3. **Sector inertness** of even tasks against odd channels.
4. **Protocol algebra of the #549 language** (affine-in-a single-root form;
   degree = depth), which is a function of the language, not of the state
   space.

## 7. What this means (one paragraph)

For P398 the only "state-dimension-like" numbers that survive inner-product
and symmetry-bookkeeping changes are **selection zeros**, the **quotient
identity**, and **per-sector inertness**; every positive magnitude (spectrum,
order, class count, rank) is either dictionary-dependent or inner-product-
dependent or language-dependent, and each of them is attached to a *declared*
object.  The catalogue is the raw material for the complexity-ladder document
(C9): an arrow labelled "unconditional" in that ladder has to pass through
this table without breaking.

Related notes: round-1 (`probe-p398-quotient-factored-balancing`),
round-2 (`probe-finite-group-response-isotypic-factorization`),
C4 (`probe-experiment-languages-formal`), C5
(`probe-markov-selection-rule-and-l2pi-fibre`), C7
(`probe-threshold-no-go-theorem`).
