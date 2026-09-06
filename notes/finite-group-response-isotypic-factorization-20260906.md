# Finite-group response selection and isotypic factorization: one theorem for #244 and #598, and why "even = orbit-constant" is a C2 accident

**Round 2 probe note.**  This note upgrades the P398 parity result (#598,
which is the finite `C2` case) and the deck-character rule (#244) into one
statement about an arbitrary finite group `K` acting on a finite Markov
state space, with a complete proof, a numerically verified non-`C2` example
(`K = C3`), and a structural correction that matters for how the P398
factorization may be generalized.

**Scope / claim boundary.**  Exact finite statements only.  `K` labels are not
continuum spins; a finite-group selection rule is not a CFT statement; the
balanced order stays task-relative.  The synthetic `C3` object is a declared
test object, not a repository system; it is used only to show the theorem is
not an artefact of `C2`.

---

## 1. Setup

Let `X` be a finite state space and `V = R^X` the space of functions on it,
with the counting inner product `<a,b> = sum_x a_x b_x`.  A finite group `K`
acts on `X` by permutations, hence orthogonally on `V` (`K` is an isometry of
the counting product).  Write the isotypic decomposition

```
V = oplus_rho V_rho ,   V_rho = m_rho copies of the (real) irrep rho.
```

Isotypic subspaces are mutually orthogonal; `V_triv = Fun(X/K)` is exactly the
space of `K`-invariant functions, canonically the function space of the orbit
quotient.

Let `G` be a continuous-time Markov generator (`Q`-matrix) on `X` that is
`K`-equivariant: `G k = k G` for every `k in K` (equivalently: transition
rates are constant on `K`-orbits of pairs).  Then `e^{tG}` and `e^{tG^T}`
preserve every isotypic subspace.

A *task* is a set of source measures `mu_a` (row vectors, but we identify
`R^X` with its dual through the counting inner product, where the `K`-action
is orthogonal and self-dual) and readout functions `f_j`, with response kernel

```
K_{a,j}(t) = <mu_a, e^{tG} f_j>.
```

A perturbation is an operator `H` on `V`.  `K` acts on operators by
conjugation `H -> k H k^{-1}`; we say `H` *transforms in the representation
`rho_H`* if it lies in the `rho_H`-isotypic part of the conjugation
representation on `End(V)`.  Concretely, in every repository example `H` is a
"field" built from functions/measures that carry `rho_H`.

## 2. Theorem 1 — response selection rule (all orders, any finite group)

**Setup.**  Fix a source `mu in V_{rho_B}`, a readout `f in V_{rho_C}`, and a
perturbation `H` transforming in `rho_H`.  The linear-response integrand is

```
I(t,s) = < e^{(t-s)G^T} mu ,  H e^{sG} f >,     0 <= s <= t.
```

**Theorem 1a (first order).**  If the trivial representation does **not**
occur in `rho_B (x) rho_H (x) rho_C` — equivalently `rho_B` does not occur in
`rho_H (x) rho_C` — then `I(t,s) = 0` for **every** `0 <= s <= t`; hence the
first-order response `d/de <mu, e^{t(G + e H)} f>|_0 = int_0^t I(t,s) ds`
vanishes identically.

*Proof.*  `K`-equivariance of `G` gives `e^{sG} f in V_{rho_C}` and
`e^{(t-s)G^T} mu in V_{rho_B}`.  Applying `H` (which transforms in `rho_H`)
to `V_{rho_C}` lands in the isotypic subspace of `rho_H (x) rho_C`.  If
`rho_B` has zero multiplicity in that tensor product, the vector
`H e^{sG} f` is orthogonal to the whole subspace `V_{rho_B}` (isotypic
subspaces of distinct irreps are orthogonal under any `K`-invariant inner
product, in particular the counting product), so `I(t,s)=0`.  The claim is
therefore pointwise in the time integrand — a quadrature that happened to
cancel would be a strictly weaker statement.  `QED`

**Theorem 1b (order ell).**  More generally, a term of order `ell` in `H`
(an `ell`-fold insertion `H ... H` at ordered times, or `ell` perturbations
transforming in `rho_{H_1}, ..., rho_{H_ell}`) vanishes identically whenever
the trivial representation is absent from
`rho_B (x) rho_{H_1} (x) ... (x) rho_{H_ell} (x) rho_C`.

*Proof.*  Same argument along the chain
`V_{rho_C} -> V_{rho_{H_ell} (x) rho_C} -> ... -> V_{rho_H^{x ell} (x) rho_C}`
and the orthogonality of isotypic subspaces.  `QED`

**Corollaries.**

* `K = C2`, `rho_B = rho_C = triv`, `rho_H = det` (the odd character): the
  product is `det`, which does not contain `triv`; first-order response is
  zero.  This is #598 Phase A, now as a one-line corollary.  (Second order:
  `det (x) det = triv`, so the even sector can be reached at order 2 — #598
  Phase B.)
* #244's deck characters: an invariant observer against a perturbation
  carrying a nontrivial deck character `chi` has `rho_B = rho_C = triv`,
  `rho_H = chi`, product `chi`, no `triv` — zero linear response.  The two
  repository rules are the same theorem with the same proof.

**Numerical verification outside C2 (`K = C3`).**  The synthetic object in
`scripts/probe/group_selection_demo.py` has six states (two `C3`-orbits),
a rotation-equivariant generator, `V = 2 triv oplus 2 W` (`W` the real 2-d
irrep of `C3`), and a diagonal perturbation field transforming in `W`.
Measured peak pointwise integrands over a `(t,s)` grid:

```
source   readout   perturbation   peak |integrand|   prediction
triv     triv      W              2.6e-17            ZERO   (triv not in W (x) triv)
triv     W         W              1.34e-2            nonzero (triv in W (x) W)
W        triv      W              7.10e-2            nonzero
W        W         W              1.54e-2            nonzero
```

The zero row is at machine precision, exactly as in #598's `C2` tables; the
other rows show the rule is selective, not vacuous.  The second-order
re-entry of the invisible direction into the *unmarked* task (`triv`/`triv`)
is also measured: log–log slope **2.006** over `eps in {0.02, 0.04, 0.08}`
(predicted 2, since `W (x) W = triv oplus W`), with the control — a
`triv`-transforming perturbation on the same task — at slope **1.022**.
Data: `results/probe-group-selection/latest.json`.

## 3. The structural correction: "even = orbit-constant" is special to `C2`

#598 Phase C factorizes P398's task through the *reflection quotient*, using
the identification "R-even functions = functions constant on R-orbits".  That
identification is true for `C2` **because** `C2` has a single nontrivial
irrep and the even subspace is the invariant subspace.  For a general finite
group it fails in a way that matters:

```
V_triv = Fun(X/K)   (orbit quotient; what an unmarked observer sees)
V_rho  = nontrivial isotypic fibres (invisible to an unmarked observer)
```

`V_rho` is *not* a function space on `X/K`.  An isotypic fibre is a vector
bundle over the orbit space (a `K`-equivariant bundle with typical fibre
`rho`); its sections are exactly the functions that transform in `rho`, and
two states in the same `K`-orbit can differ in `V_rho` (e.g. the two points of
a `C2` orbit carry opposite odd values, or the three points of a `C3` orbit
carry the three phases).  Consequently:

1. **The correct general first factor is the isotypic decomposition, not the
   orbit quotient.**  #598 Phase C's quotient factorization is the `C2`
   special case in which "project to `V_triv`" and "average over the group"
   coincide.
2. **An unmarked observer cannot even in principle resolve `V_rho`.**
   If every source and readout is `K`-invariant, all response data factor
   through `Fun(X/K)` — the non-trivial fibres are *unobservable in the sense
   of identifiability*, not merely hard to measure.  This is the exact,
   finite-group content of "complexity is indexed by the experiment language
   together with its symmetry": two states in one `K`-orbit are the same
   object to every symmetric task.
3. **Exposing a fibre requires a marked channel.**  Theorem 1a says precisely
   when a perturbation `H` in `rho_H` is first-order visible to a task with
   sources in `rho_B` and readouts in `rho_C`: exactly when `triv` occurs in
   `rho_B (x) rho_H (x) rho_C`.  For an unmarked task (`rho_B = rho_C = triv`)
   this requires `rho_H = triv`: the *whole* `rho_H != triv` direction is
   invisible to first order, matching #598's "marked readout" row — with the
   sharpening from round 1 that one marked *readout* alone is still not
   enough unless the perturbation flips it back into an even subspace (in
   `C2`: `H_odd e^{sG} f_odd` is even, so an odd *readout* plus an odd
   *perturbation* is required; equivalently a marked source).  The theorem
   makes the exact counting explicit: visibility ⇔ a representation-theoretic
   tensor-product condition, not a heuristic about "breaking symmetry".

## 4. Theorem 2 — isotypic factorization of the balanced task data

Let the task have all sources and readouts in one isotypic sector `rho` (this
includes the unmarked case `rho = triv`), and let `<>` be any `K`-invariant
inner product for which isotypic subspaces are orthogonal (counting measure
and the stationary `L2(pi)` product both qualify, since `pi` is `K`-invariant
when `G` is `K`-equivariant and the chain is irreducible).  Then:

**(a) sector decomposition of the response.**  `e^{tG}` preserves `V_rho`;
the response kernel of the task is carried entirely by `V_rho`, and the task
Hankel/Gramian data of the full system equals the data of the system
restricted to `V_rho`, with all other sectors contributing exactly zero.

**(b) quotient case (`rho = triv`).**  `V_triv = Fun(X/K)` and the orbit
partition is a strong lumping (K-equivariance of the rates), so the restricted
system on `V_triv` is exactly the Markov chain on the orbit quotient `X/K`.
The full-space and quotient balanced singular data coincide — provided the
quotient is given the induced (orbit-weighted) inner product.  This is the
round-1 verification on P398 widths 4–8 (agreement to ~1e-16) as a special
case.

**(c) nontrivial case.**  For `rho != triv` the correct "quotient" is not a
Markov chain on `X/K` but the isotypic fibre system; the balanced data equals
the data of that fibre system.  Constructing it requires choosing a marked
(basis) channel in `rho` — which is exactly the operation of *adding an
observer with nontrivial character*.

*Proof sketch.*  (a) is the block-diagonality of `e^{tG}` and of the
source/readout matrices on isotypic subspaces, plus orthogonality.
(b) is the strong-lumping fact + the induced inner product (round-1 note,
Section 3).  (c) is (a) read in fibre coordinates.  `QED`

## 5. What this buys the repository

* **#244 and #598 are one theorem.**  Deck characters (`#244`) and reflection
  parity (`#598`) are both instances of Theorem 1a; P398's phases A/B are the
  `C2` corollary with exact finite verification, and the new `C3` demo shows
  the statement is not `C2`-specific.
* **The factorization has a correct generalization target.**  "Symmetry
  quotient, then task compression" (#599 high-value item 4) should be stated,
  in general, as *isotypic decomposition, then per-sector task compression*;
  the orbit quotient is the `triv` sector and coincides with "even functions"
  only for `C2`.  Any attempt to generalize P398's Phase C to a larger
  symmetry group via orbit functions alone would silently discard the
  non-trivial fibres.
* **A cheap prospective discriminator.**  Theorem 1a gives, for any declared
  finite symmetry `K` of a microscopic family, an exact zero-test for candidate
  "symmetry-breaking" perturbations before any sampling: compute the
  multiplicities in `rho_B (x) rho_H (x) rho_C`.  If zero, the direction is
  first-order invisible *by representation theory*, and any measured
  first-order signal must be blamed on a broken equivariance or a misdeclared
  character — the same logic #598 used to find `halves_linked`.

## 6. Limits and open ends

* The theorem is about linear (and order-`ell`) response of *finite* systems.
  Nothing here upgrades "invisible at first order" to "irrelevant at finite
  coupling": at order `ell` the criterion is the tensor power
  `rho_H^{x ell}`, and for `rho_H != triv` the *first* `ell` with
  `rho_H^{x ell} contains triv` gives the order at which the hidden sector can
  re-enter an unmarked task.  For `C2`, `ell = 2` (#598 Phase B measured slope
  2.0006–2.0010).  For `C3` with `rho_H = W`, `W^{x 2} = triv oplus W`, so
  `ell = 2` again; for `rho_H` of order `d` in the representation ring the
  re-entry order is the order of `rho_H` in the subgroup of the character
  group.  Computing that re-entry order on a repository family with a genuine
  cyclic deck is the natural next exact check.
* Balanced-data claims hold for `K`-invariant inner products.  A
  non-`K`-invariant convention (round 1 showed the plain Euclidean product on
  the quotient changes P398 spectra by 37–61%) breaks the statement, and that
  breakage is diagnostic, not a counterexample to the theorem.

## Files

`scripts/probe/group_selection_demo.py` — C3 demo (deterministic);
`results/probe-group-selection/latest.json` — its numbers;
round-1 note `probe-p398-quotient-factored-balancing-20260906.md` — the C2
verification on widths 4–8.
