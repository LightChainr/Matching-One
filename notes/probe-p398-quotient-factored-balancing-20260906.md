# P398: the task balanced realization factors exactly through the reflection quotient

**Status:** exact finite computation + structural theorem (widths 4–8).
**Closes the open Phase C of #598** and supplies the inner-product diagnosis that
#598's decision table flagged as the thing to check before quoting any balanced
order.
**Scope:** P398 finite connectivity process only. Not percolation, not a field
statement, not a claim that balanced order 3–4 is a physical state dimension.

---

## 1. What was open

#598 established, on the widths 4–8:

1. the reported `r_positive = 10, 26, 76, 232, 750` is the orbit count of the
   cut reflection `R : i -> (w-1) - i` (and *every* reflection has the same
   fixed-point count `C(w, floor(w/2))`, so the number alone is not
   identifying);
2. the coarsest exact strong lumping admissible for the D0 dictionary and
   stable for `J` and `D` separately equals the `R`-orbit partition
   **block for block** (Gate-1 branch A);
3. a parity selection rule: for `R`-even sources and readouts the first-order
   response to the odd part `H_odd` of the local out-of-pencil tilt `H_single`
   vanishes pointwise, while `halves_linked` (odd width) is the one declared
   readout that is *not* `R`-even and therefore fires at first order.

Left open in #598 was **Phase C**: factor the finite-horizon balanced
input/output description through the symmetry quotient and compare

* full microscopic state space,
* exact reduction to the reflection-quotient / even sector.

The stated reason it was left open is exactly the right one: the balanced
construction needs an inner-product convention, and a symmetry-incompatible
convention silently changes the answer.  This note pins the convention down,
proves why the factorization is exact for any compatible convention, verifies
it numerically to machine precision on widths 4–8, and shows quantitatively
what goes wrong when the quotient is given the wrong (unweighted) inner
product.

---

## 2. Object and conventions

State space: canonical circular noncrossing connectivity states on `w`
boundary points (`Catalan(w)` states).  Generator pencil `G_eta = G0 + eta H`
with `G0 = J + D` (uniform join/detach rates 1, minus exit diagonal).  Four
declared sources (delta on all-singletons, all-in-one-block, wrapped pair,
uniform).  Eight declared readouts (blocks, singletons, wrap, max_block,
linked_pairs, halves_linked, covering_depth, boundary_span).  Reflection
`R : i -> (w-1)-i` acting on states.

Task objects.  A task is a pair (sources, readouts).  The response kernel is

```
K_{a,j}(t) = <mu_a, e^{tG0} f_j>
```

with `mu_a` a source (probability measure) and `f_j` a readout function.
Readouts are centred by subtracting their counting mean (this makes every
observability direction orthogonal to the constant mode, which is exact
because `e^{tG0} 1 = 1`).

Balanced data.  With `V` the `n x m` source matrix (columns `mu_a`) and `U` the
`n x r` readout matrix (columns the centred `f_j`), the finite-horizon task
Gramians on `[0,T]` are

```
P_c(T) = ∫_0^T e^{tau G0} V V^T e^{tau G0^T} dtau,
P_o(T) = ∫_0^T e^{tau G0^T} U U^T e^{tau G0} dtau,
```

and the **task singular values** are the singular values of `Z^T C`, where `C`
holds columns `sqrt(w_l) e^{tau_l G0} mu_a` and `Z` holds columns
`sqrt(w_l) e^{tau_l G0^T} f_j` over a Gauss–Legendre grid `{(tau_l,w_l)}` on
`[0,T]` (we used `T = 4`, 32 nodes).  This is the standard finite-horizon
balanced object for the impulse-response system `x' = G0 x, x(0+) = V u,
y = U^T x`, with the constant mode removed exactly on the output side.

Inner product.  The full-space inner product is the counting measure
`<a,b> = sum_i a_i b_i`.  It is `R`-invariant: the reflection is a permutation,
and `R` is an isometry of this inner product.

---

## 3. Structural theorem

**Setup.**  `R^2 = 1`; write `V = V^+ ⊕ V^-` for the `±1` eigenspaces of `R`
(the even/odd functions).  `dim V^+` = number of `R`-orbits, `dim V^-` =
(n - fix(R))/2.  The baseline generator satisfies `R G0 R = G0` (verified
exactly, move by move), hence `e^{tG0}` preserves both sectors.

**Lemma 1 (orbit functions = even functions).**  A function is `R`-even iff it
is constant on every `R`-orbit.  Proof: on a two-point orbit `{i, R i}` the
evenness constraint `f_i = f_{R i}` is exactly constancy; fixed points impose
nothing.  Hence `V^+` is canonically isomorphic to the function space of the
orbit quotient, and the counting inner product pulls back to the
**orbit-weighted** product

```
<f,g>_{orbit} = sum_o |o| f_o g_o
```

(weight `|o|` = orbit size), not the plain Euclidean product.

**Lemma 2 (exact quotient dynamics).**  The orbit partition is a strong lumping
for `J` and `D` separately (Gate-1, block-for-block verified on widths 4–8), so
`G0` induces an exact generator `G0^+` on orbit functions, and for orbit-constant
`f`,

```
(G0 f)|_orbit = (G0^+ f_orbit).
```

**Theorem (factorization).**  Let the task have `R`-even sources and `R`-even
readouts only (every declared source is even at all widths; at even widths
every declared readout is even; at odd widths `halves_linked` is the unique
odd readout).  Let `<>` be any inner product for which the reflection is an
isometry and for which the sector decomposition is orthogonal (the counting
and the stationary `L2(pi)` products both qualify).  Then, for every finite
horizon `T`:

1. **The odd sector is exactly uncontrollable and unobservable for the
   task**: the response kernel satisfies `K_{a,j}(t) = 0` whenever the initial
   measure lies in `V^-` and the readout is even, and likewise for even
   measures against odd readouts.  Equivalently the task Gramians have zero
   blocks on `V^-`.
2. **The full-space and the even-sector task singular values coincide**:
   because the full data `(Z, C)` only ever sees `V^+`, the nonzero spectrum is
   carried entirely by the even sector.
3. **The even sector *is* the orbit quotient, provided the orbit-weighted
   inner product is used.**  In `sqrt(|o|)`-coordinates the quotient balanced
   computation reproduces the full-space task singular values exactly.
4. **With the plain (unweighted) Euclidean product on the quotient, the
   factorization fails numerically**: this is not a dynamics effect, it is
   purely the inner-product mismatch identified in #598's decision table.

Any energy-threshold "balanced order" derived from these spectra (e.g. order
for a 90%/99% energy cut) is therefore identical on the full space and on the
quotient, and the small "order 3–4" of the P398 task description is an order of
the *quotient dynamics after task compression*, not of the microscopic space
and not of the exact positive quotient (whose size 10..750 is an orbit count).

**Proof of (2).**  `V^±` are invariant under `e^{tG0}` and orthogonal under
`<>`.  The columns of `C` lie in the span of `{e^{tau G0} mu_a}` ⊂ `V^+`
(even sources).  The columns of `Z` lie in `V^+` (even readouts, and `e^{tG0^T}`
preserves `V^+` because `R` commutes with `G0^T` too).  Hence `Z^T C` depends
only on the even-sector coordinates of every factor, and the odd block of the
Gramians is zero identically.  Claim (1) is the same statement at the level of
single channels.  Claims (3)–(4) are Lemma 1 applied to the coordinate system.

---

## 4. Numerical verification (widths 4–8)

Implementation: independent Python reimplementation of the P398 object
(`scripts/probe/p398core.py`) and of the balanced computation
(`scripts/probe/phase_c_quotient_factor.py`), exact/float as declared.
All numbers in `results/probe-p398-quotient-balanced/latest.json`.

| w | n | n_orbit | even task: full vs quotient (weighted) | quotient with PLAIN product |
|---|------|---------|------------------------------------------|-----------------------------|
| 4 | 14 | 10  | rel. diff **5.2e-16** | rel. diff **0.61** |
| 5 | 42 | 26  | rel. diff **2.6e-16** | rel. diff **0.47** |
| 6 | 132| 76  | rel. diff **2.4e-16** | rel. diff **0.37** |
| 7 | 429| 232 | rel. diff **2.3e-16** | rel. diff **0.40** |
| 8 | 1430|750 | rel. diff **2.3e-16** | rel. diff **0.44** |

Leading task singular values (full space = quotient):

| w | spectrum head |
|---|---------------|
| 4 | 0.4196, 0.1104, 0.0399, 0.0082, ... |
| 5 | 0.4496, 0.0936, 0.0433, 0.0195, ... |
| 6 | 0.5067, 0.0904, 0.0453, 0.0331, ... |
| 7 | 0.5605, 0.0842, 0.0507, 0.0415, ... |
| 8 | 0.6406, 0.0861, 0.0658, 0.0438, ... |

The relative differences are at the level of the quadrature/uniformization
round-off (1e-16).  The **plain-product column is the diagnosis**: giving the
quotient its orbit weights is not cosmetic, it changes the task singular data
by tens of percent.  This is the quantitative form of #598's warning that the
balanced construction must fix a symmetry-compatible convention before the
"order 3–4" is quoted.

**Anti-invariant sector.**  Peak response of the odd readout `halves_linked`
(odd widths, where it is the unique odd readout) from the four even sources:
5.0e-18 (w=5) and 2.4e-18 (w=7) — i.e. zero to machine precision.  From a
signed odd source (a delta difference on one two-point orbit, which is *not* a
probability measure and therefore a deliberately marked channel) the same
readout responds at 0.029 (w=5) and 0.0063 (w=7).  This is the exact finite
realization of #598's decision-table row "odd sector becomes visible only
after adding a marked readout" — with the precision that *a single marked
readout is not enough*: an odd readout is invisible to even sources at every
time.  Visibility requires a marked channel on at least one side, or a
parity-odd perturbation as in #598 Phase A.

---

## 5. What this does to the wider picture

The P398 complexity stack now reads cleanly as three successive, exactly
separated factors:

```
microscopic Catalan state n = 14..1430
   --[exact symmetry quotient / strong lumping]-->  orbit space (10..750)
   --[finite-horizon task compression, R-compatible]-->  balanced data (order 3-4)
```

and each step is either exact (lumping = orbit partition, Gate-1) or verified
to machine precision (balanced factorization).  In particular:

* `r_positive` (10..750) is an *orbit count*, i.e. a symmetry-quotient size,
  and is not a second, independent "compression"; the Gate-1 partition
  equality carries that claim, not the count formula.
* The balanced order is a *task-relative* number attached to the quotient
  dynamics: identical on the full space and the quotient (this note), which is
  the precise sense in which "750" and "3–4" are not competing estimates of
  one state dimension (they live on different rungs of the factorization).
* The inner product must be chosen on the full space first; the quotient then
  inherits orbit weights automatically.  #598 Phase C is thereby closed under
  the counting-measure convention, with the plain-product mismatch quantified
  as the diagnosed failure mode.

**Claim boundary.**  All statements are exact finite statements about the P398
process and the declared dictionary.  Reflection parity is a `C2` label, not a
continuum spin.  The balanced order remains readout/source-relative.  Nothing
here is a site-percolation statement.

Files: `scripts/probe/p398core.py`, `scripts/probe/phase_c_quotient_factor.py`,
`scripts/probe/verify_gate1_independent.py`,
`results/probe-p398-quotient-balanced/latest.json`.
