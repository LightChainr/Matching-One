# The state was never fragile: P398's positive quotient is a reflection orbit, and the breaking direction is first-order invisible

**Date:** 2026-09-06
**Issues:** #598 (Phases A and B); closes #596's Gate 1 on branch **A**; revises the headline of #580 and #594
**Evidence type:** exact rational arithmetic and exact partition refinement, plus a uniformization response ladder. No sampling, no new width, no fit.

## What #580 reported and what it actually was

#580 found that a localized out-of-pencil tilt (`single_point_join`) destroys P398's coarsest exact positive lumping at every width — the admissible partition collapses to the identity — while the frozen finite-horizon input/output description barely moves. That was written up as *the positive state is fragile*.

It is not fragile. It is a **symmetry quotient**, and symmetry quotients are destroyed by any symmetry-breaking perturbation as a matter of course. What needed explaining was never the collapse; it was why the task did not notice.

Both halves now have exact answers.

## Half one — the quotient is the orbit partition of one identified reflection

The owner's observation in #598 is that the reported block counts match Burnside for a single reflection:

```text
w   Catalan   C(w,floor(w/2))   [C+b]/2    r_positive reported by #580
4      14            6              10          10
5      42           10              26          26
6     132           20              76          76
7     429           35             232         232
8    1430           70             750         750
```

Five for five. **But this identifies nothing on its own**, and that has to be said before the result is used: *every* reflection of the cycle fixes exactly `C(w, floor(w/2))` noncrossing states, so every reflection gives this same orbit count. The formula is a necessary condition, not a fingerprint.

What names the reflection is the **task**. Of the eight declared readouts, seven are fully dihedral-invariant; `wrap` = `[state[0] == state[w-1]]` is not, and it is preserved exactly by the reflection that fixes the cut between points `w-1` and `0`:

```text
R = ( i -> (w-1) - i ),   the unique reflection preserving the declared D0 dictionary
                          at every width 4..8.
```

The decisive check is then the one #596's Gate 1 actually asks for, and it is stronger than cardinality — two partitions can both have 750 blocks and share none. Computing the coarsest strong lumping admissible for `D0` and stable for `J` and `D` separately, by exact partition refinement:

```text
w   initial blocks   coarsest   R-orbits   coarsest == orbit partition
4         8             10         10              YES
5        12             26         26              YES
6        18             76         76              YES
7        26            232        232              YES
8        37            750        750              YES
```

Not the same size — **the same partition**. Gate 1 resolves on branch **A (symmetry orbit)**.

A consistency check that could have failed and did not: at odd widths `halves_linked` is *not* `R`-invariant, so demanding the lumping carry all eight readouts must destroy the quotient. It does — at `w = 5` the coarsest admissible partition is the identity, all 42 blocks. The quotient exists exactly when the dictionary is `R`-invariant, which is what "task-compatible" has to mean if it means anything.

## Half two — the parity selection rule, checked pointwise

Write the localized tilt in its parity parts. `R` maps the join at point `p` to the join at `(w-2) - p`, so the mirror of `join@0` is `join@(w-2)` — verified exactly at the level of the move maps, not just the assembled operators:

```text
H_even = (H_single + R H_single R)/2 = ( join@0 + join@(w-2) ) / 2,
H_odd  = (H_single - R H_single R)/2 = ( join@0 - join@(w-2) ) / 2.
```

Exact rational arithmetic confirms `R G_0 R = G_0`, `R J R = J`, `R D R = D`, `R H_even R = +H_even`, `R H_odd R = -H_odd`, and that the two halves sum back to `H_single`. **The odd part is not small**: `||H_odd|| / ||H_single|| = 0.553, 0.559, 0.560, 0.560, 0.560` for widths 4–8. Over half the breaking direction is odd.

The claim is that an invariant task cannot feel it at first order. The first derivative of a response is

```text
d/de [ mu^T e^{t(G + e H_odd)} f ]_(e=0)  =  int_0^t  mu^T e^{(t-s)G} H_odd e^{sG} f  ds,
```

and the theorem is that **every integrand vanishes**, not merely the integral:

```text
mu even  ->  e^{(t-s)G^T} mu  stays even     (G is R-equivariant)
f  even  ->  e^{sG} f          stays even
H_odd    :   even -> odd
<even measure, odd function> = 0.
```

Checked pointwise over a declared grid (`t` in {0.5, 1, 2}, `s/t` in {1/4, 1/2, 3/4}), for all four declared sources — all of which are `R`-even — against every readout:

```text
   w  states       R      protected pairs   worst relative integrand   exposed pairs
   4      14   i -> 3-i         288               6.57e-15                  0
   5      42   i -> 4-i         252               2.14e-16                 36
   6     132   i -> 5-i         288               1.12e-15                  0
   7     429   i -> 6-i         252               1.05e-15                  0 (36 exposed)
   8    1430   i -> 7-i         288               5.53e-16                  0
```

**The test can fail, and it does, exactly where it must.** At `w = 5` and `w = 7`, `halves_linked` is odd-containing, so the rule *requires* a nonzero first-order response there. Observed relative integrands: **0.292** at `w = 5`, **0.107** at `w = 7`. A run in which everything vanished would have been a broken test, not a stronger theorem. This is also the positive control #598's decision table asks for, and it came free from the existing dictionary rather than from a readout invented for the purpose.

## Phase B — the consequence, at second order

If the rule holds, turning on the pure breaking direction must move an invariant task at `O(e^2)`, while the even part moves it at `O(e)`. Over the declared ladder `e = 1/16, 1/8, 1/4`, log-log slopes on the even readouts:

```text
   w   odd vs none   single vs even   even vs none
   4      2.0006         1.9512          0.9800
   5      2.0007         1.9485          0.9812
   6      2.0008         1.9460          0.9820
   7      2.0009         1.9438          0.9825
   8      2.0010         1.9421          0.9828
```

The `even vs none` column is the control that the machinery can see a first-order effect at all — without it, "everything is small" would be indistinguishable from a bug.

The magnitudes are the quantitative version of #580's puzzle. At `e = 1/4`, width 8: the odd direction moves the invariant task by `9.9e-4` relative, the even direction by `2.55e-2`. **A factor of 26.** That is why a frozen realization survived an intervention that annihilated the exact quotient — it never had to represent the part that did the annihilating.

## What this changes

**#580's headline is wrong and should be replaced.** Not "exact positive realization is fragile while signed descriptions are robust". Rather:

> P398's exact positive quotient *is* the orbit partition of the one reflection the declared dictionary preserves. Any symmetry-breaking intervention destroys it discontinuously. But over half of that intervention lies in the odd sector, and an invariant finite-horizon task is exactly blind to the odd sector at first order — so the exact quotient can die while the task-relative description does not move.

That is a transportable finite theorem, and it is not special to P398. The general statement is the `C2` case of a character selection rule: with a `K`-equivariant baseline, sources and readouts in irreps `rho_B`, `rho_C`, and a perturbation in `rho_H`, the first-order response is allowed only when the trivial representation occurs in the relevant tensor product. That is the same shape as the project's deck-character rule in #244 — an analogy at the level of representation theory, *not* an identification of the two physical systems.

**#594's Q1 is answered, and the answer is not the conjecture in its body.** The asymmetry between "three bounded descriptions survive, one closed condition dies" does not need a codimension argument about lumpability-preserving varieties. The closed condition is a symmetry, the surviving descriptions are invariant, and parity does the rest. The owner's `ker L_P` tangent-space computation would still be worth doing to characterize *which* interventions a lumped model may transport in general, but it is no longer needed to explain P398.

**#593 changes character.** It was "extend the memory-degree table to widths 9 and 10". It is now a falsification test of a closed form, with the prediction stated before the computation:

```text
r_positive(9)  = (4862 + 126)/2  = 2494
r_positive(10) = (16796 + 252)/2 = 8524.
```

If DeepSeek's exact refinement returns those, the orbit account holds two widths beyond where it was fitted. If not, the account is wrong and the state semantics need repair before anything else proceeds. Either way that is worth far more than two more rows.

## Phase C — closed by #600 / PR #603

Phase C (does the balanced realization factor through the quotient?) was not run here: it
needs finite-horizon Gramians and Hankel SVDs on a 1430-state operator twice, which is the
wrong job for a pure-Python environment. It was handed to #600 and came back affirmative.

Protected `D0`, microscopic vs quotient construction, `w = 4..8`:

```text
w   states -> orbits   max |A-B| Hankel spectrum   response rel-Frobenius
4      14 ->  10             1.5e-16                    9.7e-16
5      42 ->  26             3.3e-16                    1.3e-15
6     132 ->  76             3.0e-16                    2.6e-15
7     429 -> 232             3.4e-16                    4.2e-15
8    1430 -> 750             2.7e-15                    1.1e-14
```

Balanced order and numerical rank preserved at every width. The odd sector is **exactly
inert** — reach/observe energy `1.7e-17` to `5.7e-17` across `dim = 4/16/56/197/680` — which
is why the reduction is exact rather than accurate.

The declared trap fired. `halves_linked` is `R`-even only at even widths, and the agreement
breaks at exactly `w = 5` (`1.8e-2`) and `w = 7` (`8.7e-3`), recovering to `1e-16` under
symmetrization. The break is the `R`-odd part and nothing else.

So the factorization theorem is established in the form #598 proposed:

```text
exact task symmetry quotient  ->  task-relative balanced reduction,
```

the first an exact factor, the second a lossy compression of the quotient dynamics. That
is a much cleaner statement than comparing `r_positive` and balanced order as two
unrelated dimensions.

## Claim boundary

- This is an exact finite symmetry statement about P398. It is **not** a percolation-threshold result.
- Reflection parity here is a `C2` label on a finite connectivity process. It is not a continuum spin label.
- The small balanced order remains task-relative; nothing here promotes it to a system state dimension.
- The connection to #244 is a representation-theoretic analogy, not a claim that the two systems are the same.
- The orbit-count identity by itself identifies nothing. The load is carried by the partition equality and the selection rule.

## Reproduction

```bash
python3 scripts/p398_reflection_parity.py
```

writes `results/p398-reflection-parity/latest.json` (~3 min, widths 4–8). Tests:
`python3 -m unittest tests.test_p398_reflection_parity` (17 tests).
