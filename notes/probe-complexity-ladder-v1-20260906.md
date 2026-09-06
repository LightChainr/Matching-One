# C9: complexity-ladder theorems v1 — every arrow labelled

**Status:** consolidated v1 over C1–C8 + rounds 1–2.  Definitions and proofs
live in the cited notes; this document is the *ladder*: which complexity
notion bounds which, on which object, under which assumptions.  Arrow labels:
**P** = proved exact, **F** = false (explicit counterexample), **C** = true
under stated conditions, **U** = unknown (with cheapest falsification).

## 0. The notions (all finite, all object-attached)

| symbol | meaning | where defined/measured |
|---|---|---|
| n_micro | microscopic state count | P398: Catalan(w)=14..1430 |
| n_orb | symmetry-orbit count (C2) | P398: 10..750 (= r_pos, Gate-1) |
| r_pos | exact strong-lumping / positive-quotient size | P398: = orbit partition (block equality) |
| r_lin | signed linear realisation dimension | P398 w4–7: 10,26,72,218; w8: open (definition-blocked, see rlinear note) |
| r_bal | finite-horizon balanced order | P398: 3–4 (task-relative, convention-dependent values, quotient-stable) |
| r_mem | memory/MZ order | P398: effective 2–4; numerical 4,9,12,13,14 (w4–8, open at w9) |
| kappa_L | predictive/Nerode class count of a language | #549: k+1; P398 dictionaries: 8..209 static |
| r_L | response rank of a language | #549 single-root: ≤2; depth-d: ≤d+1 |
| d_shape | full-law finite-size tangent dims | #582/#584: 1 dominant + small structured remainder (smooth-null caveat) |

## 1. Main ladder (on the repository objects)

| arrow | label | statement / counterexample |
|---|---|---|
| n_orb <= n_micro | P | orbits refine states |
| r_pos = n_orb | P (C2, D0-admissible, J/D-stable) | Gate-1 partition equality, w4–8; both directions of the caveat recorded (orbit-count formula alone not identifying) |
| r_lin <= r_pos | C | holds w4–7 (10≤10,26≤26,72≤76,218≤232); w8 open; FALSE as a theorem: w6 strict (72<76) already shows they are different functionals |
| r_pos - r_lin grows | U | w4–7 differences 0,0,4,14; w8 needs exact r_lin |
| r_bal <= min(r_pos, r_lin) in spirit | C | balanced order 3–4 while r_pos ≥10: task-relative compression after the symmetry factor; values convention-dependent (C5) |
| full spectrum = quotient spectrum (even task) | P | counting + L2(pi), w4–8, ~1e-16; plain-Euclidean quotient product F (37–61% mismatch) |
| odd sector inert for even tasks | P | responses ≤1e-17; visible only with marked channels or odd perturbation |
| r_mem bounded in w at fixed accuracy | C | energy order 99/99.9% saturates at 3–4; numerical order (tol 1e-6) grows 4,9,12,13,14 => "bounded" is accuracy-dependent (#588) |
| r_mem(w) growth law | U | w9/w10 needed (cost probe: feasible with numpy; gate = convention port, C3) |
| kappa_L > r_L possible | P | #549: kappa=k+1, r≤2 (affine obstruction) |
| kappa_L = r_L attainable | P | depth-k root language: rank k+1 (C2 abstract exact) |
| r_L grows with language depth d | P | degree(P)=d, tests d=0..k => rank k+1; clone width at fixed depth does not help |
| complete-survival language rank small while branching rank grows | P | unbranched = depth-0 sector (rank ≤1 in a); fork depth-1 rank ≤2; #435 identical S(z) |
| any K-invariant language separates K-orbits | F | impossible (Theorem 1/2 round 2) |
| character-marked perturbation visible iff triv in rho_B⊗rho_H⊗rho_C | P | Schur-based, pointwise; C2 (#598) + C3 demo + C5 formal |
| even = orbit-constant for general K | F | holds for C2 only (C2-accident); general K needs isotypic fibres |
| finite even task constrains hidden-sector threshold | F | C7 no-go: task data = 1e-16 constant in g, gap g·m -> 0 |
| 99% one-direction quantile flow => 1-d RG state | F | smooth Taylor null (generic rank-1 local diffs) |
| bounded unbranched task order constrains p_c | F/C | #435/#549 unbranched-identical, fork-different; no bridge known |
| E_p[X]=0 implied by low task order | U | no implication found (homological-balance program separate) |

## 2. Counterexample archive (each with its sharpest number)

| # | counterexample | key datum | file |
|---|---|---|---|
| X1 | r_lin ≠ r_pos on P398 (w6: 72 vs 76) | exact mod-p elimination vs lumping | rlinear note / #593 |
| X2 | orbit-count formula not identifying | every reflection fixes C(w,⌊w/2⌋) states | verify_gate1 / #598 |
| X3 | quotient breaks under plain inner product | 37–61% spectrum change | round-1 |
| X4 | even task blind to odd sector at any time | response ≤5e-18 (w5/7 halves_linked) | round-1/C5 |
| X5 | odd readout alone not enough; marked channel needed | 0.029 (w5) from odd source | round-1 |
| X6 | kappa=k+1 with rank=2 (#549 single-root) | F_{k,a} affine in a; gap 1/[2k(8k-1)²] | C2 |
| X7 | rank growth needs depth not width | degree(P)=d; h_A=h_B degenerates to rank 1 | C2 |
| X8 | task identical, hidden-sector gap movable | task rel diff ≤7e-16; gap g·m -> 0 | C7 |
| X9 | unbranched identical, fork different | #435: 95/196 vs 93/196; S(z) equal | #435/#549 |

## 3. Unknown arrows and their cheapest falsification (ranked)

1. **Exact r_lin at w=8 (and does r_lin > r_pos occur?)** — needs the #593
   definition of r_linear; one mod-p elimination then settles it.
2. **r_mem at w9/w10: saturation vs slow growth** — cost-probe shows the
   numpy path is cheap; the gate is porting the frozen projected-memory
   construction (continuity gate vs w8), decision table in the C3 note.
3. **D2-static signature gap vs orbit count closes dynamically?** — run the
   D2 signature under the P398 dynamics (Krylov-style refinement): cheapest
   exact test of "language sharpens with time" on P398.
4. **Compositional depth vs algebraic degree for general tree languages**
   (#550 witness): formal extension of C2's root-count algebra; needs a
   definition of compositional depth first.
5. **p_c constrained by any current finite task?** — after C7 the burden is
   a positive-sector bridge; cheapest check is whether the repository's
   declared readouts see the homology/topological sector that E_p[X]=0
   targets (this is the #581/#584 typed-labels direction).

## 4. Boundaries

* P-arrows are exact finite statements on the named objects; none is a
  percolation-threshold statement; C2 parity is not continuum spin; balanced
  and memory numbers are convention/task-relative by construction (C5/C8).
* The ladder is a *v1*: the U arrows are the acquisition queue of the probe.

Related: all probe notes C1–C8 (list in `PROGRESS.md`).
