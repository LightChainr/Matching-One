# Second-wave research: beyond the known matching estimator

Status: speculative research program, with literature corrections and falsifiable experiments.

This note deliberately separates three things:

1. results that are already known in the literature;
2. deductions that follow from those results with modest assumptions;
3. aggressive conjectures that are worth trying because they can be killed cheaply.

The target is not another decimal obtained by a fragile fit. The target is a structural mechanism that either removes finite-size fields predictively, constrains a putative exact solution, or yields a new certified statement.

## 0. Literature correction: do not rediscover the 2016 result

Mertens and Ziff (2016) already derived the exact finite-size matching identity on an `L x L` torus. For the square site lattice and its NN+NNN matching lattice, define

\[
M_L(p)=N_L(p)-\widehat N_L(1-p)-L^2\chi(p),
\qquad
\chi(p)=p-2p^2+p^4.
\]

They proved that this equals a difference of primal and matching wrapping probabilities (for several equivalent wrapping conventions). Thus the finite-torus identity itself is not a new target; our task is to reproduce it as a regression/property test and then go beyond it.

They also found empirically:

- the root `M_L(p)=0` approaches `p_c` approximately as `L^-4`;
- assuming the leading critical behavior `M_L(p_c) ~ L^-13/4`, the condition

  \[
  L^{13/4}M_L(p)=(L-1)^{13/4}M_{L-1}(p)
  \]

  produced an error close to `L^-7` on their available small systems.

So a generic proposal to "use matching to cancel the leading correction" is already behind the literature. The new question is:

> **What finite-size operator produces the residual after the known `L^-7` acceleration, and can it be removed without fitting the answer?**

References:

- S. Mertens and R. M. Ziff, *Percolation in finite matching lattices*, Phys. Rev. E 94, 062152 (2016), arXiv:1603.07289, DOI 10.1103/PhysRevE.94.062152.
- C. R. Scullard and J. L. Jacobsen, *Bond percolation thresholds on Archimedean lattices from critical polynomial roots*, Phys. Rev. Research 2, 012050(R) (2020), arXiv:1910.12376. They find two observed critical-polynomial correction classes, `(6,7,8)` and `(4,6,8)`, which is a warning that the correction spectrum is observable/geometry dependent.

## 1. P0 hypothesis A: orientation is an operator projector

### Evidence

Feng, Deng, and Blöte (2008) found the second thermal scaling dimension

\[
X_{t2}=4
\]

and evidence for a correction of the form

\[
(b+d\log L)L^{-2}
\]

in scaled-gap observables. More importantly for this project, the non-logarithmic power-law amplitude depends strongly on the orientation of the square lattice relative to the cylinder.

For exactly solved square **bond** percolation they report, in their lattice length conventions,

| orientation | power amplitude `C` | log amplitude `A` |
|---|---:|---:|
| parallel | `+0.0306(1)` | `-0.0054(1)` |
| diagonal | `-0.0205(1)` | `-0.0027(1)` |

They explicitly note that the factor two between the log amplitudes is consistent with the `sqrt(2)` difference in length units, while the power-law amplitude changes with orientation.

This is exactly the pattern one would expect if the leading correction contains an anisotropic square-lattice component plus an isotropic/logarithmic component. In CFT language, a natural suspect is a spin-4 irrelevant field: under a rotation by `pi/4`, a pure spin-4 harmonic changes sign.

This spin assignment is a **hypothesis**, not a claim proved by the 2008 data.

Reference:

- X. Feng, Y. Deng, H. W. J. Blöte, *Percolation transitions in two dimensions*, Phys. Rev. E 78, 031136 (2008), arXiv:0901.1370, DOI 10.1103/PhysRevE.78.031136.

### Aggressive conjecture A1: square-site leading bias is mostly anisotropic

For a properly normalized dimensionless estimator `E_theta(L)` on a square geometry rotated by angle `theta`, write schematically

\[
E_\theta(L)=p_c + L^{-\Delta}
\left[a_0+a_4\cos(4\theta)\right]
+L^{-\Delta}\log L\,a_{\log}+\cdots.
\]

For `theta=0` and `theta=pi/4`, the spin-4 term changes sign. Therefore

\[
E_+(L)=\tfrac12(E_0+E_{\pi/4})
\]

projects away the pure spin-4 contribution after physical-length normalization, while

\[
E_-(L)=\tfrac12(E_0-E_{\pi/4})
\]

isolates it.

**Strong version:** after combining orientation projection with the matching observable, the leading surviving correction is smaller by at least two powers of `L` than in either geometry alone.

**Why plausible:** the exactly solved square-bond benchmark already shows a large orientation-sensitive amplitude. After the length-unit conversion suggested by Feng et al., its two power amplitudes are not equal and opposite, so a scalar admixture remains, but the anisotropic part appears large enough to be worth projecting.

**Falsification:** on square bond percolation at exact `p_c=1/2`, the orientation-projected estimator does not improve held-out finite-size convergence. If it fails on the exact benchmark, do not apply it as an accuracy trick to square site.

### Experiment A

Use three controls before touching the unknown threshold:

1. square bond, `p_c=1/2`: positive control for square anisotropy;
2. triangular site, `p_c=1/2`: high-rotational-symmetry control; Feng et al. found the logarithmic amplitude small/absent there;
3. square site: target.

For each, compute the same observable on axis-aligned and diagonal/diamond periodic geometries. Normalize `L` by the actual physical circumference before comparing amplitudes.

Do not merely compare fitted `p_c`. Fit the **amplitude transformation law**. The valuable output is a quantitative statement such as

\[
a_4(\theta+\pi/4)=-a_4(\theta)
\]

within an error budget, or its rejection.

### Extension A2: angular tomography

Do not stop at two orientations. Use integer-period/sheared tori to sample several modular shapes/orientations and fit angular harmonics

\[
C(\theta)=C_0+C_4\cos4\theta+C_8\cos8\theta+\cdots.
\]

A stable harmonic decomposition would identify the lattice-spin content of the correction directly. This is more informative than fitting a list of unrelated inverse powers.

## 2. P0 hypothesis B: annihilate the residual after the known `L^-7` estimator

Let

\[
S_L(p)=L^{13/4}M_L(p).
\]

The Mertens-Ziff two-size condition is simply

\[
S_L(p)-S_{L-1}(p)=0,
\]

a first discrete annihilator of an `L`-independent critical amplitude.

This viewpoint suggests a systematic generalization.

Assume, only as an asymptotic model to be tested,

\[
S_L(p_c)=A_0+A_1L^{-q_1}+A_2L^{-q_2}+\cdots
\]

and near criticality

\[
\partial_p S_L(p_c)\asymp L^4,
\]

because `M'_L(p_c) ~ L^(1/nu)=L^(3/4)` and `1/nu=3/4`.

Choose weights `w_i` on consecutive sizes so that

\[
\sum_i w_i=0,
\qquad
\sum_i w_i L_i^{-q_j}=0
\]

for one or more preregistered correction exponents. Then define the pseudo-critical point by

\[
\sum_i w_i S_{L_i}(p)=0.
\]

This is a finite-difference/Richardson filter applied to the **matching function**, not a free polynomial extrapolation of the final threshold.

### Aggressive conjecture B1: there is a stable second annihilation

The observed `~L^-7` behavior is evidence that after the leading `L^-13/4` term in `M_L(p_c)`, the scaled function `S_L(p_c)` has a reasonably simple next correction. We should estimate `q_1` from exact/very-high-statistics data, freeze it to a small rational candidate only if the evidence supports that, and then test a 3-size annihilator on held-out sizes.

Do **not** assume `q_1=3` merely because a naive adjacent-size argument maps a `L^-7` root to that value. Small-`L` effective exponents are notoriously deceptive; Mertens and Ziff explicitly warn that adjacent-size slopes converge slowly.

Candidate families to test in preregistered order:

- integer `q`: `2, 3, 4, 5, ...`;
- logarithmic companions `L^-q log L` motivated by the `c=0` correction structure;
- geometry-dependent classes suggested by the `(4,6,8)` versus `(6,7,8)` critical-polynomial spectra.

### Success criterion B

A second annihilator counts as real only if all are true:

1. exponent/filter chosen without the largest sizes;
2. weights are frozen;
3. the filtered root improves prediction on at least three withheld size steps;
4. improvement appears in both exact/symbolic small systems and an independent Monte Carlo or transfer-matrix family;
5. numerical conditioning is reported. Huge alternating weights that amplify noise are not an improvement.

A helper script `scripts/correction_filter.py` in this branch generates such weights and reports their noise-amplification factor.

## 3. P0/P1 synthesis: matching x orientation may remove different sectors

The strongest conjecture in this note is that **matching symmetry and lattice rotation project different irrelevant sectors**.

If so, a four-channel measurement

- primary square, axis geometry;
- matching NN+NNN, complementary occupation, axis geometry;
- primary square, diagonal geometry;
- matching NN+NNN, complementary occupation, diagonal geometry;

could separate at least:

- matching-even versus matching-odd amplitudes;
- rotation-even versus rotation-odd (`spin 4 mod 8`) amplitudes;
- a residual scalar/log sector.

This is a small representation-theory experiment masquerading as finite-size scaling. The objective is to build projectors, not to add fit parameters.

If the amplitude table approximately factorizes into matching parity and rotation parity, that is a structural discovery even before it improves `p_c`.

## 4. P1 hypothesis C: search for a "magic" torus shape where the leading amplitude is zero

Finite-size amplitudes are not universal constants; they depend on shape and boundary conditions. Treat the torus modular parameter (aspect ratio plus shear) as a control variable.

### Conjecture C1

There exists a rectangular or sheared periodic geometry `tau_*` for which the leading correction amplitude of the matching estimator crosses zero:

\[
a_1(\tau_*)=0.
\]

If found, `tau_*` is an **improved geometry**: no cross-size fitted cancellation is needed.

### Search

1. At modest `L`, scan aspect ratios and integer shears.
2. Estimate the signed leading amplitude using exact square-bond `p_c=1/2` first.
3. Locate sign changes in geometry space.
4. Freeze a rational/integer geometry near the zero.
5. Test at larger sizes.

This is GPU-friendly because many shapes/sizes can be evaluated with the same many-replica kernel, but discovery can start on CPU.

## 5. P1 hypothesis D: geometry diversity is more valuable than `n=25` alone

The present threshold dispute demonstrates that one long finite-width sequence can support an overconfident extrapolation. A better information strategy may be to compute **several independent geometries at moderate widths** and impose a shared `p_c` with geometry-dependent amplitudes.

Before spending large memory to extend a single transfer-matrix sequence from width 24 to 25+, compare the expected information gain from:

- axis cylinder;
- diagonal cylinder;
- helical cylinder;
- square torus;
- sheared torus;
- matching-function roots;
- exact finite-square median/cell estimators.

If two geometries couple differently to the leading fields, width `18` in a new geometry can be more diagnostic than width `25` in an old one.

## 6. P1 hypothesis E: exact-solvability obstruction by polynomial factors

Mertens and Ziff make an important algebraic observation: for exactly solvable cases such as triangular bond percolation or the martini lattice, the finite matching polynomial/function has a common low-degree factor whose physical root is the exact threshold, for every finite size in the applicable construction.

That suggests an exact-computation program:

1. construct exact integer matching/critical polynomials for small square-site bases;
2. factor them over `Z[p]` (or factor modulo several primes and reconstruct);
3. compute pairwise and all-size gcds;
4. search for stable nontrivial factors and stable low-degree resultants across embeddings/orientations.

### What a negative result means

If the stabilized gcd is `1`, this does **not** prove `p_c` transcendental or non-algebraic. It does rule out a much narrower and more meaningful class: an exact threshold produced by a basis-independent finite local factor of this matching/critical-polynomial family.

This is a better negative result than blind PSLQ because it tests a mechanism known to explain actually solvable lattices.

### Strong positive signal

A factor or subresultant pattern that recurs across sizes and predicts another observable would be extremely important, even if its exact physical root is not yet established.

## 7. P2 hypothesis F: a no-go lemma for ordinary independent-bond site gadgets

A Bernoulli site interacting with four external neighbors has an extreme terminal-connectivity distribution:

- with probability `1-p`, it connects none of the terminals through the site;
- with probability `p`, it connects all incident terminals through the site.

There are no partial terminal partitions generated by the site itself.

A finite gadget built from ordinary independent bonds with every edge probability strictly between 0 and 1 generically assigns positive probability to partial terminal connections. This suggests a simple no-go theorem:

> No nondegenerate finite independent ordinary-bond gadget with four distinct terminals can exactly reproduce the two-state terminal partition law of a single Bernoulli site for all `p`.

This needs a careful graph-theoretic proof and careful statement of exceptions (zero/one-probability edges, identified terminals, hyperedges, correlated bonds).

If true, it explains why searching naive finite bond replacements for an exact square-site solution is structurally misguided. Exact mappings must introduce correlation, hyperedges/multispin terms, or a more global transformation.

After proving the no-go statement, gadget search should switch from "exact representation" to:

- stochastic domination and rigorous bounds;
- approximate moment/connectivity matching;
- correlated/hyperedge self-dual constructions.

## 8. P2 hypothesis G: deform square site into an exactly solvable extended model

A more radical route is to enlarge the model rather than force the one-parameter site problem into a bond formula.

Historical duality/decimation transformations map site percolation to Potts-like models with multispin/checkerboard interactions. The square-site point may therefore be viewed as a point in a higher-dimensional coupling space where self-dual/integrable manifolds can exist even though the original one-parameter slice is not self-dual.

Research question:

> Is there a low-dimensional correlated plaquette/hyperedge family containing square site percolation for which a duality map is explicit and the critical manifold can be traced from an exactly solved point to the square-site slice?

The useful output need not be a closed form. A convergent analytic continuation, rigorous bracketing path, or identified obstruction would already be progress.

This is high-risk and should not consume large compute until the finite-size/operator program is mature.

## 9. PSLQ remains last

The 2024 literature contains competing ultra-high-precision extrapolations near the last few quoted digits. Therefore ordinary PSLQ against a rounded decimal is especially dangerous.

Only resume closed-form searches after:

- independent estimator families agree;
- the interval accounts for extrapolation model uncertainty;
- any basis of constants is structurally motivated;
- interval arithmetic and false-positive controls are used.

The strongest algebraic search before then is the finite-polynomial factor/GCD program in section 6.

## 10. Distinguish correction exponents that describe different observables

Do not mix the following merely because all are called "corrections to scaling":

- `X_t2=4` and associated `L^-2` / `L^-2 log L` corrections in finite-size CFT observables;
- `Delta=4,6,8` or `6,7,8` empirical convergence powers of particular critical-polynomial threshold sequences;
- the exact cluster-size correction exponent `Omega=72/91` derived for 2D percolation cluster-size distributions (Xu et al., Phys. Rev. E 111, 034108, 2025).

They live in related scaling theory but apply to different quantities. Any claimed identification must be derived, not made by matching numbers.

Reference:

- Y. Xu, T. Chen, Z. Zhou, J. Salas, Y. Deng, *Correction-to-scaling exponent for percolation and the Fortuin-Kasteleyn Potts model in two dimensions*, Phys. Rev. E 111, 034108 (2025), arXiv:2411.12646.

## 11. Compute priority for the available hardware

### On the 8-core / 16 GB Debian host: start here

1. reproduce the Mertens-Ziff exact identity for exhaustively enumerable tiny tori;
2. implement/validate `M_L` with shared primal/matching configurations;
3. generate moderate-size CPU Monte Carlo for axis and diagonal geometries;
4. calibrate orientation projection on square bond `p_c=1/2`;
5. fit amplitude parity, not just `p_c`;
6. test correction annihilators with frozen weights;
7. exact-factor/gcd experiments for the smallest computable polynomials.

### RTX 5090-class GPU: rent only after the CPU reference passes

Best use:

- thousands of independent replicas in parallel;
- bit-packed occupancy and wrapping tests;
- narrow `p` scans with common random numbers;
- many `L`, orientation, aspect-ratio, and shear points;
- covariance-rich four-channel matching/orientation measurements.

Poor first use:

- irregular transfer-matrix state hashing with rapidly growing connectivity state spaces.

For transfer matrices near or beyond width 24, memory capacity and state representation are likely more important than raw GPU FLOP/s. Profile bytes/state and transition fan-out before renting hardware.

## 12. Decision gates

### Gate A — orientation

Advance if the exact square-bond control shows a stable orientation-odd amplitude and a frozen projection improves held-out convergence. Otherwise demote the spin-4 projection idea.

### Gate B — second annihilator

Advance if a preregistered multi-size filter beats the known two-size `~L^-7` construction out of sample without pathological weight amplification.

### Gate C — combined projector

Advance to large GPU runs only if matching parity and rotation parity produce a reproducible amplitude table on CPU.

### Gate D — exactness mechanism

Advance algebraic/gadget routes only when they produce a stable factor, a proved no-go statement, a certified bound, or a prediction beyond `p_c` itself.

## 13. Current ranking

| rank | path | probability of useful result | probability of exact `p_c` mechanism | compute fit |
|---|---|---:|---:|---|
| 1 | orientation/spin-4 tomography + matching | high | low-medium | excellent |
| 2 | multi-size annihilator beyond known `L^-7` | high | low | excellent |
| 3 | magic torus shape / zero amplitude | medium-high | low | excellent, GPU later |
| 4 | multi-geometry joint scaling | high | low | good |
| 5 | exact polynomial factor/GCD obstruction | medium | medium if positive | good |
| 6 | correlated/self-dual gadget approximants | medium | medium | good |
| 7 | rigorous substitution/gadget bounds | medium for theorem | very low for exact value | CPU/high-RAM |
| 8 | extended multispin/integrable deformation | low | potentially high | theory first |
| 9 | frontier transfer matrix width only | high for digits | low | high-RAM bottleneck |
| 10 | elementary PSLQ | very low | very low | cheap but misleading |

The key strategic shift is this: **stop asking which inverse powers fit the same sequence best, and start constructing transformations of observables/geometry that force selected amplitudes to vanish.**
