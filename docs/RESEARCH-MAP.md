# Research Map

**Updated:** 2026-08-29

This is a navigation layer. It does not grant permission, lock tasks, require branch consolidation, or declare lower-priority work false. `STATUS` owns claim language; `ROADMAP` ranks attention; `NEXT-TARGETS` is the fast decision board.

## The project in one page

Matching One now has enough evidence that the main uncertainty is no longer “does a finite-size signal exist?” The useful questions are:

```text
Q1. Which continuum sector produces the global matching H4 signal?
Q2. Why does the lattice matching observable select that sector?
Q3. Can any of this explain or constrain the microscopic value of p_c itself?
```

Everything else in the repository should be readable as evidence, a discriminator, a control, or a tool serving one or more of these questions.

## Q1. Which continuum sector produces the global matching H4 signal?

### What is already durable

- independent square-site blocks strongly disfavor global zero;
- the leading orientation dependence transfers as H4 across the tested Gaussian lineages;
- H12/H8 aliases tested by norm-5 are strongly disfavored relative to H4;
- N145->290 shows that the full thermal response is not one scalar multiplier;
- pure `P4[S']`, scalar-width and simple one-coordinate shell descriptions fail;
- the live response is low-rank/multicomponent rather than one clean scalar correction.

### Live field/mechanism classes

1. **Thermal Q4 / inherited Jordan candidate**

   The exact Virasoro/LCFT construction gives an `x=21/4`, spin-4 rank-2 candidate and exact torus-shape relations. This is a real representation-theory bridge, but lattice overlap remains unproved.

2. **Lower four-leg spin-4 competitor**

   The exact loop spectrum contains `V_(2,+/-2)` with `x=17/4`, spin 4. It is lower-dimensional and naturally related to four-leg geometry. Its physical `Q=1` multiplicity/representation and overlap with the global matching observable are active questions.

3. **Ordinary low-rank mixing**

   Two or more diagonalizable finite-size fields can mimic part of a log/Jordan sequence over the current size range. This remains a first-class competitor.

4. **Cover/defect/topological memory**

   Smith/deck/combinatorial-map information may create a small state that is not one local bulk field. The exact deck-character selection rule makes this possibility sharply testable with the correct charged/twisted observables.

5. **Higher-rank structure**

   Keep this live only where a basis-invariant rank-2 description fails or higher-point/second-derivative observables require it. Do not use rank 3 as a free rescue of one-insertion data.

### Highest-value discriminators

- Potts/projector/multiplicity and charged-insertion work for the `x=17/4` competitor;
- covariance-aware minimal realization before naming latent directions;
- modulus/shape fingerprints rather than another radial exponent;
- exact scale/cover composition under norm 4;
- independent positive controls where the expected field structure is already known or tunable.

## Q2. Why does matching select that sector?

This is the main conceptual bridge still missing.

### Exact lattice side

The project has strong finite-volume semantics:

- Sykes-Essam / Mertens-Ziff matching structure;
- exact wrapping-channel maps;
- finite Russo/pivotal identities;
- exact cover/deck arithmetic and selection rules;
- configuration/topological side programs.

### Missing continuum side

The empirical S/D or matching-even/odd decomposition does **not** by itself prove a local CFT involution with `+/-` eigenoperators.

The relevant routes are:

- **#61:** matching action on the RG tangent space;
- **#114:** FK/Potts `Q -> 1` topological-sector definition of the matching observable;
- **#233:** local defect/interface or derivative-defect construction, including a useful obstruction if no bounded-locality object exists;
- **#120:** transfer-matrix operator spectroscopy and matrix-element selection rules.

A partial formula, exact zero, small-width intertwiner, or precise no-go can all advance this question. This work should run alongside numerical discrimination rather than act as a prerequisite.

## Q3. Why is the square-site threshold that number?

The operator program and the threshold-value program should stay connected but distinct.

CFT/operator identification explains critical corrections and universal structure. It does not automatically determine the microscopic location

```text
p_c(square site) = 0.592746...
```

The direct threshold-origin line includes:

- correlated-hyperedge/self-dual or Yang-Baxter-compatible embeddings (#123);
- exact decorated-cell / critical-polynomial structure (#3/#13);
- finite matching-polynomial factors, recurrences and topology (#17/#29 and related exact work);
- post-leading annihilator structure (#47) when used to identify the microscopic correction hierarchy rather than only improve decimals.

A rigorous obstruction to a finite local exact mechanism is also useful: it would explain why the threshold is structurally harder than self-dual neighbors.

## Independent controls: the project should use them more aggressively

Controls are the fastest way to tell whether the analysis machinery recognizes known physics rather than merely organizing the target data.

### Triangular-site logarithmic pair

#246/#234 provides an externally defined percolation energy/logarithmic-pair construction. This is the strongest calibration for Jordan/radial/two-cutoff diagnostics.

### Exact C4 controls

- square-bond `p_c=1/2` duality control (#42);
- C4 self-matching site control (#44/#155).

These separate generic square anisotropy from matching/duality-odd structure.

### Tunable exactly-critical anisotropy

#106 is a high-value improved-action laboratory. A controllable zero crossing of `A4(lambda)` would test the H4 pipeline and expose subleading sectors without relying on the target lattice alone.

### C6/E6 mirror program

#165 tests whether the square-lattice `C4 -> E4/Gaussian` picture has a genuine symmetry/modular analogue or is model-specific.

## Data/analysis infrastructure

### Threshold-rank and thermal coordinates

The threshold-rank archive supports values, roots, derivatives, Krawtchouk/Hermite modes, quantiles and covariance-aware full-curve analysis. Different projections of one histogram are coordinates, not independent evidence.

### Gaussian-cover category

Norm-2, norm-5, norm-4 and norm-10 covers provide scale, phase, Smith/deck and composition information. The exact theorem that unmarked deck-invariant global observables have zero linear response to nontrivial deck characters should shape the observable design, not shut down the charged-sector program.

### Local pivotal / annulus sector

Russo connects slope to pivotal mass, and marked four-arm H4 observables are measurable. Current small-torus/local rows do not support a one-coordinate shell law; this motivates better state-space/radial diagnostics and independent controls, not a blanket rejection of local mechanisms.

### Primitive square-bond sector

The primitive square-bond H4/KdV response remains scientifically separate from the square-site thermal H4 line. Its next value is in shape/module identification, not another same-purpose sign repetition.

## Current decision frontier

```text
representation/selection:  x=17/4 competitor vs global Q4 candidate
state dimension:            Jordan vs ordinary low-rank vs cover-enriched memory
positive controls:          known log pair / self-matching / tunable anisotropy
orthogonal new evidence:    same-N coalescence / norm-4 / modulus shape
observable bridge:          matching -> RG/FK/defect/TM selection rule
threshold origin:           self-dual/hyperedge/exact microscopic mechanism
```

See [`NEXT-TARGETS.md`](NEXT-TARGETS.md) for the default attention order.

## Priority philosophy

A priority reduction means only that the next unit of effort is expected to change the research picture less. It does not close a branch, reject a mechanism, or prevent cheap exploratory work.

Raise a task whenever it gains a new independent observable, a stronger adversarial prediction, a better positive control, or a much cheaper path to the same scientific decision.