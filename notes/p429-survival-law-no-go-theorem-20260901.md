# Survival-law no-go theorem for topological continuation states

Date: 2026-09-01  
Paper track: continuation representations / predictive-state nonclosure  
Primary sources: PR #435 (`ffd91ebd819f7893cbee84aeb3f40da14a700a7b`) and PR #491 (`ab90201e88409310632812727e0138c56b455644`).

## Claim status

**EXACT_THEOREM, finite embedded scope.**

This note extracts a manuscript-level no-go theorem from the already exact N16 branching witness in PR #435 and connects it to the rank-one cut-network representation of PR #491.

It introduces no new Monte Carlo, no new topology oracle, no new continuum claim, and no new independent evidence block. The exact numerical/topological witness remains the one independently verified in PR #435.

## Setup

Consider the finite occupation-growth process on a fixed embedded torus graph. From a current rank-one microscopic configuration `x`, future vacant sites are activated uniformly without replacement until ambient homology rank reaches two.

Let

```text
L(x) = (s_0(x), s_1(x), ..., s_d(x))
```

be the **complete unbranched survival law**, where `s_m(x)` is the probability that rank remains one after `m` uniformly chosen future activations.

A coarse state `Z(x)` is called **survival-law measurable** if there is a deterministic map `F` such that

```text
Z(x) = F(L(x)).
```

This class is deliberately very broad. It contains every scalar/vector descriptor computed only from the complete unbranched survival distribution: all survival moments, hazards, cumulative hazards, finite transforms, fitted parameters, and arbitrary nonlinear functions of the full survival vector.

Now define a delayed-fork experiment `B(x)`:

1. make one common uniformly random safe/absorbing activation from `x`;
2. if the successor remains rank one, clone that microscopic successor into two copies;
3. evolve each clone by one independent uniformly random activation;
4. record whether both clones remain rank one.

This is the smallest shared-prefix branching experiment used by the exact witness in PR #435.

## Exact witness imported from PR #435

On the `4 x 4` square torus (`N=16`), PR #435 gives two rank-one configurations `A` and `B` with the same primitive ambient line `(1,0)` and the same complete survival count vector

```text
b = (1, 7, 18, 20, 8, 0, 0, 0, 0),
```

hence the same complete survival law

```text
L(A) = L(B)
     = (1, 7/8, 9/14, 5/14, 4/35, 0, 0, 0, 0).
```

The same PR exhaustively verifies that all ordinary unbranched future rank traces have the same distribution.

But after one shared update and a one-step independent fork,

```text
B(A) = 95/196,
B(B) = 93/196,
B(A) - B(B) = 1/98 != 0.
```

The witness is exact; it is not a sampling fluctuation.

## Theorem 1 — complete unbranched survival is not branching-sufficient

**Theorem.** There is no function `G` such that

```text
B(x) = G(L(x))
```

for every rank-one state in any state class containing the two PR #435 witnesses `A` and `B`.

**Proof.** `L(A)=L(B)` exactly, while `B(A) != B(B)` exactly. Therefore any function of `L` must take the same value on `A` and `B`, but the delayed-fork probability does not. Contradiction. QED.

This is stronger than saying that a particular scalar summary fails. The entire complete unbranched survival distribution is insufficient for this branching prediction.

## Corollary 1 — no survival-law-measurable compression can close delayed branching

**Corollary.** Let `Z(x)=F(L(x))` be any survival-law-measurable coarse state, of arbitrary finite or infinite dimension. Then `Z` cannot be sufficient for the delayed-fork observable `B` on a state class containing `A` and `B`.

**Proof.** Since `L(A)=L(B)`, one has `Z(A)=Z(B)` for every deterministic `F`. If `B` were determined by `Z`, then `B(A)=B(B)`, contradicting the exact `1/98` gap. QED.

Consequences include immediate rejection, for branching closure, of any hierarchy that only adds more functions of the same complete unbranched survival law. In particular, adding further survival moments, fitted hazard coefficients, or exact transforms of `L` cannot repair this witness.

This does **not** reject those summaries for predicting the original unbranched rank trace; `L` is complete for that narrower target by construction.

## Theorem 2 — survival-equivalence is not strongly lumpable for the microscopic growth kernel

Let two microscopic rank-one states be equivalent when they have the same declared labels and the same complete survival law. Consider the ordinary one-site growth kernel with rank-two configurations sent to a cemetery state.

**Theorem.** On the N16 state space of PR #435, the partition induced by complete survival-law equivalence is not strongly lumpable.

**Proof.** PR #435 computes the coarsest strong Markov refinement and finds four survival classes that split at N16 (`210` survival classes versus `214` strong Markov classes). It also gives a direct witness: under the actual uniform-growth law, states/cohorts with the same recomputed survival signature have different probabilities of transitioning to a successor signature with three one-site exits; e.g. `1/6` versus `2/11`, a difference `1/66`. Thus transition probabilities to a target coarse class are not constant inside the survival-equivalence block. By the finite-chain strong lumpability criterion, the partition is not strongly lumpable. QED.

The imported Markov-chain theorem is classical: a partition is strongly lumpable iff the total transition probability from any two states in one source block to each target block is identical. The Matching-One contribution is the exact embedded-percolation witness, not the abstract lumpability criterion.

## Corollary 2 — linear trace realization and autonomous branching state are different problems

Any representation whose state is designed only to reproduce the family of unbranched survival functions

```text
P^m 1_alive
```

need not be closed under shared-prefix branching experiments, because branching introduces pointwise products of successor observables before further propagation.

For successor survival observables `f_i`, the relevant object includes

```text
P(f_i f_j),
```

whereas an unbranched linear realization only requires the linear Krylov family generated by repeated `P` action.

The exact PR #435 witness proves that these two closure questions differ in the physical occupation-growth process; this is not only a generic systems-theory possibility.

## Connection to the cut-network representation

PR #491 supplies a constructive state that is rich enough for the declared embedded rank-one continuation problem:

1. cut along a simple occupied essential cycle;
2. contract old occupied components;
3. retain the resulting two-terminal planar vertex network;
4. future ambient rank reaches two iff the two cut boundaries become connected.

Keeping the initial cut fixed and updating/contracting activated sites gives an update-closed network-valued continuation representation within the PR #491 scope.

The two results therefore fit together as a theorem pair:

```text
NO-GO:
complete unbranched survival law
    is not sufficient for delayed branching / autonomous closure

CONSTRUCTIVE:
rank-one cut network
    is update-closed for the declared unmarked continuation problem
```

This is a substantially sharper paper claim than a sequence of scalar descriptor failures.

## What this theorem does not prove

It does not prove:

- that the PR #491 cut-network representation is minimal;
- a lower bound on continuum field dimension;
- that every fixed-`k` scalar summary fails on every graph size;
- that survival-law equivalence fails lumpability for every HNF or every embedded graph;
- that the number of predictive states grows at a particular asymptotic rate;
- a CFT/Jordan/operator identification;
- a new population-level experimental effect from #334.

The exact quotient counts reported in PR #435

```text
N16: survival classes 210, strong Markov classes 214
N17: survival classes 346, strong Markov classes 390
```

are finite census results, not an asymptotic lower-bound theorem.

## Manuscript consequence

The paper can now make a clean, bounded statement:

> **We show that complete unbranched survival information does not determine branching continuation in a finite topological birth process, and we give an exact cut-network state that restores update closure in the rank-one embedded setting.**

The natural manuscript structure is:

1. ambient-homology birth process and continuation language;
2. complete survival law and exact N16 no-go theorem;
3. strong-lumpability failure under the actual growth kernel;
4. cut-network representation theorem;
5. pair-trigger bipartiteness and genuine triple cooperation;
6. exact real N425 mechanism examples as applications, not theorem evidence;
7. discussion of predictive equivalence, approximate compression, and limits of continuum interpretation.

## Highest-value next theorem

The next useful strengthening is **not** another survival descriptor.

Try to prove a genuine growth statement for predictive complexity, for example one of:

```text
(a) number of strong continuation-equivalence classes grows unboundedly with graph size;
(b) a declared family requires an unbounded number of distinguishable network states;
(c) no fixed finite family of named scalar graph invariants is update-closed on that family.
```

A valid theorem may be restricted to a deliberately constructed embedded family. It need not cover all square tori.

### Stop rule

If no growing family can be produced with an explicit injective/distinguishing invariant or exact branching experiment, stop at the present no-go + constructive representation pair. Do not respond by extending the scalar ladder (`H2`, `W2`, `c3`, `c4`, ...).

## Provenance and evidence role

- PR #435: exact/control; open PR / unmerged at the 2026-09-01 snapshot; exact witness and finite quotient census.
- PR #491: exact/control; open PR / unmerged at the 2026-09-01 snapshot; cut-network theorem and real-checkpoint deterministic reconstruction.
- PR #492: exact/control; open PR / unmerged; complementary dual-cycle certificates.
- Issue #334 prospective 600k-prefix intervention: primary experimental result, separate dependency family; useful motivation but not used to prove the theorems above.

No p-values or independent-evidence counts are created by this synthesis.
