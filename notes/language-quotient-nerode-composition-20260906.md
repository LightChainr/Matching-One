# Programs A + D + G + J: language semantics, canonical quotient, Myhill–Nerode, and compositionality

This note states what is now precise and what remains open for the four
structural programs.  Evidence tags are enforced throughout.

## Program A — a compositional experiment language (precise for the parallel family)

The parallel-gadget experiments admit the following exact semantics, which is
the cleanest special case of a rooted experiment tree:

```
experiment E ::= 1 (trivial terminal)
              |  fork-group g  (one #549 fork on a group of g gadgets)
              |  E1 . E2       (disjoint union of future regions: independent
                                product of success probabilities)
```

Two hidden networks agree through branching depth `d` iff every experiment
built from at most `d` fork-groups (products included) gives the same success
probability.  This defines `L_0 ⊆ L_1 ⊆ ... ⊆ L_k` and, by the rank theorem,

```
~_{L_0}  strictly coarser than  ~_{L_1}  strictly coarser than ... ~_{L_k},
```

with class counts `1 = r_0 < r_1 = 2 < ... < r_k = k+1`.  The hierarchy is
**strict at every step and stabilises at `L_k`** on this finite family — a
fact, not a conjecture, for the parallel case.

**OPEN (Program A, general case).**  The tree semantics with *sequential*
shared updates inside one group (non-product trees) is not yet reduced to a
closed algebra; doing so is the prerequisite for a general Myhill–Nerode
statement.

## Program D — canonical quotient of the cut network

**FACT (#550).**  The frozen summary `{S(z), n, H2, b2, radius-1 terminal
neighbourhood}` is *not* sufficient for a depth-2 compositional language: an
explicit planar two-terminal `n = 7` witness has identical summary but
different depth-2 outcome probabilities (`937/1050` vs `313/350`, gap `1/525`).
This is a summary-insufficiency certificate, not yet a minimality theorem.

**OPEN — coarsest predictive quotient.**  For each candidate summary `S`
(reliability polynomial; all vertex-deletion reliability polynomials; SPQR;
block-cut tree; Tutte-like terminal polynomial; boundary response), the single
test is `S(G)=S(H) => G ~_L H?`.  Known: `S(z)` alone is **false** (#435);
the #550 frozen summary is **false** at depth 2.  Conjectured simplifying
classes where a true `S` may exist: series-parallel, outerplanar, bounded
terminal pathwidth, cactus/block networks.  This is the main remaining
theorem target of the program.

## Program G — Myhill–Nerode / weighted-tree-automaton correspondence

**CONJECTURE (standard, to be verified against the tree semantics).**  Let the
branching experiment semantics be made compositional (Program A).  Then, for
the associated weighted tree language, finite branching-Hankel rank is
equivalent to the existence of a finite-dimensional linear predictive
representation, in exact analogy with ordinary weighted word languages
(Hankel rank = minimal linear-realisation dimension) and weighted tree
automata (the clone/fork is the tree-arity constructor).  Positivity then
selects the nonnegative/positive-realisation subclass, which is the content of
Program E.

**FACT (direction already established).**  On the parallel family the
branching-Hankel rank `r_d = d+1` is exactly the dimension of the minimal
linear predictive representation of the depth-`d` language (`span{1,a,...,a^d}`),
so the correspondence holds on the only family where the semantics is currently
closed.

## Program J — compositionality / syntactic congruence

**FACT (parallel composition is a congruence).**  For future-vertex-disjoint
sums, success probabilities multiply, so if `G ~_{L_d} H` then
`G ⊔ K ~_{L_d} H ⊔ K` for every disjoint `K` — disjoint union is a congruence
for every depth.  This is why the parallel family closes cleanly.

**OPEN.**  Series composition, terminal gluing, and vertex substitution are
not yet shown to preserve `~_L`; if they fail, the missing state is exactly
the counterexample the program asks for.  The long-term "syntactic algebra of
rank-one continuation networks" should be built from whichever of these
operations is a congruence.

Files: `scripts/rank_lower_bound.py`, `results/rank-vs-depth/latest.json`;
repository objects #435/#549/#550.
