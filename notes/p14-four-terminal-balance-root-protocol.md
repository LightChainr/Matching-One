# P14 bounded four-terminal balance-root screen

Protocol frozen before reading the 27-candidate score table.

## Input

- Parent: `exact/issue-14-bounded-terminal-reliability-corpus@0295d68`.
- Frozen corpus: `results/terminal-reliability/bounded-four-terminal-corpus.json`.
- One homogeneous independent bond probability `p` on every gadget edge.

## Primary scalar

For every candidate define the permutation-invariant necessary balance defect

```text
B_G(p) = P_G(0123; p) - P_G(0|1|2|3; p).
```

Enumerate every real root in the open unit interval with exact polynomial
coefficients and interval-certified numerical isolation. Rank the nearest root
to the descriptive square-site reference `0.5927460507896`. The ranking is an
exploration coordinate, not a critical-threshold estimator.

## Structural gates

Report, without optimizing them after the roots are seen:

1. bridge and articulation status of the finite graph;
2. terminal-automorphism orbit size and whether the terminals are transitive;
3. whether at least one cyclic terminal order has identically zero crossing
   partition probability. This is only a necessary outer-face-planarity gate;
4. edge count, internal degree and exact minimal polynomial of each root.

Primary attention goes to candidates that pass the no-bridge,
no-articulation, terminal-transitive and crossing-zero gates. If none pass,
that emptiness is the result; do not silently relax the gates.

## Interpretation boundary

`B_G(p)=0` alone is not a four-terminal self-duality theorem. The corpus has no
embedding, periodic tiling, stochastic-domination map or published comparison
baseline. No root from this screen is a rigorous percolation bound or a claim
about the square-site threshold. Its purpose is to decide whether the frozen
bounded family contains a structurally credible target for a later real
comparison construction.

## Pre-score algebra addendum

The certificate always reports the exact primitive square-free balance
polynomial. It labels a factor as the root's minimal polynomial only when an
exact rational factorization certificate is supplied; otherwise it makes no
irreducibility claim. This narrows the wording above without changing any
candidate, root, reference value or structural gate.
