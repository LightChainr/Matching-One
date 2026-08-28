# A no-go lemma for exact finite independent-bond replacements of a Bernoulli site

Status: elementary theorem candidate. The statement is intentionally narrow; it rules out one mechanism, not exact solvability in general.

## Motivation

A single occupied square-lattice site acts as a four-terminal connectivity device. Ignoring connections through the rest of the lattice, its local terminal partition has only two possibilities:

- site vacant: the four incident terminals are mutually disconnected through the site;
- site occupied: all four incident terminals are connected through the site.

Thus the local partition law has support only on

\[
1|2|3|4
\quad\text{and}\quad
1234.
\]

It assigns exactly zero probability to every partial partition such as

\[
12|3|4,\qquad 123|4,\qquad 12|34,
\]

for every `0 < p < 1`.

A natural exact-mapping idea is to replace the site by a finite graph of ordinary independent bonds whose edge probabilities depend on `p`. The following observation rules out the nondegenerate version of that idea.

## Lemma

Let `H=(V,E)` be a finite undirected graph with at least three distinct labelled terminals `T subset V`. Suppose every edge is independently open with probability

\[
0<q_e<1.
\]

Assume that the event "all terminals are connected" has positive probability.

Then at least one **proper nontrivial terminal partition** also has positive probability: there is positive probability that some two terminals are connected while not all terminals are connected.

Consequently, no such finite independent-bond gadget can have terminal-partition support only on

- all terminals separate, and
- all terminals connected.

In particular, with four distinct terminals it cannot exactly reproduce the local terminal-connectivity law of one Bernoulli site at a nontrivial occupation probability `0<p<1`.

## Proof

Because all terminals can be connected with positive probability, the underlying graph contains a connected subgraph spanning all terminals. Choose an inclusion-minimal tree `S` spanning all terminals; equivalently, prune every nonterminal leaf from any spanning tree until every remaining leaf is a terminal.

`S` has at least two edges. Remove any edge `e` of `S`. Since `S` is a minimal terminal-spanning tree, each of the two components of `S-e` contains at least one terminal.

Now consider the following exact bond configuration of the full finite graph `H`:

1. every edge of `S` except `e` is open;
2. `e` is closed;
3. every edge of `E \ S` is closed.

Its probability is

\[
\left(\prod_{f\in S\setminus\{e\}} q_f\right)
(1-q_e)
\left(\prod_{g\in E\setminus S}(1-q_g)\right)>0,
\]

because every `q_f` lies strictly between zero and one.

In this configuration the terminals split according to the two components of `S-e`, so not all terminals are connected. Since there are at least three terminals, at least one of the two components contains at least two terminals. Those terminals are connected by the open edges of that tree component.

Therefore a proper partial terminal connection occurs with positive probability. QED.

## Four-terminal corollary

For four distinct terminals, any finite nondegenerate independent-bond gadget with positive all-connected probability necessarily gives positive probability to at least one partition of type

- `3+1`, or
- `2+2`, or
- `2+1+1`.

A Bernoulli site gives all such partitions probability exactly zero. Hence the two local laws cannot be identical.

## Deterministic edges

Edges with probability exactly zero can be deleted and edges with probability exactly one can be contracted before applying the lemma.

For an exact site replacement, deterministic contractions are harmless only if they do not identify distinct terminals in the site-vacant state. After this reduction, if at least three terminal classes remain distinct and all remaining stochastic edges satisfy `0<q_e<1`, the obstruction applies.

This also covers edge probabilities `q_e(p)` depending on the site parameter: for any fixed `p` in an open interval on which the reduced gadget is nondegenerate, the argument is pointwise.

Pathological constructions that deliberately switch edges between probability zero and one as a discontinuous function of `p` are outside the intended ordinary-bond mapping class and would need to be treated separately.

## What this does not rule out

The lemma does **not** rule out:

1. correlated bonds — e.g. four spokes controlled by one shared Bernoulli variable reproduce the site law immediately;
2. hyperedges or explicit multi-terminal interactions;
3. Potts/multispin/checkerboard reformulations;
4. infinite gadgets or limits of finite gadgets;
5. transformations that preserve only the critical point rather than the full local connectivity law;
6. stochastic domination gadgets used for rigorous upper/lower bounds;
7. approximate gadgets whose partial-partition probabilities tend to zero along a sequence.

These exceptions are not technical annoyances; they identify the model classes in which an exact structural mapping must live.

## Research consequence

The automated gadget search should be split into two programs.

### Exact mapping search

Do **not** spend time enumerating ordinary finite independent-bond gadgets hoping to reproduce a Bernoulli site's four-terminal law exactly. Search instead among:

- correlated edge bundles;
- hyperedge / random-cluster cells;
- multispin Potts interactions;
- transformations preserving a critical manifold without preserving the full local law.

### Approximation / theorem search

Ordinary independent-bond gadgets remain useful for:

- stochastic domination;
- certified replacement bounds;
- matching a selected subset of terminal-connectivity probabilities;
- constructing exactly solvable approximants whose local law converges to the site law.

Here the unavoidable partial partitions become an approximation error to minimize or bound, rather than something to pretend can vanish exactly.

## Stronger quantitative question

The qualitative no-go lemma suggests a useful optimization problem.

Fix the desired all-connected probability `p` and all-separated probability `1-p`. Among finite independent-bond gadgets of bounded complexity, how small can the total probability mass on forbidden partial partitions be?

Define

\[
\epsilon(H,p)=1-P_H(1|2|3|4)-P_H(1234).
\]

For each complexity budget `m`, seek

\[
\epsilon_m(p)=\inf_{|E(H)|\le m}\epsilon(H,p)
\]

subject to a calibration such as `P_H(1234)=p`.

Questions:

- Is `epsilon_m(p)` bounded away from zero for fixed `m`?
- What is its asymptotic decay as `m -> infinity`?
- Can a sequence minimizing `epsilon_m` be embedded into a self-dual critical-manifold construction?
- Can lower bounds on `epsilon_m` be certified by semialgebraic/interval methods?

This turns the failed exact gadget idea into a quantitative approximation theory.
