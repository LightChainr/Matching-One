# Predictive class count is not response-matrix rank: the affine-response obstruction

**Scope.** Exact finite statements around the cut-network material #435/#491/#549
and the general experiment-language question #599 (program "Hankel bridge").
The only repository object used as a *black box* is the closed-form probability
published with #549; everything else is elementary linear algebra on response
matrices.

---

## 1. The question and the trap

#549 exhibits, for every future size `k`, one complete-survival equivalence
class that splits into at least `k+1` exact branching-predictive classes
(`a = 0..k` copies of an A-type gadget inside a parallel bank of `k`).  The
declared fork test (one shared update, then two independent clone
continuations) has the closed form

```
F_{k,a} = [343 k^3 - 182 k^2 + 25 k + 4a] / [8k (8k-1)^2].
```

#599 program 2 asks whether those `k+1` classes can be converted into a growing
*response-matrix rank*: find future tests `E_j` with `M_{a,j} = P(E_j
succeeds | class a)` and `rank M >= k+1`.

The trap is the word "classes": a partition into distinguishable classes does
not produce directions in a linear space.  The theorem below shows that in the
#549 family it *cannot*, for the declared single-round fork language, no matter
how many fork tests are packed into the language.

## 2. Elementary lemma: response rank is bounded by polynomial degree

Let states be indexed by a parameter `theta in Theta subset R`, and let every
test in a future language `L` have success probability

```
P(E | theta) in span{1, theta, ..., theta^d}
```

as a function of `theta`.  Then for any finite set of states `{theta_0, ...,
theta_k}` and any finite `L`,

```
rank( M_{a,j} = P(E_j | theta_a) )  <=  d + 1.
```

Proof: every column of `M` lies in the span of the columns
`(1, theta_a, theta_a^2, ..., theta_a^d)_a`, a space of dimension at most
`d+1`.

**Corollary (affine obstruction).**  For the #549 fork test `F_{k,a}` is
affine in `a` (coefficient `4 / [8k(8k-1)^2]`).  Hence *any* finite language
built from single-round shared-update double-clone fork tests has
`rank M <= 2` on the `k+1` predictive classes, while the class count is
`k+1`.  The separation

```
predictive class count   vs   response rank
        k + 1                        <= 2
```

is exact and can be made arbitrarily large.  Distinct predictive classes do
*not* imply linear independence, and in this family they provably cannot
produce it from the single-round language.

**What would produce rank growth.**  A language containing a test whose success
probability is genuinely of degree `>= 2` in the class coordinate.  Concretely,
if tests `E_j` give probabilities `theta^j` for `j = 0..k`, the response matrix
is the Vandermonde matrix on `{theta_0, ..., theta_k}` and has rank `k+1`
whenever the `theta_a` are distinct.  So the question "does branching rank
grow?" is really the question "does the future language contain interactions
whose probability is a high-degree function of the class coordinate?"  For
#549, the fork test reaching degree 1 in `a` is exactly the same mechanism that
makes its *derivative in `a`* positive (`1/[2k(8k-1)^2] > 0`): the language is
first-order in the class label by construction.

Related exact repository anchor: #401's step-3 formula
`q3 = 1 - C(t,3)/C(d,3) + [e(t-2) - w + z + c3]/C(d,3)` with `t = d - x`
contains a term cubic in `x` (through `C(d - x, 3)`); a depth-3 compositional
language therefore has the *capacity* for rank 4 in an `x`-indexed family.
That is the expected mechanism by which depth upgrades rank, matching #550's
"r=1 bounded summaries are insufficient for a frozen depth-2 compositional
language": the witness lives precisely in the gap between summary coordinates
(degree ≤ 1 in the neighbourhood data) and the depth-2 outcome.

## 3. Where the honest lower bound stands

For the #549 family and the declared fork protocol we cannot exhibit a growing
rank without inventing a new protocol, and we do not do that here.  The precise
open statements we leave (for a future round or for whoever owns the #549
protocol file) are:

1. In the *actual* #549 sampling protocol, is the two-round probability (two
   sequential shared-update double-clone forks on disjoint vertex sets) a
   polynomial of degree exactly 2 in `a`?  If yes, rank 3 is reachable with
   three suitable tests and the class-vs-rank gap closes to `k+1` vs 3.
2. What is the minimal future-language depth at which the response rank must
   exceed any fixed bound, as `k` grows?  (Open; the polynomial-degree lemma
   makes this equivalent to the minimal depth producing a degree-(k+1)
   probability in `a`.)

## 4. Boundary

This note is linear algebra applied to one published closed form and to #401's
exact step-3 identity.  It is not a claim about the physical cut-network
percolation threshold, and it does not assign the "polynomial degree" any
continuum meaning.

Demo script: `scripts/probe/rank_vs_classes_demo.py` (affine family with many
classes and rank ≤ 2; Vandermonde family with rank = k+1).
