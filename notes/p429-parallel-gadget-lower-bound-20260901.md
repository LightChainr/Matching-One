# Parallel-gadget lower bound for branching predictive complexity

Date: 2026-09-01  
Paper track: continuation representations / predictive-state nonclosure  
Depends on: PR #435 exact N16 witness, PR #491 cut-network theorem.

## Claim status

**EXACT_THEOREM in the finite two-terminal cut-network continuation class.**

This note strengthens the single N16 no-go witness into a parameterized family with an unbounded number of branching-distinguishable predictive classes inside one and the same complete-unbranched-survival class.

It adds no Monte Carlo and no new descriptor. The only numerical input is the exact PR #435 witness, already independently verified there.

---

## 1. Base gadgets

PR #435 gives two rank-one N16 microscopic states, call them `A` and `B`, with eight future vacancies each and identical safe-subset counts

```text
b = (1, 7, 18, 20, 8, 0, 0, 0, 0).
```

After cutting along the occupied essential cycle, PR #491 maps the continuation problem to a planar two-terminal vertex network in which ambient rank reaches two exactly when the two cut boundaries connect. The N16 branching witness is included in the PR #491 exact controls.

For either gadget define the safe-subset polynomial

```text
S(z) = 1 + 7 z + 18 z^2 + 20 z^3 + 8 z^4.
```

The first-step data from PR #435 are:

```text
one immediate absorbing vacancy per gadget,
seven safe first vacancies per gadget.
```

For a safe first vacancy `v`, let `x(v)` be the number of the remaining seven vacancies whose activation would then create terminal connection / rank two.

The exact successor distributions are

```text
A: x=1 occurs 3 times, x=2 occurs 2 times, x=3 occurs 2 times,
B: x=1 occurs 1 time,  x=2 occurs 6 times, x=3 occurs 0 times.
```

Hence both gadgets have the same first moment

```text
sum_safe x = 13,
```

but different second moments

```text
q_A := sum_safe x^2 = 3*1 + 2*4 + 2*9 = 29,
q_B := sum_safe x^2 = 1*1 + 6*4       = 25.
```

The equality of the first moments is exactly why the ordinary two-step unbranched survival probability agrees; the second-moment gap `q_A-q_B=4` is what the delayed fork sees.

---

## 2. Parallel composition

Fix `k>=1`. Take `k` pairwise future-vertex-disjoint copies of the cut-network gadgets and place them in parallel between the same two deterministic terminal boundaries. Each component is independently chosen to be type `A` or type `B`.

Let `a` be the number of type-A components, so `0<=a<=k`.

The components share only the deterministic terminal boundaries. Therefore the global network is safe exactly when every component is safe. A future subset of total size `m` is safe iff its restriction to each gadget is safe.

Consequently the safe-subset generating polynomial of the entire `k`-gadget network is

```text
S_k(z) = S(z)^k,
```

independent of `a` and independent of the A/B ordering.

Thus for every `m`,

```text
b_m^(k,a) = [z^m] S(z)^k,

s_m^(k,a) = b_m^(k,a) / binom(8k,m),
```

so all `k+1` composition classes `a=0,...,k` have exactly the same complete unbranched survival law.

Because rank/terminal connection is absorbing and monotone, this also fixes the entire ordinary unbranched exit-time law.

---

## 3. Delayed-fork observable

Use the same branching language as PR #435:

1. choose one of the `8k` future vertices uniformly and activate it;
2. if the resulting network is already connected, score zero;
3. otherwise clone the safe successor into two copies;
4. in each clone independently activate one uniformly chosen remaining vertex;
5. score one iff both clones remain safe.

Call the resulting probability `F_{k,a}`.

Suppose the common first activation is safe and occurred in gadget `j`, with local successor exit count `x`.

Each untouched gadget still has exactly one immediately absorbing vacancy. Therefore among the `8k-1` remaining vertices, the total number of one-step absorbing choices is

```text
(k-1) + x.
```

Hence each clone has exactly

```text
(8k-1) - [(k-1)+x] = 7k-x
```

safe second choices, and the conditional two-clone survival probability is

```text
((7k-x)/(8k-1))^2.
```

Averaging over the `8k` possible common first activations, with absorbing first activations contributing zero, gives

```text
F_{k,a}
 = 1/[8k(8k-1)^2]
   * sum_over_gadgets sum_over_safe_first_v (7k-x(v))^2.
```

For one gadget of type `T`,

```text
sum_safe (7k-x)^2
 = 7*(7k)^2 - 14k*sum_safe x + sum_safe x^2
 = 343 k^2 - 182 k + q_T.
```

Using `q_A=29` and `q_B=25`, the exact closed form is

```text
boxed:
F_{k,a}
 = [343 k^3 - 182 k^2 + 25 k + 4 a]
   / [8 k (8k-1)^2].
```

In particular,

```text
F_{k,a+1} - F_{k,a}
 = 4/[8k(8k-1)^2]
 = 1/[2k(8k-1)^2] > 0.
```

So `a=0,1,...,k` gives `k+1` distinct delayed-fork probabilities despite one identical complete unbranched survival law.

For `k=1`, this reduces exactly to the PR #435 values

```text
F_{1,0}=93/196,
F_{1,1}=95/196,
difference=1/98.
```

---

## 4. Theorem — unbounded splitting of one survival-law class

**Theorem.** For every `k>=1`, the parallel cut-network family above contains `k+1` states with identical complete unbranched survival law but pairwise distinct delayed-fork probabilities.

**Proof.** The common survival law follows from `S_k(z)=S(z)^k`. The fork probabilities are given by the boxed formula and differ strictly by `1/[2k(8k-1)^2]` as `a` increases by one. QED.

This is a growing-family statement, not another finite census.

---

## 5. Predictive-state lower bound

Any exact autonomous coarse state that is sufficient for the delayed-fork branching language must assign different states whenever the delayed-fork probabilities differ.

Therefore, within a single complete-survival equivalence class, the exact branching predictive partition has at least

```text
k+1
```

classes on a network with

```text
8k
```

future random vertices.

Equivalently,

```text
number of exact branching predictive classes >= n_future/8 + 1
```

on this family.

A finite-state exact representation therefore needs at least

```text
ceil(log2(k+1))
```

bits to distinguish these states.

This is a lower bound on **state cardinality / exact predictive equivalence classes**, not on Euclidean state dimension, continuum field count, or the number of real-valued scalar coordinates. One real number can encode arbitrarily many discrete labels; this theorem does not claim otherwise.

---

## 6. Embedded realization lemma

The construction is not restricted to an abstract product Markov chain.

Each PR #435 rank-one witness has, by PR #491, a planar two-terminal cut-network realization. Put `k` copies in disjoint strips of a cylinder, sharing only the two deterministic cut-boundary components. Parallel terminal connection occurs iff at least one strip connects.

Gluing the two cylinder boundaries back together produces an embedded torus graph with one deterministic old essential occupied cycle; the gadget interiors remain disjoint. The initial occupied subgraph has ambient rank one, and future rank reaches two exactly when the parallel cut network connects.

Thus the growing lower-bound family lies inside the same **finite embedded-graph / rank-one continuation category** as the PR #491 theorem.

This realization does **not** assert that every member is a nearest-neighbour square-site quotient or one of the repository's Gaussian HNFs.

---

## 7. Red-team checks

### Cross-gadget paths

Could a first update in one gadget and a second update in another jointly create a connection even though neither gadget connects individually?

No. The gadget interiors are disjoint and share only the terminal boundaries. Any L-to-R path lies entirely inside at least one gadget between its first departure from L and first arrival at R. Therefore the parallel union connects iff one component connects.

### Fixed-cardinality sampling

The argument does not replace the repository's without-replacement continuation law by independent Bernoulli sampling. The coefficient identity

```text
[z^m] S(z)^k
```

counts fixed-cardinality safe subsets exactly, and division by `binom(8k,m)` gives the actual uniform-`m` survival probability.

### Is the lower bound only another function of survival law?

No. Every member has exactly the same full survival polynomial `S(z)^k`. The separating observable uses a shared update followed by a product/fork, which probes a second moment of the successor hazard distribution that the original survival law does not determine.

### Does this prove cut-network minimality?

No. It proves only that the exact branching predictive partition has unbounded cardinality on one explicit family. It does not prove the full PR #491 network object is minimal.

---

## 8. Manuscript consequence

The cut-network/nonclosure paper can now make a quantitative growing-family statement:

> **There exist rank-one embedded continuation systems with linearly many exact branching-predictive classes inside a single complete-unbranched-survival class.**

Together with PR #491, the theorem pair becomes

```text
negative:
complete unbranched survival can collapse k+1 distinct branching states
into one class;

positive:
the cut-network representation remains update-closed in the declared scope.
```

The paper no longer needs an asymptotic interpretation of the finite N16/N17 census to claim unbounded predictive complexity.

---

## 9. Next theorem and stop rule

The next useful question is now **minimality up to network equivalence**, not another growth census.

A high-value target is to identify a canonical quotient of the cut network that preserves all unmarked shared-prefix/fork experiments, or prove that a natural proposed quotient (for example two-terminal reliability polynomial plus bounded local terminal data) still loses information.

**Stop rule:** do not build larger parallel products merely to increase the `k+1` count. The unbounded lower bound is already exact. Further work must either reduce the constructive state or prove a stronger minimality/no-quotient statement.
