# Programs H + I: memory vs branching, and the successor-moment algebra of depth

## Program H — bounded one-branch memory, unbounded branching complexity (exact)

The #435 pair already carries the cleanest statement, now sharpened against
the #549 hierarchy:

* **FACT (#435).**  Two N16 torus configurations have the *identical complete
  unbranched survival vector* `(1, 7/8, 9/14, 5/14, 4/35, 0, 0, 0, 0)`, i.e.
  identical law of every one-branch future rank trace.  In the language of
  Program A this is `x ~_{L_0} y` — the unbranched predictive class is one.
* **FACT (#549).**  On the parallel-gadget version of the same pair, the
  `L_d` hierarchy splits that single `L_0` class into exactly `min(d+1, k+1)`
  classes, reaching `k+1` at full depth.

So the single `L_0` class contains `k+1` branching classes — an exact,
repository-native instance of

```
bounded (in fact unit) unbranched memory/state  +  unbounded branching
predictive complexity.
```

The reason is conceptual and now precise: shared-prefix branching probes the
correlation of *counterfactual futures conditioned on the same hidden
microstate*, whereas a one-branch marginal time series integrates that
correlation out.  #435's non-Markovianity witness (the survival-signature
process is not Markov under the actual permutation law) is the same phenomenon
read from the marginal side: the marginal is not the state.

## Program I — successor-hazard moments and what depth does to them

**FACT (algebra, exact).**  A single #549 fork on a group of `g` gadgets has
probability

```
F(g, a_g) = [ sum over the group's vertices of (7g - x_v)^2 ] / [8g(8g-1)^2],
```

which — because `sum_v x_v` is *equal* across the two gadget types (13) — is
an affine function of the **second** successor-hazard moment `sum_v x_v^2`
(29 for A, 25 for B).  Depth-1 therefore reads only the second moment.

**What higher depth does (the subtle point).**  A grouped experiment of depth
`j` is a product of `j` independent single forks, so it does **not** need
third/fourth/… successor-hazard moments; it is the `j`-th tensor power of the
*second*-moment difference.  Projected onto the class coordinate this is
precisely why `P_j(a)` has degree `j` in `a`: the algebraic degree grows, the
moment order stays at 2.

So the correct "reliability jet" statement is:

```
branching depth d  <=>  degree d in the second-moment tensor power,
                       not  order d in the hazard-moment hierarchy.
```

This answers Program I's question ("which polynomial-derivative hierarchy
corresponds exactly to branching depth") for the parallel family, and it
predicts the analogous tensor structure to look for in a general terminal
network: depth is an exponent in a commutative (multiplication of independent
forks) algebra, not a derivative order.

Files: `scripts/rank_lower_bound.py`, `results/rank-vs-depth/latest.json`.
