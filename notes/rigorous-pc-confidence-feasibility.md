# Feasibility of certified finite-size bounds for square-site percolation

Status: Phase 0 method and cost gate for Issue 112. No new threshold bound is claimed.

## Outcome

Classification: **methodologically interesting but weak**.

The Riordan--Walters construction applies directly to independent site percolation on the square
lattice. Its statistical certificate is simple, exact, and inexpensive. The hard part is geometric:
at a parameter close to `p_c`, the block scale must be large enough that a demanding cluster event
has probability above a universal one-dependent threshold. The repository does not yet implement
that open-boundary event, and Riordan--Walters already reported a 99.9999% square-site interval
`[0.5925,0.5930]` in 2007.

The recommendation is therefore not to launch production. First build and independently validate
the event evaluator, use exploratory samples only to freeze the geometry and parameters, and then
decide whether a fresh certification run would improve the historical interval enough to matter.

## The finite event and deterministic implication

Partition the plane into disjoint `s x s` cells `S_v`, indexed by `v in Z^2`. For a nearest-neighbor
block edge `e=uv`, let `R_e=S_u union S_v`, a `2s x s` rectangle. Declare `e` open when:

1. the open subgraph induced by `S_u` has a unique largest cluster;
2. the open subgraph induced by `S_v` has a unique largest cluster;
3. those two clusters belong to the same open cluster in `R_e`.

The event depends only on sites in `R_e`. Rectangles belonging to vertex-disjoint block edges are
disjoint, so their event states are independent. The renormalized bond process is therefore
1-independent. An infinite open block path glues the selected clusters into an infinite cluster in
the original site model.

Riordan--Walters used the deterministic statement that every 1-independent bond model on `Z^2`
with edge marginal at least `0.8639` percolates. Balister--Johnston--Savery--Scott later proved the
stronger bound

```text
p_max(Z^2) <= 0.8457.
```

Thus a strict certificate that the block-event probability exceeds `0.8457` implies that the
original parameter lies at or above its critical probability.

## Exact statistical certificate

Let `pi` be the unknown block-event probability. Run `N` independent block simulations and observe
`M` successes. Under the null `pi<=p0`, monotonicity gives the exact conservative p-value

```text
P[Bin(N,p0) >= M].
```

No asymptotic normal approximation is needed. The included oracle evaluates this rational exactly.

Following the 2007 allocation, allow at most three adaptive attempts for an upper bound and three for
a matching-lattice lower bound. Bonferroni allocation of family-wise error `10^-6` gives

```text
alpha_run = 10^-6 / 6 = 1/6000000.
```

At `N=400`:

| deterministic threshold | minimum successes | exact-tail decimal |
|---:|---:|---:|
| `0.8639` | `378` | `1.148990352894...e-7` |
| `0.8457` | `373` | `9.514515135926...e-8` |

The first line reproduces the published 2007 calculation. The modern deterministic constant lowers
the required successes by five, but 400 trials still have low power unless the true event probability
is high:

| true `pi` | power at `M>=373` |
|---:|---:|
| `0.90` | `0.01494` |
| `0.92` | `0.20566` |
| `0.94` | `0.77376` |
| `0.95` | `0.95201` |

The final block scale must therefore make the event probability roughly mid-0.9, just as in the
historical design. More trials improve statistical power, but do not remove the need for a large
geometric block.

## Two-sided square-site interval

An upper-bound run on square-site `Z^2` at parameter `p` yields `p_c<=p` when successful.

For the lower side, run the upper-bound construction on the square site-matching graph, obtained by
adding same-face adjacency (the diagonal shell). If it certifies `p_c(matching)<=q`, the site matching
identity gives

```text
p_c(square) >= 1-q.
```

The family-wise allocation must include both sides and every permitted final attempt. Exploratory
parameter selection may be adaptive, but final samples must be fresh and independent of exploration.

## Heuristic cost model, not part of the theorem

Riordan--Walters used the planning ansatz that the critical window scales as

```text
delta p ~ s^(-3/4).
```

If per-trial work scales with rectangle area, shrinking the interval width by a factor `r` suggests

```text
linear scale multiplier = r^(4/3),
area-work multiplier    = r^(8/3).
```

Relative to width `5e-4`:

| target width | linear multiplier | area-work multiplier |
|---:|---:|---:|
| `2.5e-4` | `2.520` | `6.350` |
| `1e-4` | `8.550` | `73.100` |
| `5e-5` | `21.544` | `464.159` |
| `1e-5` | `184.202` | `33930.220` |

These are heuristic planning ratios. They are neither rigorous complexity bounds nor evidence about
the true finite-size constants.

## Deterministic versus statistical rigor

The implication from the block event probability to a critical-probability bound is deterministic.
The Monte Carlo conclusion is a high-confidence random interval, conditional on independent trials
from genuine Bernoulli site fields. A pseudorandom generator can be reproducible, well tested, and
domain separated without being a mathematical source of genuine randomness.

Deterministic exact enumeration would sum over exponentially many configurations in `O(s^2)` sites;
a frontier method remains exponential in width `s`. Neither is credible at the scales needed for a
narrow interval. The realistic deliverable is therefore a rigorous statistical confidence interval,
not a deterministic new theorem about the decimal value of `p_c`.

## Production gate

Before production:

1. implement the exact open-boundary unique-largest-cluster event;
2. test it on tiny rectangles by exhaustive enumeration;
3. test square/matching graph complement and orientation conventions;
4. use exploratory streams to freeze `s,p,N`, maximum attempts, and alpha allocation;
5. archive code, seeds/randomness provenance, compiler/runtime versions, and the raw success counts;
6. draw fresh final streams and evaluate only the predeclared exact binomial test.

Primary sources: Riordan--Walters, `arXiv:math/0702232`; Balister--Johnston--Savery--Scott,
`arXiv:2206.12335`.
