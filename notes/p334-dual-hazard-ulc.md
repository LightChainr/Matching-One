# The dual-hazard reduction of fixed-line ULC

Let `F_k` be the size-`k` configurations whose homology image is one fixed
primitive line, and let `A_k=|F_k|`. For `S in F_k`, split its Boolean boundary:

- `b(S)`: occupied deletions that fall to rank zero;
- `d(S)=k-b(S)`: deletions that stay in the fixed line;
- `x(S)`: absent insertions that jump to rank two;
- `u(S)=N-k-x(S)`: insertions that stay in the fixed line.

After summing, write

`beta_k = sum_(S in F_k)b(S)/(k A_k)`,

`xi_k = sum_(S in F_k)x(S)/((N-k) A_k)`.

The internal edge count between `F_k` and `F_(k+1)` can be read from either
side. Therefore

`(N-k)A_k(1-xi_k)=(k+1)A_(k+1)(1-beta_(k+1))`,

and hence

`q_(k+1)/q_k=(1-xi_k)/(1-beta_(k+1))`.

This is the source/sink/current form of the ULC question. Alexander complement
duality sends a primal fixed-line state at layer `k` to the matching fixed-line
state at layer `N-k`, with `b_P=x_M` and `x_P=b_M`. Consequently

`beta_P(k)=xi_M(N-k)`, `xi_P(k)=beta_M(N-k)`.

We obtain an exact conditional lemma: **if the exit hazard `xi_k` is
nondecreasing on both complementary carriers, then `beta_k` is nonincreasing,
the adjacent q ratios are nonincreasing, and q is ULC.**

Rank monotonicity and order-convexity nearly prove the hypothesis. Along every
internal edge `S subset T`, every exit-pivotal absent site of `S` remains
exit-pivotal for `T`, while every birth-pivotal site of `T` was already
birth-pivotal for `S`. Thus normalized hazards are monotone under the
edge-coupled measure.

The exact obstruction is the change of measure from internal-edge weighting to
uniform layer weighting. If `h_x=x/(N-k)`, then

`xi_(k+1)-xi_k = Delta_edge`

`  + Cov_k(u,h_x)/E_k[u] - Cov_(k+1)(d,h_x)/E_(k+1)[d]`.

The edge slack `Delta_edge` is nonnegative by the pointwise nesting lemma. The
two covariance terms have no forced combined sign. In the existing atlas their
sum is negative in 484 of 984 primal/matching adjacent-layer pairs. The worst
case is `-1/12`; an edge slack of `5/12` still leaves uniform increment `1/3`.
This is the precise missing inequality, rather than a vague absence of a
matching.

The bounded result is strong: all 83 honest quotients through `N=12`, all 240
fixed lines on each carrier, and all 59,922 complementary rank-one states pass
the duality and uniform-hazard gates. But it remains finite evidence. A full
proof now needs an aggregate two-step path or boundary double count showing

`Delta_edge + degree_bias >= 0`

without requiring a pointwise layer matching, which the `N=11` dead-end witness
already forbids.
