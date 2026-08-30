# A two-carrier abstract theorem for fixed-line ULC

## The abstract class

A boundary-regular two-carrier system consists of a monotone three-rank Boolean system, its exact complement-dual carrier, a contiguous rank-one sector, and two aggregate moment axioms:

`BA: A_(k+1) sum b x >= (sum b)(sum x)`,

`TM: (N-k) A_k sum d x >= (N-k-1) I_k sum x`.

BA is nonnegative aggregate birth/exit association on the upper layer. TM says the upper edge-weighted exit hazard dominates the lower uniform exit hazard.

## Exact theorem

BA moves from the upper edge marginal to the upper uniform layer; TM moves from the lower uniform layer to the upper edge marginal. Hence exit hazard is nondecreasing. If BA and TM hold on both complementary carriers, complement duality also makes birth hazard nonincreasing, and the exact hazard-ratio identity proves ULC.

## Minimality of the new axioms

Exhaustive symbolic enumeration contains 5 admissible systems at N=2, 111 at N=3, and 7,076 at N=4. No N<=3 system violates BA, TM, exit-hazard monotonicity, or ULC. At N=4 the basic axioms are no longer enough.

BA cannot be dropped: the first independence witness has sector `[1, 6, 9, 10, 12, 14]`. At layer 1, TM passes `6>=4`, but BA fails `32<35` and exit hazard decreases by `-1/24`.
TM cannot be dropped: sector `[1, 2, 6, 10, 14]` has BA equality `4=4`, TM failure `12<16`, and exit-hazard decrement `-1/6`.
The first direct ULC counterexample is also N=4, sector `[1, 5, 9, 13, 14]`, with failing normalized layers `[{'layer': 2, 'q_previous': '1/4', 'q': '1/3', 'q_next': '1/2', 'margin': '-1/72', 'pass': False}]`.

## Topological status

All 984 existing torus carrier-layer pairs satisfy BA and TM exactly. Thus the new theorem captures every checked topology, but the derivation of BA and TM from homology geometry remains open. The useful proof targets are now integer moment inequalities, not a global stochastic coupling.
