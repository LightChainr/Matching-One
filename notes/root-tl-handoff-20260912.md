# Completed-analysis handoff for the existing issue channels

2026-09-12. Prepared comments, NOT posted by this delivery. No new issue is needed
for the work completed here. Preserve earlier artifacts and author review holds.

## #715 / #713: replace the negative process-identification paragraph

The asserted P398/TL distinction is false after an explicit change of state representation.
For each w-point noncrossing block, pair the right endpoint of each point with the left
endpoint of its cyclic successor, using 2w endpoints. This bijection Phi obeys
`Phi D_i=e_(2i) Phi`, `Phi J_i=e_(2i+1) Phi`, with loop multiplier one.
P398 G is exactly the periodic identified-connectivity (IC), zero-defect O(1) TL
stochastic generator on 2w endpoints (transpose the standard column-state convention).
The DC/through-line representations are NOT interchangeable with this one.

See `notes/p398-is-periodic-tl-20260912.md`. The all-width proof is local reconnection;
31,040 operations at w=2..8 and the entire archived w=4/5 G/H matrices were checked.
The half-step endpoint rotation also exchanges J/D and conjugates eta to -eta,
with necessary source/readout changes. Cantini--Sportiello's established RS theorem
supplies the eta=0 FPL-pushforward stationary law. The particular #709 pulse and
observer certificate remains a model-specific result; this correction neither
invalidates it nor certifies its novelty. No extra width table is commissioned.

Retain #715's realization-theory/Volterra references. Revise only the process map and
its dependent framing. For deeper retrieval compare the exact S/F/H and IC boundary
condition, not labels such as partitions vs matchings or join-only vs join/detach.
This is the continuation of #713, not a second broad prior-art ticket.

## #276 / #613 / #716: distinguish root consistency from full-law concentration

For square-site axis rectangles m>=w>=2, the matching root satisfies
`lim_(w->infinity) sup_(m>=w) |p_(w,m)-p_c|=0`.
The proof uses site one-arm exponential decay on NN and its matching partner,
Harris association, disjoint slabs, and existing digital-Alexander rank duality.
The key bounds are `P0>=(1-a_R)^(wm)` and
`P2<=(w*a_R)^floor(m/(R+1))`, with R=floor(w/8).
They compare exponential RATES even when both P0 and P2 are tiny.

The same claim holds for compact interior quantiles of `H=P2/(P0+P2)`;
H is not the birth-mixture CDF F. On w=j,m=ceil(exp(j^2)) the root still tends
to pc while the F distribution tends to half mass at each endpoint.
Thus #613's full-law condition stays intact, but it was unnecessarily strong
for the root alone on rectangles. No claim for arbitrary tilted quotients,
no numerical pc, no width-uniform critical exponent, and no all-width pTL map.

See `notes/aspect-uniform-balance-consistency-20260912.md` for the proof and
`rectangular-rank-odds-controls.json` for finite controls. Re-running thin lengths
or the old #613 root proof is not the next task. Literature should compare this
precise all-aspect conditional-ratio theorem with published estimator consistency.
This retrieval can stay in the existing #276/#613/#716 discussion.

## #594 / #598 / #600 / #610: use the correctly identified process

The original reflection quotient and first-order parity rule are unaffected.
They now live inside a named periodic IC TL process. The local H from #709 is
`e_1-e_(2w-3)` in the same operator algebra, although outside span{J,D}.
Do not call P398 a new unknown process; do not call it square-site percolation.
The half-step complement is NOT generally an involution and the original
three-readout dictionary is NOT invariant under it. It therefore supplies
transport identities, not an automatic larger quotient for the unchanged task.

The stationary FPL theorem gives exact baseline benchmarks at all widths without
a new enumeration program. Dynamic inserted matrix elements matching #709 are
still a distinct question for the ongoing primary retrieval. No new GPU or
large-CPU allocation follows from this handoff.
