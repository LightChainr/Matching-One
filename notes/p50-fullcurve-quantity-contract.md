# P50 full-curve quantity-family contract

This contract precedes a typed operational wrapper for
`score_p50_fullcurve_n290.py`. It does not modify or replay the score.

Four exact topology anchors are registered:

- the raw matching-even orientation contrast for the three `X_even` levels;
- the raw matching-function value underlying the intrinsic center, mean slope,
  and implicit roots;
- angular-normalized matching-even `P4_S`;
- angular-normalized matching-odd `P4_D`.

The intrinsic levels, mean derivative, signed root gap, P4 derivatives,
finite-size ratios, and model corrections are response/model coordinates.
They are deliberately not added to `ObservableDescriptor` as if they were
new exact topology identities.

The contract also freezes N145→N290 representation order, both positive stored
lineage signs, independent random streams, the nine-feature order, the five
frozen scoring stages, and the residual-covariance rule
`Cov_child + ratio^2 Cov_parent`.

All nine features at one size reuse the same two-orientation histogram block.
They therefore remain correlated views, not nine independent evidence rows.
The numerical migration and any scientific interpretation remain separate
work; Issue #146 stays open.
