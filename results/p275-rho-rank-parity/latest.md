# Issue #275 rho-child rank-parity C3 split

## Decision

**RANK_EVEN_ONLY_C3_CORRECTION_AT_CURRENT_PRECISION.**  In the same 100 aligned
rho-child batches, the matching-odd `A_top=P2-P0` nontrivial C3 coordinate is
compatible with zero, while the Alexander-even `E_top=P0+P2` coordinate is strongly
resolved.  They are correlated transforms of one stream, not independent votes.

```text
A_top r1 = +4.04166666667e-05 -3.18986023727e-05 i
Cov(A_top r1) = [[2.50507488075e-08, -2.328295569e-09],
                 [-2.328295569e-09, 2.52022430556e-08]]
chi2=0.0969051/2, p=0.952703

E_top r1 = -0.00174463073628 -4.71983845063e-05 i
chi2=142.199/2, p=1.32382e-31
```

At Q=1 the continuum `A_top` baseline is exactly zero for every modulus because
Arguin gives `P2=P0`; no fitted continuum subtraction is used for that row.  The result
is finite square-bond fixed-p rank parity.  It does not supply a B-source insertion,
pooled-root p-jet, or original-U normalizer.

## Immediate use

The underlying archive retains rank0 and rank2 counts for all three children in every
batch, so its six-coordinate covariance is already sufficient to score any frozen
rank-restricted theory vector.  The missing object is the typed observer/source
transport column, not more rho-child sampling.

## Boundaries

- `A_top_and_E_top_are_same_stream_correlated_observer_transforms_not_independent_votes`
- `this_is_fixed_p_square_bond_rho_child_rank_parity_not_a_B_source_response`
- `no_pooled_root_p_jet_slope_or_original_U_normalizer`
- `a_zero_A_top_character_does_not_identify_the_nonzero_E_top_character_as_a_continuum_energy_field`
- `zero_new_Monte_Carlo_GPU_cloud_or_server`
