# Scorer descriptor-adoption inventory

At main commit `11ba86892fabd3349e322d198f23ef24d8cbe828`, the exact
`scripts/*score*.py` corpus contains 35 files. Python-AST inspection finds four
direct typed entrypoints importing both `ObservableDescriptor` and
`map_observable` from `wrapping_channels`. Three frozen kernels are covered by
those entrypoints through explicit semantic-manifest wrapper relationships.
One path, `kappa3_half_score.py`, is classified as a generic utility for which
a descriptor is not applicable. One path, `score_angular_root_amplitude.py`,
is confirmed channel-bearing and requires a typed migration. The remaining 26 files are outside a
registered typed path.

This closes an inventory gap: the repository now has a deterministic list,
including Git blob identities, rather than an informal claim that all scorers
have or have not adopted descriptors. The audit fails closed if a direct typed
import is added or removed without updating the manifest, or if a declared
wrapped kernel disappears.

`outside_registered_typed_path` is deliberately a triage label. Some of those
26 scripts may not compare channel-bearing quantities, some may be historical,
and some may need a future typed wrapper. Static membership alone does not prove
a semantic bug, data-provenance failure, or need for migration. Each candidate
must be reviewed before modification; Issue #146 therefore remains open.

The `kappa3_half_score.py` exception is narrow. It defines only exact Bernoulli
likelihood-score polynomials at `p=1/2` and aggregates an opaque caller-supplied
`D(C)`; it intentionally does not define a wrapping, homology, or other
topology observable. Geometry-specific callers still must type their own
`D(C)` before comparison.

`score_angular_root_amplitude.py` is not eligible for the generic-utility
exception. Its current blob reconstructs rank-2 cross `K_minus/K_plus`
observables, forms signed orientation differences and root gaps, normalizes by
signed `DeltaCos4`, and scores `A_p` across sizes. A later migration must decide
explicitly whether the descriptor schema needs a root-location quantity; this
inventory result does not retrofit or reinterpret the completed P45 score.
