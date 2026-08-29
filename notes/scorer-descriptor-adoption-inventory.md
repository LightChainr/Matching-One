# Scorer descriptor-adoption inventory

At main commit `11ba86892fabd3349e322d198f23ef24d8cbe828`, the exact
`scripts/*score*.py` corpus contains 35 files. Python-AST inspection finds four
direct typed entrypoints importing both `ObservableDescriptor` and
`map_observable` from `wrapping_channels`. Three frozen kernels are covered by
those entrypoints through explicit semantic-manifest wrapper relationships.
The remaining 28 files are outside a registered typed path.

This closes an inventory gap: the repository now has a deterministic list,
including Git blob identities, rather than an informal claim that all scorers
have or have not adopted descriptors. The audit fails closed if a direct typed
import is added or removed without updating the manifest, or if a declared
wrapped kernel disappears.

`outside_registered_typed_path` is deliberately a triage label. Some of those
28 scripts may not compare channel-bearing quantities, some may be historical,
and some may need a future typed wrapper. Static membership alone does not prove
a semantic bug, data-provenance failure, or need for migration. Each candidate
must be reviewed before modification; Issue #146 therefore remains open.
