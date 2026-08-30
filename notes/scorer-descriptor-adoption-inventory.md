# Scorer descriptor-adoption inventory

In this checked snapshot, the exact `scripts/*score*.py` corpus contains 45
files. Python-AST inspection finds fourteen
direct typed entrypoints importing both `ObservableDescriptor` and
`map_observable` from `wrapping_channels`. Thirteen frozen kernels are covered by
those entrypoints through explicit semantic-manifest wrapper relationships.
One path, `kappa3_half_score.py`, is classified as a generic utility for which
a descriptor is not applicable. Seventeen paths are confirmed channel-bearing and
require typed migrations: `score_angular_root_amplitude.py`,
`score_axis_pair_annihilator.py`, its stable entrypoint,
`score_issue43_secondary.py`,
`score_norm4_production.py`,
`score_norm4_thermal_jet.py`, `score_norm5_thermal_jet.py`,
`score_p159_pell_hex_filter.py`,
`score_p231_vacuum_kdv_sector.py`,
`score_p49_fullcurve_doubling.py`,
`score_p50_fullcurve_n290.py`,
`score_prequential_evidence.py`,
`score_rank_gap_boundary_targets.py`, `score_threshold_rank_root_doubling.py`,
`threshold_score_modes.py`,
`score_v14_fixedp_scalar_projector.py`, and
`score_v14_scalar_root_projector.py`. No scorer remains outside a registered
audit class.

This closes an inventory gap: the repository now has a deterministic list,
including Git blob identities, rather than an informal claim that all scorers
have or have not adopted descriptors. The audit fails closed if a direct typed
import is added or removed without updating the manifest, or if a declared
wrapped kernel disappears.

The manifest explicitly excludes `scorer_descriptor_adoption_audit.py` from
its own `scripts/*score*.py` glob. The exclusion is checked against the live
corpus, and the regression test now recomputes the complete audit before
comparing it with the checked result. This prevents the meta-tool from being
silently counted as an unclassified scorer while keeping the scorer partition
unchanged.

`outside_registered_typed_path` remains a deliberate fail-closed triage label,
but no current scorer has that status. Static membership alone would not prove
a semantic bug, data-provenance failure, or need for migration; any future
corpus addition must still be reviewed before modification. Issue #146 remains
open because inventory classification is not the same as implementing the 17
required typed migrations.

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

`score_prequential_evidence.py` is also migration-required. Its scored-block
gate currently compares free-form `channel.source` and `channel.target` strings
and treats any truthy `exact_map` field as sufficient. A typed migration must
parse source/target descriptors and verify the registered affine transform,
while preserving historical ledger rows and chronology.

`score_c4_self_matching_n26.py` is now a covered frozen kernel. Its typed
entrypoint validates an exact identity map for the frozen `either` matching-odd
value before replaying the kernel and annotates the output with both
descriptors and the applied transform. The semantic gate binds the prediction
file hash and declared channel. The fact that all five wrapping-channel
Bernstein vectors coincide on this finite control remains only a reported
numerical identity, not permission to exchange channel labels. The frozen
hypotheses, scoring order, enumeration artifacts, and
stop-without-generalized-fit rule are unchanged.

`score_c4_tangent_orthogonal_holdout.py` is now a covered frozen kernel. Its
typed entrypoint validates an exact identity map for the same `cross/primal`
event at N=130 and N=170 before replaying the kernel. The semantic gate also
freezes ordered `(t, lambda)` response coordinates, the source-frozen
`lambda/t` projection, synchronized delete-one batches, and the fixed
`N^(3/8)` thermal map. These response coordinates are not reinterpreted as
topology channels or continuum spin eigenfields; frozen numerics and the
interpretation rule are unchanged.

`score_matching_odd_synthesis.py` is now a covered frozen kernel. Its typed
entrypoint validates the registered exact identity `D_either=D_cross` for both
frozen matching-odd orientation-contrast blocks before synthesis. It records
the two source/target descriptors and applied maps without changing the frozen
block selection, already-primary scores, distinct raw-data-group requirement,
block-diagonal calculation, or the rule that derived output cannot become a
new primary evidence row.

`score_rank_gap_boundary_targets.py` is migration-required. It scores the
orientation-pooled paired observable `G=K_plus-K_minus` in rank units from a
frozen source fit against target sizes, while binding its meaning only through
ordered representation arrays and metadata. A typed migration must represent
the paired rank-gap quantity, rank units, orientation pooling, and cross-size
identity. It must preserve the fixed `5/8` exponent, source/target chronology,
production-metadata gate, and covariance calculation.

`score_threshold_rank_root_doubling.py` is migration-required. It tests a
fixed `-1/4` doubling ratio for two Gaussian threshold-rank root-gap lineages,
reversing the stored child sign to follow multiplication by `1+i`. A typed
migration must represent the root-gap quantity, genealogy, stored-versus-lineage
orientation order, sign map, and cross-size relation. It must preserve both the
full-covariance score and the diagonal-only diagnostic.

`score_issue50_n290.py` is now a covered frozen kernel. Its typed entrypoint
validates an exact identity map for the fixed-p `either/matching-odd`
first-minus-second orientation contrast before scoring, and records both
descriptors plus the applied transform. The semantic gate freezes N=290 and
ordered Gaussian lineages `(13,11)` then `(17,1)`. The prospective target, run
provenance, counter-range validation, sampling-error calculation, zero-control
evidence reuse, and numerical outputs remain unchanged.

`score_p231_vacuum_kdv_sector.py` is migration-required. It concatenates
ordered `C_nontrivial_real`, `Q_reflection_null`, and `S_scalar` coordinates
for two designs, then compares them with a theory vector normalized per unit
`g4`. A typed migration must bind every sector coordinate, design order, and
theory normalization. It must preserve the block covariance, one-amplitude
GLS, non-scalar diagnostic, and the rule that the retrospective reuse is not
new independent evidence.

`score_p48_sprime_frozen.py` is now a covered frozen kernel with its own typed
entrypoint and semantic manifest. The gate proves that this chronological
contract uses the same cross/even, angular-normalized `P4_S_prime` descriptor
as the separate prospective path, and validates an exact identity map before
scoring. Target sizes, leading power, four-model chronological order,
source-plus-target covariance, target independence, and no-refit semantics are
unchanged.

`score_p49_fullcurve_doubling.py` is migration-required. One entrypoint emits
matching tails and slopes, thermal `X_even/X_odd` contrasts, signed lineage
root gaps, and four normalized P4 projectors across two doubling lineages. It
also distinguishes the unnormalized H4 sign reversal from the positive ratio
for size-normalized P4 quantities. A migration needs separate typed quantity
families plus explicit lineage, sign, and normalization maps; a single channel
string is insufficient. Numerical full-curve, covariance, frozen-model, and
report contracts remain unchanged.

`score_p50_fullcurve_n290.py` is migration-required. It carries thermal-even
DeltaM coordinates, a mean slope, a signed lineage root gap, and four P4
diagnostics through independent N145/N290 streams. `FEATURE_ORDER`, ordered
Gaussian representations, `LINEAGE_SIGN`, and frozen ratios jointly define the
semantics. A migration needs distinct typed quantities plus explicit lineage
and independent-stream maps. It must preserve size-local jackknifes, numerical
covariance-rank handling, frozen prediction order, and provenance.

`score_p50_sprime_n290.py` is now a covered frozen kernel. Its typed
entrypoint validates the registered exact identity map for the cross/even,
angular-normalized `P4_S_prime` observable before replaying the kernel. The
semantic gate separately freezes N=290, the q2-before-Jordan order, both source
prediction identities, and the fact that this score reuses the P50 raw block.
The scalar variance calculation, prediction hashes, chronology, numerical
scores, and declared decision text are unchanged.

`score_intrinsic_quantile_center_n145_n290.py` is now a covered frozen kernel. Its typed entrypoint validates an exact identity map for the cross/matching scalar value before replaying the frozen N145-to-N290 coordinate score. The semantic gate freezes u={0.025,0.05}, feature and residual order, Q and width normalization powers, independent RNG domains, and the zero cross-size covariance contract. It does not create independent evidence, change delete-one reconstruction, refit the 2^{-3/4} target, or alter numerical results.

The Issue #43 full-curve base scorer and locked wrapper are now covered as one
operational unit. Their typed entrypoints validate exact N185-to-N265 identity
maps for the `cross`, matching-odd `DeltaM` and matching-even `DeltaS` raw
first-minus-second contrasts before scoring. The semantic gate freezes
`p_ref`, sector order, both kernel identities, prediction hash, source-error
correlation, independent target streams, and the exact production allocation.
The locked entrypoint still activates the original metadata and joint-moment
validators; neither entrypoint refits amplitudes or changes reconstruction,
covariance, numerical scoring, or the no-refit rule. Both sectors come from
the same runs and are not declared independent evidence. The separately typed
cross/either correction remains an erratum entry rather than a replacement for
this two-sector contract.

`score_issue43_secondary.py` is migration-required. It consumes the primary
Issue #43 score and constructs a frozen ordered ledger containing reused
`DeltaM/DeltaS`, an x17 radial competitor, a zero benchmark, an intentionally
not-scorable H4+H12 stage, and an optional P48 `P4_S_prime` stage. The stage
names, ordering, hashes, and numeric payloads act as the semantic contract.
A migration must type the quantities, model/stage maps, excluded/not-scorable
states, and evidence reuse. It must preserve the no-target-refit and raw-data
boundaries. Its direct call to the untyped P48 kernel is not covered by the
separate prospective typed wrapper.

The norm-4 scalar production scorer and thermal-jet scorer are both
migration-required and share one evidence boundary. The scalar path orders
`U`, `P4_D`, and root-gap coordinates across two lineages and six sizes, then
applies frozen q2/Jordan and secondary transforms. The jet path orders ranks
2--6 after Hermite--Krawtchouk projection and width normalization, then applies
one cocycle multiplier to both lineages. A migration must type every quantity,
mode, normalization, lineage/size transform, and model order. It must also
record that scalar and jet outputs reuse the same histograms and cannot be
added as independent evidence. Existing covariance blocks, multipliers,
provenance, delete-one construction, and numerical scores remain fixed.

`score_norm5_thermal_jet.py` is migration-required. It recomputes the intrinsic
center, ranks 2--6 Krawtchouk thermal jet, canonical rank-gap width, and three
point lineage residuals before applying width-collapse and frozen q2/Jordan
cocycles. Counter-identical runs use synchronized delete-one covariance groups,
while disjoint groups contribute independent blocks. A migration must type the
modes, normalization, lineage maps, model order, and covariance-group relation.
It must also preserve that width and cocycle diagnostics reuse the same raw
curves and are not additive evidence. Prediction chronology, multipliers,
numerical-rank cutoff, jackknife construction, and scores remain unchanged.

The two v14 scalar projectors are migration-required. The fixed-p path selects
a free-form matching-function channel, applies cos4 weights to cancel H4, and
reports a `N^(25/8)`-scaled scalar at explicit `p_ref`. The full-curve path
instead solves two orientation roots, forms H4-null scalar and H4 coordinates,
then scores fixed `beta=7/2` GLS and two ordered norm-2 lineages. A migration
must type sector, channel, orientation weights, fixed-p versus root quantities,
units/scalings, and lineage maps. It must preserve covariance reconstruction,
root solving, synchronized delete-one batches, the conditional parity boundary,
and all limitations and numerical results.

`score_p159_pell_hex_filter.py` is migration-required despite serializing an
`observable_descriptor` dictionary in its output. That dictionary is free-form
provenance, not a registered `ObservableDescriptor` checked through
`map_observable`. The score orders three primitive rank-1 homology lines in a
transported positive-rho basis, transforms their continuum-subtracted
probabilities into C/Q/S character coordinates, compares two Pell designs
through E4 phase transport, and explicitly reuses PR #222 evidence. A migration
must type the primitive-sector channel, line-basis transport, character
coordinates, fixed-p and amplitude normalizations, cross-design phase map, and
non-independent evidence relation. It must preserve the exact oracle, frozen
continuum baselines, covariance transform, post-reveal boundary, and all
existing gate conclusions.

`score_p48_new_geometry_channels.py` is now a covered frozen kernel, separate
from the typed wrappers for the prospective `P4_S_prime` scorer. Its typed
entrypoint validates four exact cross-size identities: angular-normalized
matching-even `P4_S/P4_S_prime` use the frozen second-minus-first order, while
matching-odd `P4_D/P4_D_prime` use first-minus-second. The semantic gate also
freezes value-versus-first-derivative coordinates, the four exact N powers,
scaled keys, N=65/85/130 source, independent N=185/265 targets, canonical input
hashes, and shared-source covariance contract. It preserves every numerical
score, the no-target-refit rule, and the distinction from fixed-coordinate
P31/P43 `either/even` DeltaS. The four summaries share source and target blocks
and are not promoted to four independent evidence rows.

`threshold_score_modes.py` is migration-required rather than a generic
Krawtchouk helper. Although it exposes reusable basis functions, its operational
entrypoint reconstructs the intrinsic center inside synchronized delete-one
replicates, pairs two ordered orientations, forms matching S/D sectors, divides
by `DeltaCos4`, emits four P4 value/derivative views, and assigns a parity-tower
scaling convention. A migration must type the threshold-rank channel, intrinsic
center, orientation pair, S/D sector, Krawtchouk order, angular normalization,
derivative coordinate, and parity-tower units. It must preserve the exact mode
0/mode 1 identities, aligned jackknife covariance, and the guard that those
views reuse existing evidence rather than creating independent blocks.

The axis-pair annihilator score path is migration-required as one operational
unit. The base scorer reconstructs the cross matching function from
`K_minus/K_plus` rank histograms, solves implicit ordinary and annihilator
roots, and compares an `L^(13/4)`-rescaled adjacent-size combination. The
stable entrypoint corrects batch-reader ordering and delegates to that base
without adding a semantic gate. A future migration must type both the implicit
root-location quantity and the ordered adjacent-size relationship, while
preserving the stable reader contract; this audit does not change its frozen
numerics or reinterpret existing score artifacts.
