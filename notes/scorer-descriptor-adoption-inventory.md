# Scorer descriptor-adoption inventory

In this checked snapshot, the exact `scripts/*score*.py` corpus contains 51
files. Python-AST inspection finds twenty
direct typed entrypoints importing both `ObservableDescriptor` and
`map_observable` from `wrapping_channels`. Nineteen frozen kernels are covered by
those entrypoints through explicit semantic-manifest wrapper relationships.
One path, `kappa3_half_score.py`, is classified as a generic utility for which
a descriptor is not applicable. Eleven paths are confirmed channel-bearing and
require typed migrations:
`score_issue43_secondary.py`,
`score_norm4_production.py`,
`score_norm4_thermal_jet.py`, `score_norm5_thermal_jet.py`,
`score_p159_pell_hex_filter.py`,
`score_p231_vacuum_kdv_sector.py`,
`score_p49_fullcurve_doubling.py`,
`score_p50_fullcurve_n290.py`,
`score_prequential_evidence.py`,
`threshold_score_modes.py`,
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
open because inventory classification is not the same as implementing the 11
required typed migrations.

`score_axis_pair_annihilator.py` and its stable-reader entrypoint are now
covered frozen kernels. Their two typed entrypoints share one semantic gate
and validate the exact identity for a `cross` matching-odd raw value before
delegating. The gate separately freezes the adjacent upper/lower roles,
zero-based aligned batches, `p_ref`, the `L^(13/4)` annihilator, q/w ordering,
and source-fit/held-out-no-refit boundary. The implicit root and candidate
operator interpretation remain response/model coordinates rather than newly
asserted topology identities; both historical numerical kernels are unchanged.

The `kappa3_half_score.py` exception is narrow. It defines only exact Bernoulli
likelihood-score polynomials at `p=1/2` and aggregates an opaque caller-supplied
`D(C)`; it intentionally does not define a wrapping, homology, or other
topology observable. Geometry-specific callers still must type their own
`D(C)` before comparison.

`score_angular_root_amplitude.py` is now a covered frozen kernel. Its typed
entrypoint applies the registered raw-to-angular-normalized map separately for
the signed N=65 and N=85 orientation designs before replaying the scorer. The
gate freezes first-minus-second order, signed `DeltaCos4`, the implicit-root
response, the additional `-N^2` transform, aligned common-stream batches, full
cross-size covariance, and the frozen prediction with zero target refits. This
types the topology contrast without reinterpreting root location or finite-size
normalization as a new exact topology identity; every completed P45 number is
unchanged.

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

`score_rank_gap_boundary_targets.py` is now a covered frozen kernel. Its typed entrypoint validates the exact cross/primal identity used by the same paired rank-gap observable across source and target sizes before replaying the frozen scorer. The semantic gate freezes `G=K_plus-K_minus`, rank units, equal first/second orientation pooling, source and target order, the unfitted `5/8` exponent, and `(A,B)` parameter order. It preserves pre-reveal chronology, disjoint counters, covariance, predictions, and all numerical results; rank units are not reinterpreted as probability values.

`score_threshold_rank_root_doubling.py` is now a covered frozen kernel. Its typed entrypoint uses the registered raw orientation-order reversal to obtain the exact stored-child transform `(-1,0)` before evaluating the two `1+i` lineages. The semantic gate freezes sizes `[65,85,130,170]`, lineage order `(65,130)` then `(85,170)`, the threshold-rank root-gap quantity, fixed `-1/4` ratio, and full-versus-diagonal covariance order. The full covariance remains primary, the diagonal score remains diagnostic, and no numerical output changes.

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

`score_v14_fixedp_scalar_projector.py` is now a covered frozen kernel. Its typed entrypoint locks the basis-dependent `direction_1/matching` value, fixed-p coordinate, `p_ref=0.592746050790`, H4-null orientation weights, covariance reconstruction, and `N^(25/8)` normalization before delegating row projections. It remains a retrospective discovery/power diagnostic and neither proves the V_<1,4> assignment nor removes common thermal displacement. The separate full-curve `score_v14_scalar_root_projector.py` remains migration-required because it solves implicit roots and scores two lineages.

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

The axis-pair annihilator base and stable-reader paths are now covered frozen
kernels, not migration-required adoption paths. Their shared semantic gate
types the `cross` matching-odd topology value, adjacent upper/lower roles,
orientation order, angular normalization and stable batch-reader contract.
The implicit ordinary and annihilator roots, the `L^(13/4)` adjacent-size
relationship and the candidate `q -> w` interpretation remain bespoke
response/model metadata. They are a second-layer contract gap, not evidence
that these two paths remain untyped in the current adoption inventory; the
wrappers change neither frozen numerics nor existing score interpretation.
