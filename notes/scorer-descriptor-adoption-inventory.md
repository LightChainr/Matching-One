# Scorer descriptor-adoption inventory

At main commit `11ba86892fabd3349e322d198f23ef24d8cbe828`, the exact
`scripts/*score*.py` corpus contains 35 files. Python-AST inspection finds four
direct typed entrypoints importing both `ObservableDescriptor` and
`map_observable` from `wrapping_channels`. Three frozen kernels are covered by
those entrypoints through explicit semantic-manifest wrapper relationships.
One path, `kappa3_half_score.py`, is classified as a generic utility for which
a descriptor is not applicable. Ten paths are confirmed channel-bearing and
require typed migrations: `score_angular_root_amplitude.py`,
`score_axis_pair_annihilator.py`, its stable entrypoint,
`score_c4_self_matching_n26.py`, `score_c4_tangent_orthogonal_holdout.py`,
`score_issue50_n290.py`,
`score_matching_odd_synthesis.py`, `score_prequential_evidence.py`,
`score_rank_gap_boundary_targets.py`, and
`score_threshold_rank_root_doubling.py`. The remaining 17 files are outside a
registered typed path.

This closes an inventory gap: the repository now has a deterministic list,
including Git blob identities, rather than an informal claim that all scorers
have or have not adopted descriptors. The audit fails closed if a direct typed
import is added or removed without updating the manifest, or if a declared
wrapped kernel disappears.

`outside_registered_typed_path` is deliberately a triage label. Some of those
17 scripts may not compare channel-bearing quantities, some may be historical,
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

`score_prequential_evidence.py` is also migration-required. Its scored-block
gate currently compares free-form `channel.source` and `channel.target` strings
and treats any truthy `exact_map` field as sufficient. A typed migration must
parse source/target descriptors and verify the registered affine transform,
while preserving historical ledger rows and chronology.

`score_c4_self_matching_n26.py` is migration-required even though its committed
N=26 control finds all five wrapping-channel Bernstein vectors identical. It
selects the scored observable from the prediction's free-form
`geometry.wrapping_channel` string, then compares that exact law with frozen
hypotheses. Numerical equality on one finite self-matching geometry does not
establish semantic interchangeability. A typed migration must bind the frozen
prediction and result to the exact channel descriptor while preserving the
pre-target scoring order and stop-without-generalized-fit rule.

`score_c4_tangent_orthogonal_holdout.py` is migration-required. It selects the
`cross` response through a free-form channel string, freezes `lambda/t` on
N=130, and applies that projection to aligned N=170 batches before scoring
orthogonal and thermal residuals. A typed migration must bind both sizes to the
same exact channel semantics and type the response-coordinate pair, without
changing batch alignment, jackknife construction, or the frozen interpretation
rule.

`score_matching_odd_synthesis.py` is migration-required. It selects two frozen
ledger blocks by requiring free-form `channel.source` and `channel.target`
strings to equal `matching_odd`, then combines their already-primary scores.
A later migration must parse both typed descriptors and verify their registered
identity maps. That semantic gate must not change the frozen block selection,
distinct raw-data-group requirement, block-diagonal synthesis, or the rule that
the derived output cannot become a new primary evidence row.

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

`score_issue50_n290.py` is migration-required. It computes the
`either/matching_function` contrast between ordered Gaussian lineages
`(13,11)` and `(17,1)` at fixed p, with channel, sector, lineage order, and
subtraction sign encoded as strings and constants. A typed migration must bind
the exact observable and ordered lineage map while preserving the prospective
target, run provenance, counter-range validation, and sampling-error score.

The axis-pair annihilator score path is migration-required as one operational
unit. The base scorer reconstructs the cross matching function from
`K_minus/K_plus` rank histograms, solves implicit ordinary and annihilator
roots, and compares an `L^(13/4)`-rescaled adjacent-size combination. The
stable entrypoint corrects batch-reader ordering and delegates to that base
without adding a semantic gate. A future migration must type both the implicit
root-location quantity and the ordered adjacent-size relationship, while
preserving the stable reader contract; this audit does not change its frozen
numerics or reinterpret existing score artifacts.
