# Completed work and corrections for the existing issue channels

2026-09-12. Prepared handoff; not yet posted to GitHub. No duplicate issue,
new production, hardware allocation, merge, or STATUS promotion is requested.

## #711 / #717: dictionary correction, not another bibliography search

Withdraw #717 §2.3's categorical distinction between critical-polynomial 2D/0D
and the same finite site ensemble's rank2/rank0. Read Mertens--Ziff 2016,
arXiv:1603.07289v2 introduction and §II Eqs. (12),(20),(21) with the following
paragraph. The finite event identity is already explicitly connected there.
Use `notes/finite-critical-polynomial-dictionary-20260912.md` for the exact
periodic-lift and site-weight dictionary, including the rank-one spiral trap.
Preserve the distinct correct boundary: no all-width B5/B15/B16 matrix
intertwiner or closure-weight theorem follows. Also do not attribute uniqueness
to the intermediate value theorem alone.

The remaining literature task should compare actual local site weights,
embedding, topological closure and scalar spectral decompositions, not just
whether the papers use identical names. Missing small-block formulas in a
bounded reading are not proof of novelty.

## #636: source compatibility and rare-sector conditioning delivered

On the existing PR708 finite automaton, all seven D4 types of column-probability
groupings are now classified exactly. Common strong lumpings have 94/303/179/
262/509/303/509 states, respectively; every one equals the corresponding
colour-preserving D4 orbit partition blockwise. Two independently addressed
adjacent columns already require all 509 deterministic rank classes within
this lumping class. No further width-four state enumeration is needed.

A genuine square-site two-row dipole source has zero linear response but
nonzero mixed response. The 4x4 coefficients are 327/1024 and 633/2048 at
p=1/2; the finite root's mixed shifts have strict rational sign certificates.
Four-sign extraction of M is exact at finite amplitudes, by separate degree<=2
multiaffinity. This does not inherit P398's continuous-time response formula.

An exact integer backward sampler conditions on any final rank. The complete
nonuniform 4x3 conditional law is verified configurationwise; the 4x128 rare
case has exact normalizers and fixed-seed algorithm controls. Do not commission
billions of direct snapshots merely to see rare sectors that this oracle can
condition on. Conversely do not claim the fixed-width oracle scales cheaply
to large circumference.

Files are additive on PR708, no dependence on unmerged #710 or #716.
All-width structural closure weights and width-uniform spectral estimates
remain different, uncompleted questions. This deliverable does not authorize
a next-width scan or a new transfer engine.

## #275 / #337: what this does and does not supply

The full multivariate rank event polynomial fixes occupation-source derivatives
with the correct physical normalizer. It does not fix a generic-q continuation,
a homology-line-marked source, or the candidate-specific restricted-trace and
moving-root original-U map. Do not use the extra site sources as a post-hoc
rescue of a frozen candidate or claim continuum identification from a nonzero
mixed finite response. The source-safe exact oracle can verify a specified
future forward map once that map is actually supplied.

## Verification boundary

The completed archive includes proofs, executable standard-library-only code,
full labels/actions, rational source/root certificates, exact conditioned path
probabilities, and a minimal-workspace patch application/reproduction check.
The repository's complete CI suite was NOT executed for these new files.
Past successful CI on #708/#710 is not a substitute.
