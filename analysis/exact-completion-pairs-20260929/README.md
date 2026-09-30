# Exact conditional averaging of the existing safe-insertion block

This is a post-readout estimator refinement of the fixed
[safe-insertion target](../safe-insertion-production-20260929/contract.json),
not a new acquisition, changed source, or independent confirmation.
It replays exactly the same 14 x 10,000 permutations and replaces each
four-probe average by its exact conditional expectation.

## Geometry and target

The [projected component-pair proof](../../notes/projected-component-pair-criterion-20260929.md)
gives a necessary and sufficient criterion for two individually safe vacancies
to complete ambient rank two together. Integer transverse potentials detect
inconsistent shared-component or physical-edge constraints. Component-pair
buckets generate only true synergy edges, then deduplicate them. No quadratic
all-vacancy-pair scan is used in production. The four small physical controls
do use direct insertion for their reference pair sets.

For one rank-one prefix, let m=N-b, c=nu_b, s=m-c and e be the number of
unordered safe synergy pairs. For s>0:

    E[nu(A+V)-nu(A) | A, V uniform safe] = 2e/s.
    P(survive two more insertions | A) = [s(s-1)-2e]/[m(m-1)].

Conditioning on the full original prefix sample also fixes J1, directions,
cell membership and overlap weights. Thus replacing the actual probe mean
by 2e/s preserves the finite-sample conditional target and removes probe
randomness; it does not remove configuration variation. This argument does
not assume the current pair graph closes recursively.

The proof's generic expected bound O(N*z^2+z^2*e) uses merged contact lists
in its adjacency scan. The implementation instead uses bounded nested lists,
O(N*z^3+z^2*e); z is fixed at four (square) or six (controls), giving the same
expected O(N+e) fixed-lattice bound. Hash costs are expected, not worst-case;
large e can still make explicit enumeration expensive.

## Computation and inference

`engine.cpp` reuses the unchanged lifted union-find, bounded integer RNG and
full Fisher-Yates permutation stream. It builds transverse occupied-component
potentials separately, computes c and e, and saves one row per old prefix.
`run.py` takes all seeds and sample counts from the original manifest. It
records source hashes, compiler, commands, timing and the original dependency
groups. No new random prefix is generated; the first-128 timing sample occurs
again in batch zero and is not extra data.

`analyze.py` requires row-for-row equality of `(J1,rank_b,dx_b,dy_b,nu_b)` across
all 140,000 old/replayed rows. It retains original early/late definitions,
exact (D,c) support, weights n_early*n_late/(n_early+n_late), and the original
direction-only secondary target. Every whole-batch deletion recomputes
support and weights. Four actual old probes are averaged before inference.

Eight jointly scored coordinates include the exact and old probe contrasts,
their paired difference, the exact edge/two-step contrasts, cohort means,
and the unchanged secondary. Derived consumer factors are applied inside
cells; the two-step contrast is not another independent measurement.
The JSON contains integer sufficient sums, all 14 deletion vectors and full
8x8 covariance. The measured SE ratio is descriptive, not a guarantee that
every empirical jackknife SE must shrink. No cutoff, feature, sign or size
search is included.

## Reproduction

Owner routing preference updated 2026-09-30: medium/large replays go to
Huawei first; local work is for small calculations, scoring and editing.
This run began locally, was paused by the owner, then retained batches0..7
and transferred only the missing batches8..13 to TV2N0X. `resume.py` records
the first local restart; `finish_on_cloud.py` is this specific six-batch
handoff, not a generic scheduler. No stopped partial batch is scored.

To rescore the saved data without overwriting the report:

```sh
python3 analysis/exact-completion-pairs-20260929/analyze.py \
  --output-dir analysis/exact-completion-pairs-20260929/reproduction
```

To reproduce geometry from seeds, use an isolated checkout of the replay
source commit named in `data/run.json` (before replay outputs existed) on
the selected compute host; use `--compiler g++` on the Huawei image:

```sh
python3 analysis/exact-completion-pairs-20260929/run.py \
  --source-commit 8b5eee0d47b6c220c2ba46a5bae08d2d19b392a4 --workers 4
```

The runner deliberately does not overwrite an existing replay directory.
This is output preservation, not a research permission gate. The original
production and its published primary readout remain unchanged.
