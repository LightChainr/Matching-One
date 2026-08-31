# P334/P403/P429 — exact two-white-cycle trigger bipartition

- **Lifecycle:** branch-only / exact topology theorem with explicit rank-duality
  input / deterministic replay / no new Monte Carlo / no PR or issue comment.
- **Changes:** minimal two-site trigger graphs on a black-NN rank-one torus
  state are bipartite. A fixed deletion-compatible white site/face spine lies
  in an annulus; a proved integer face-distance packing lemma gives two
  essential white cycles disjoint on every singleton-safe site. Every trigger
  pair hits one site of each cycle. Trigger triangle count is therefore exactly
  zero, rather than merely zero in the 22 selected graphs.
- **Input boundary:** square NN/matching ambient-rank duality, not ordinary
  graph cycle rank. The honest-cell model is unconditional within the existing
  digital-Alexander theorem; arbitrary integer periods use the repository's
  unrestricted-rank extension `73d4960`. No integral-saturation premise is
  needed. The annular packing step is proved in the note, not assumed.
- **Observer/sector/source/geometry:** unordered minimal pair additions on
  occupied rank-one states; white matching essential-cycle deletion dual;
  N325/N425 first/second Gaussian configurations, plus bounded exact quotients.
- **Existing source:** `1b5a9dea07e1c62f69798fddbf4899ff986c0b72`, selected
  `trigger_graph_raw/` files under `results/local-20260831/P334-cooperative-closure/`.
  Same `p334-cooperative-N325-20260831` / `p334-cooperative-N425-20260831`
  dependency groups. The 22 graphs are not a fresh replication block.
- **Artifacts:** `archived_two_cycles.json` saves actual white site/face walks,
  deck addresses and each trigger edge's crossing of the cycle sides; all 22
  pass. `tiny_census.json` covers 36 quotients, 3,388 rank-one states and 1,804
  nonempty pair graphs, with no counterexample. Two focused walk-certificate
  tests pass; no old Ferrers, triple-census or full-repository tests rerun.
- **Does not prove:** Ferrers/chain structure (already false), sufficiency of
  the labelled pair graph for triple survival, continuum field identity, an
  exponent, intrinsic temporal memory or a population causal allocation.
- **Next useful consequence:** model the observed side-capacity and overlap
  using annular two-sided separators; the exact three-step identity retains
  genuine triple nonfaces but no trigger-triangle term. Further raw sampling
  is not requested by this result.
- **Execution:** managed local `/Users/lc/python-envs/research-py311/bin/python`,
  integer circulation for certificate construction; no server connection.
