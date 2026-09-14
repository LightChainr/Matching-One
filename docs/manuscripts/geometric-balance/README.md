# Geometric balance: one manuscript, one new lemma

Read [manuscript.md](manuscript.md). It combines the root theorem of #735 with
#736's full-law question, and supplies the missing arbitrary-direction
staircase corridor. The resulting all-period full-law criterion is
`log N / ell -> 0`; balance-root consistency only needs `ell -> infinity`.

The new proof uses axis-aligned NN rectangles and disjoint translates in the
finite quotient group. It does not rotate the physical interaction, need an
ambient-primitive shortest vector, or import a fixed-width continuum limit.
No unrelated source/Jordan calculation is a dependency.

From the repository root:

```sh
python -m unittest discover -s tests -p 'test_oblique_winding_corridor.py' -v
python scripts/oblique_winding_corridor.py --output /tmp/oblique-corridor-new.json
```

The script refuses to overwrite an existing result. The committed result is
`results/research-control-20260913/oblique-corridor-controls.json`.
Python standard library only. Five local tests and 135,168 tiny configurations
were executed; full repository CI was not run. No Monte Carlo was performed.

The manuscript states the external RSW, site-sharpness and matching inputs,
contains the entire root and concentration arguments, and gives a bounded
closest-source comparison. It is an author-supplied proof; no independent
publication acceptance or originality certification is claimed. Existing
proofs, data, frozen designs and research branches are retained unchanged.

## 2026-09-14 integrated structural continuation

The stacked continuation on PR #771 collects consequences that reduce several
previously separate analysis directions without changing the parent #739
acceptance status.

### Core structural notes

- [`structural-consequences-20260914.md`](structural-consequences-20260914.md): persistent 4/8 birth reflection; dual-even/odd birth coordinates; exact same-parameter `(rank,K)` count reduction; alternating black/white barrier Palm identities; marked-Poisson transport; complementary Palm score constraints; convex loop/branch frontier; the `D^{-1}=partial_yy tau` consistency relation; and deterministic directional separation in exponential elongation.
- [`rare-charge-balance-20260914.md`](rare-charge-balance-20260914.md): exact finite factorization of the matching-root slope into rare topological-charge susceptibility and conditional endpoint-odds slope, making `balance without concentration` algebraically explicit.
- [`projective-homology-gas-20260914.md`](projective-homology-gas-20260914.md): exact rank-one state `(slope,K)` and the resulting projective hard-core-gas language for multi-direction competition.

### Alternating black/white geometry

- [`poisson-tessellation-consequence-20260914.md`](poisson-tessellation-consequence-20260914.md): explicit total-variation/Markov-kernel derivation of the Exp component-Palm gap, Gamma(2,1) stationary-location gap, uniform relative position, and higher Gamma spacings.  The coarse two-colour limit is a Poisson interval tessellation, not two independent Poisson clouds.
- [`supercritical-white-slab-bulk-20260914.md`](supercritical-white-slab-bulk-20260914.md): exact infinite-cluster boundary-density identity `beta(q)=(1-q)theta(q)/q` and a concrete bulk LLN/CLT target for the huge complementary white component.

### Prefactor and directional mass

- [`matrix-sewing-unit-residue-20260914.md`](matrix-sewing-unit-residue-20260914.md): finite-state matrix cyclic-sewing theorem.  A simple Perron band with finite Markov memory still gives exactly `exp(-kappa w)/sqrt(2 pi D w)` with unit logarithmic residue; finite local memory alone cannot explain an anomalous power or amplitude.
- [`matching-enhancement-mass-gap-20260914.md`](matching-enhancement-mass-gap-20260914.md): author-level two-terminal enhancement proof, using Grimmett--Li plus the corrected square-lattice rerouting theorem of Balister--Bollobas--Riordan, giving `kappa_8(p)<kappa_4(p)` for `p<p_c(G8)` and hence the strict centre inequality `a(d)+b(d)>1`.
- [`research-frontier-20260914.md`](research-frontier-20260914.md): the homological-free-energy variational picture, a small-`p` actual-SITE transfer/local-CLT programme for the complete-component `w^{-1/2}` prefactor, and the recommended order for the remaining common-window work.

Finite controls for the persistent reflection are in
`scripts/persistent_alexander_birth_reflection.py`, with committed L=3 and L=4
outputs in `results/geometric-consistency/`. The L=3 run exhausts all 512
configurations and all 362,880 strict site orders; the L=4 run exhausts all
65,536 configurations and checks 20,000 fixed-seed site orders. A separate
zero-cost calculation of the constrained-Poisson/topology constants is in
`scripts/poisson_topology_constants.py`.

These notes deliberately separate deterministic consequences, author-level
proofs requiring independent review, deductions that use existing #739
probability inputs, and conjectural proof programmes. They do not identify a
continuum field, claim a full all-subcritical SITE sewing theorem, or turn the
near-critical crossover into an accepted square-site theorem.
