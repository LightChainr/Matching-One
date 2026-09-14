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

Two follow-up notes collect consequences that reduce several previously
separate analysis directions without changing the main theorem's acceptance
status:

- [`structural-consequences-20260914.md`](structural-consequences-20260914.md): persistent 4/8 birth reflection; dual-even/odd birth coordinates; exact same-parameter `(rank,K)` count reduction; alternating black/white barrier Palm identities; marked-Poisson transport; complementary Palm score constraints; convex loop/branch frontier; the `D^{-1}=partial_yy tau` consistency relation; and deterministic directional separation in exponential elongation.
- [`research-frontier-20260914.md`](research-frontier-20260914.md): the homological-free-energy variational picture, a small-`p` actual-SITE transfer/local-CLT programme for the closed-component `w^{-1/2}` prefactor, and the recommended order for the remaining common-window work.

Finite controls for the persistent reflection are in
`scripts/persistent_alexander_birth_reflection.py`, with committed L=3 and L=4
outputs in `results/geometric-consistency/`.  The L=3 run exhausts all 512
configurations and all 362,880 strict site orders; the L=4 run exhausts all
65,536 configurations and checks 20,000 fixed-seed site orders.  A separate
zero-cost calculation of the constrained-Poisson/topology constants is in
`scripts/poisson_topology_constants.py`.

These notes deliberately separate deterministic consequences, deductions that
also use the author-level probability results already on this branch, and
conjectural proof programmes.  They do not identify a continuum field, certify
an OZ sewing theorem, or turn the near-critical crossover into an accepted
square-site theorem.
