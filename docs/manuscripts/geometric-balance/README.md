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
