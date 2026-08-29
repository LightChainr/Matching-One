# Issue #156 norm-2 H4/H8 production job

The first primitive-sector pilot resolves one real nontrivial C3 mode but
cannot distinguish its conjugate H4 and H8 interpretations. The active
Gaussian multiplier `1+i` rotates the embedding by `pi/4` and doubles its
area. In the existing `Z=C+iQ` convention its exact action is `-I` for H4
and `+I` for H8. With the live `C~N^-1` radial transfer this gives the two
fixed lineage nulls

```text
H4: 2*C_child+C_parent = 0
H8: 2*C_child-C_parent = 0
```

The executable job contains both Pell sides:

```text
N30  [[6,3],[0,5]]  -> N60  [[6,-2],[6,8]]
N56  [[8,4],[0,7]]  -> N112 [[8,-3],[8,11]]
```

Each child is exactly `[[1,-1],[1,1]] P_parent`. All four designs receive
different deterministic seed blocks. Parent and child are **not** common
random numbers; their cross-covariance is zero by independent-stream design.
Equal batch labels across designs do not create pairing.

The batch CSV stores all category counts needed to reconstruct sector
probabilities. The JSON stores every design's complete `C/Q/S` covariance,
the block-diagonal six-coordinate covariance for each lineage, and the joint
covariance of the H4/H8 nulls. `Q` remains a reflection/convention null and
`S` remains a scalar diagnostic.

Production command (not run in this commit):

```bash
python3 scripts/square_bond_primitive_norm2.py \
  --samples-per-design 2000000 \
  --batches 200 \
  --seed 20260829 \
  --workers 16 \
  --dps 80 \
  --output-prefix results/server-20260829/P156-norm2-h4-h8/result
```

Only parsing and tiny-sample execution belong in CI:

```bash
python3 -m unittest discover -s tests -p 'test_square_bond_primitive_norm2.py'
```
