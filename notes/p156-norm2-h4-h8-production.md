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

Production command:

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

## Huawei production result

The frozen job was run on `DevEnvC_XPk2PZ` at commit
`9a2e62b8e59b7a3599989ab1e8eff4286bbf578d`, with 2,000,000 independent
configurations per design and 200 batches.  Both previously unseen children
have a negative and individually resolved nontrivial character:

| lineage | `C_parent` | `C_child` | child/parent |
|---|---:|---:|---:|
| N30 -> N60 | `+0.00862633(319)` | `-0.00359392(291)` | `-0.416623` |
| N56 -> N112 | `+0.00457402(271)` | `-0.00200423(302)` | `-0.438177` |

For the fixed rank-4 H4 ratio `-1/2`, the two residual z scores are `2.165`
and `0.854`; their independent two-lineage score is
`chi2=5.41712/2 df`, `p=0.06663`.  The fixed even-character/H8-alias ratio
`+1/2` gives `chi2=734.487/2 df`, `p=3.22e-160`.  The secondary positive
ratios frozen before reveal are also excluded: `+1/4` gives
`chi2=465.534/2 df`, and `+1/8` gives `chi2=324.046/2 df`.

Thus the norm-2 rotation selects the negative rank-4 character.  At current
precision it is compatible with the parameter-free H4 transfer, with visible
finite-size drift on the N30 lineage, and decisively rejects every frozen
positive-phase mechanism.  The child reflection quadratures remain null
(`z=-0.381,-0.069`), so the sign result is not carried by the convention-null
coordinate.

The covariance-bearing archive and environment record are under
`results/server-20260829/P156-norm2-h4-h8/`.
