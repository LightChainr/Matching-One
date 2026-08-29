# Second norm-2 generation for the primitive C3 character

PR #238 selects the negative rank-4 phase on the unseen N60/N112 children,
but the two ratios `-0.4166` and `-0.4382` remain short of the asymptotic H4
target `-1/2`.  This prospective follow-up applies the same `1+i` cover once
more, without refitting the exponent:

```text
N60  [[6,-2],[6,8]]  -> N120 [[0,-10],[12,6]]
N112 [[8,-3],[8,11]] -> N224 [[0,-14],[16,8]]
```

The primary question is whether the two independent residuals
`C_child + C_parent/2` remain compatible with zero.  The positive ratios
`+1/2`, `+1/4`, and `+1/8` retain their frozen secondary order.  Comparing
the first- and second-generation distances from `-1/2` is descriptive only;
it cannot change the fixed-model score.

Production uses five million fresh configurations for each child and reads
the generation-1 parent values/covariance from the committed PR #238 result.
No cross-generation CRN covariance is claimed.

```bash
python3 scripts/square_bond_primitive_norm2_generation2.py \
  --parent-result results/server-20260829/P156-norm2-h4-h8/result.json \
  --samples-per-design 5000000 --batches 250 --seed 202608292 \
  --workers 16 --dps 80 \
  --output-prefix results/server-20260829/P156-norm2-generation2/result
```

## Revealed result

The frozen Huawei run completed at commit
`3eb423da494064f27f4474553abaaa2bf463a73d`.  Both second-generation children
again carry the alternating H4 sign:

| lineage | parent C | child C | child/parent | H4 residual z |
|---|---:|---:|---:|---:|
| N60 -> N120 | `-0.00359392(291)` | `+0.00204393(172)` | `-0.568717` | `+1.096` |
| N112 -> N224 | `-0.00200423(302)` | `+0.000817119(211)` | `-0.407697` | `-0.712` |

The fixed `-1/2` model gives `chi2=1.70765/2 df`, `p=0.42578`.  The ratios do
not move monotonically toward `-1/2` lineage by lineage, so the descriptive
"closer on both" diagnostic fails; nevertheless their opposite residuals are
jointly consistent with the parameter-free H4 law.  All frozen positive-phase
models fail: `+1/2` gives `chi2=339.467`, `+1/4` gives `282.808`, and `+1/8`
gives `225.996`, each with two degrees of freedom.

The new child C effects are independently resolved (`z=11.89` and `3.86`).
Their Q reflection controls remain null (`z=-0.782,+0.835`).  Combining this
with generation one gives two successive exact sign alternations on each
lineage, while keeping the scale-ratio conclusion conditional on visible
finite-size drift.

The four transition residuals are not independent because each generation-1
child is reused as a generation-2 parent.  Propagating the exact shared-middle
covariance, `Cov(r1,r2)=Var(C_middle)/2`, gives the combined two-lineage,
two-generation H4 score `chi2=7.33031/4 df`, `p=0.11943`.  Thus the complete
six-geometry chain remains compatible with one fixed `-1/2` transfer; simply
adding the four individual z-squares would double-count the middle points.
