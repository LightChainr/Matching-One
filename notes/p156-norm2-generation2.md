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
