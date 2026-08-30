# Issue #155: large-N gate for the local odd pivotal readout

## Frozen pilot

The N=10 oracle fixes the two observable rows as

```text
global = [cross(black) - cross(complement)] / 2
local  = [H4_pivotal(black) - H4_pivotal(complement)] / 2.
```

The production stream uses the C4 checkerboard self-matching quotients
`N=130, (a,b)=(11,3)` and `N=170, (13,1)`, local radius `R=3`, and the exact
Bernoulli scores `(S_t,S_lambda)` at `p=1/2`. Both sizes use counters
`0..199999`, seed `15520260829`, and 100 aligned batches. Site bits are prefix
coupled across sizes, retaining the declared common-random-number structure.

The pre-existing gate is applied separately at each size:

```text
abs(determinant z) >= 3 and condition number <= 50.
```

Every delete-one-batch replicate reconstructs the full response matrix,
determinant, both singular values, and condition number. A generalized
eigensystem is permitted only when both size gates pass.

## Result

Rows are `(global,local)` and columns are `(t,lambda)`:

```text
R130 = [[4.70673, 3.18759],
        [0.11404, 0.07316]]

det = -0.0191684 +/- 0.0379338, z = -0.5053
singular values = (5.68616, 0.00337106)
condition number = 1686.75

R170 = [[5.20708, 3.52986],
        [0.10393, 0.08295]]

det = 0.0650689 +/- 0.0418964, z = 1.5531
singular values = (6.29215, 0.0103413)
condition number = 608.45
```

Both matrices fail both parts of the frozen rank gate. Consequently the
analysis reports only one resolved tangent dimension and suppresses
`R170 R130^-1`, generalized eigenvalues, and effective `y` values.

This is not a null local-observable result. The local entries are individually
resolved: at N130 they are `0.11404 +/- 0.00579` and
`0.07316 +/- 0.00632`; at N170 they are `0.10393 +/- 0.00626` and
`0.08295 +/- 0.00599`. What remains unresolved is their linear independence
from the much larger global thermal-like row. The determinant changes sign
between sizes within uncertainty, and the second singular values remain too
small to support a second RG eigenvalue.

The empirical score Fisher matrices remain close to the exact `4N I_2`
controls. Thus this bounded run answers the intended decision: the exact
microscopic rank-two oracle does not survive as a statistically resolved
two-dimensional large-N response at this sample size, and no free exponent
diagnostic is licensed.

## Reproduction

```bash
g++ -O3 -std=c++17 -fopenmp src/c4_local_odd_pivotal_mc.cpp \
  -o c4_local_odd_pivotal_mc

./c4_local_odd_pivotal_mc \
  --samples 200000 --batches 100 --threads 10 --radius 3 \
  --seed 15520260829 --replica-offset 0 --git-commit 9afd7de \
  --output-prefix results/local-20260829/P155-c4-local-odd-pivotal-largeN/raw/n130_n170_200k

python3 scripts/analyze_c4_local_odd_pivotal_mc.py \
  --batches results/local-20260829/P155-c4-local-odd-pivotal-largeN/raw/n130_n170_200k.batches.csv \
  --metadata results/local-20260829/P155-c4-local-odd-pivotal-largeN/raw/n130_n170_200k.metadata.json \
  --output results/local-20260829/P155-c4-local-odd-pivotal-largeN/analysis.json
```
