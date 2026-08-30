# N170 exact angle-flip reveal

Status: fresh held-out production generated only after preregistration commit
`cf1bdf8`.  The H4-only vector, absolute predictions, pair contrast,
curvature/projective basis, sample count, seed and server were frozen first.

## Frozen question

N170 `(11+7i,13+i)` is the exact `1+i` child of N85 `(9+2i,7+6i)`.
Its reflection-even H4 covectors are

```text
c = (-4633/7225, +6887/7225),
```

the negatives of the N85 covectors.  Meanwhile the A charged source retains
the exact projective scalar `q_A^2=(u+H_F3)/2`.

The model trained only on N65/N85 froze

```text
K_H4(N170) = (+0.003935533602, -0.005850209350),
Delta_K_H4 = -0.009785742953,
H4 amplitude = -0.006137325767,
A-projective scalar = 0.
```

The primary discriminator was negative `Delta_K_A` versus scalar/common-mode
zero.  The exact two-coordinate decomposition was also frozen:

```text
A_H4 = (K_second-K_first)/(c_second-c_first),
A_scalar = 2(c_second*K_first-c_first*K_second)/(c_second-c_first).
```

Departure in `A_H4` is scale curvature; departure in `A_scalar` is a
charged/projective common mode.

## Production

Huawei `DevEnvC_HZsCM6` (`033945d8bf8b47a7acf475c595169e07`) was idle when
selected.  The exact `1714141` projective-birth engine ran 8M samples/shape in
80 aligned batches with seed `202608337170`, 16 threads and replica offset
zero.  Engine time was 17.99 seconds.  Self-test passed and stderr is empty.

The sparse archive retains `tau1`, primitive `ell`, `tau2`, direct-rank2 type
and all per-batch counts, so the complete covariance is reconstructible.

## Reveal

At the unchanged `p_ref`:

```text
first  (11+7i): W_A=0.337637047140,
                 Jminus_A=+0.009000051687,
                 K_A=+0.006434708591;

second (13+i):  W_A=0.337701145631,
                 Jminus_A=-0.015783052891,
                 K_A=-0.011282165119.
```

Therefore

```text
Delta_K_A(N170)=-0.017716873710 +/- 0.001555600655,
z versus scalar zero = -11.389.
```

The exact geometry sign flip is not marginal: it is the dominant resolved
direction.  But its magnitude is larger than frozen:

```text
pair residual to H4 target = -0.007931130757,
predictive SE              =  0.002613638683,
z                          = -3.035.
```

The full two-component frozen-vector predictive score is `10.252/2`.

## Where the residual lives

The preregistered basis gives

| coordinate | frozen | observed | uncertainty | residual score |
|---|---:|---:|---:|---:|
| H4 amplitude | -0.00613733 | -0.01111149 | predictive SE 0.00163920 | -3.035 |
| A-projective scalar | 0 | -0.00138098 | measurement SE 0.00176480 | -0.783 |

The linear-basis and original-vector predictive quadratics agree to
`7.1e-15`.  Hence the excess is not an artifact of choosing coordinates.

The scientific result is sharper than the previous N145 reading:

1. exact geometry rotation determines the sign;
2. the H4-only radial amplitude underpredicts N170 at 3.03 predictive SE;
3. the residual stays in the same H4 geometry direction;
4. the orthogonal charged/projective scalar remains a clean zero control.

So the N145 curvature-like central remainder was real enough to recur in an
exact same-lineage flip.  What is missing is a curvature law for the charged
H4 amplitude, not another harmonic vote or a new projective state.

As a secondary descriptive closure, the measured H4 amplitude changes from
`-0.0160509` at N85 to `-0.0111115` at N170, ratio `0.6923`; the fixed
`2^-13/8` radial ratio is `0.3242`.  This secondary ratio was not used to
select the preregistered model.

## Gates and boundary

- all production metadata match the `cf1bdf8` freeze;
- current continuity closes to `2.12e-14`;
- raw sufficient statistics and hashes are committed;
- N170 is independent of N65/N85/N145;
- no H4/H8 comparison and no exponent fit was performed.

The result identifies geometry sign and residual direction.  It does not yet
identify an asymptotic curvature form or continuum exponent.

## Reproduction

```bash
python3 scripts/score_n170_angle_flip.py \
  --preregistration analysis/p337_n170_angle_flip_preregistration.json \
  --births results/server-20260830/P337-N170-angle-flip/raw/n170_8m.births.csv \
  --metadata results/server-20260830/P337-N170-angle-flip/raw/n170_8m.metadata.json \
  --json results/server-20260830/P337-N170-angle-flip/score.json \
  --markdown results/server-20260830/P337-N170-angle-flip/REPORT.md

python3 -m unittest discover -s tests -p 'test_score_n170_angle_flip.py'
```
