# F2/F3 flat twists from the N65 projective-birth archive

Status: existing-data bridge from the exact finite-abelian transform in #337
(`4b17955`) to the N65 marked filtration in #334 (`1714141`).  No new field,
counter block or server run is needed.

## Why the archive is sufficient

At a fixed occupation count `K`, each archived path determines

```text
K < tau1                 -> rank 0,
tau1 <= K < tau2         -> rank 1 with primitive line ell,
tau2 <= K                -> rank 2.
```

The `DIRECT_RANK2` atom has `tau1=tau2`, so it simply has no rank-one interval.
By the #269 saturation theorem, a primitive integer line never reduces to zero
in `F_q^2`.  Therefore `ell mod q` determines the unique projective rank-one
class for every prime `q`.

For a nonzero flat twist `alpha`, the exact #337 transform is

```text
T_alpha = P0 + L_ker(alpha),
```

while `T_0=1`.  The N65 archive supplies `P0` and every `L_line`, hence all
four F2 sectors and all nine F3 sectors at arbitrary fixed `p`.

## Exact gates

Every one of the 20 batches and both orientations passes

```text
S_q = sum_alpha T_alpha = q^2 P0 + q P1 + P2,
```

the individual line recovery `T_alpha-P0=L_ker(alpha)`, `T_0=1`, and the
order-2/order-3 source inversion

```text
P0=(S3-2S2+1)/2,
P1=S2-1-3P0,
P2=1-P0-P1.
```

The largest floating residual is `2.66e-15`.  Thus the finite-twist source is
not merely an exact pocket oracle: the current production-format event archive
already implements it.

## Balanced F3 H4 alias

`P^1(F3)` has four lines split into two D4 orbits:

```text
axes:     (1,0), (0,1),
diagonal: (1,1), (1,-1).
```

The minimal balanced projector is

```text
C_F3,H4
 = 1/2[(T_kernel-axis0 + T_kernel-axisInf)
       -(T_kernel-diag+ + T_kernel-diag-)],
```

where the common `P0` cancels.  It is exactly half the rank-one
axes-minus-diagonals probability.

For comparison, F2 has three projective lines and the unit-norm zero-sum
projector

```text
C_F2,H4=(T_axis0+T_axisInf-2 T_diag)/sqrt(6).
```

At `p_ref=0.592746050790`, the shared-batch same-modulus orientation contrasts
are

| projector | contrast | SE | absolute z |
|---|---:|---:|---:|
| F3 balanced H4 | 0.00110298 | 0.00150621 | 0.732 |
| raw phase-aligned chi4 | 0.00219347 | 0.00300739 | 0.729 |
| F2 H4 alias | 0.00090614 | 0.00163771 | 0.553 |

So F3 is nominally the sharpest H4 estimator, but only by about 0.4% in z:
it is effectively tied with raw `chi4`, not a material variance reduction.
Its value is semantic and representational: it is the smallest balanced flat-
twist projector and automatically removes `P0`.

## A sharper sector that raw chi4 discards

The full F3 character basis also retains a reflection-odd diagonal split

```text
C_F3,diag-odd = [L_(1,1)-L_(1,-1)]/sqrt(2).
```

Its same-modulus orientation contrast is

```text
0.00207756 +/- 0.00104063,  |z|=1.996.
```

Together with the F3 axis-odd contrast, its covariance-aware diagnostic is
`5.903 / 2 df`.  This is sharper than raw H4 `chi4`, but it is a different
reflection-odd projective sector, not a better estimator of the same H4
amplitude.  That distinction is the scientific gain: the flat-twist vector
retains modular line information that a single physical `chi4` contraction
collapses.

## Claim boundary and next discriminator

- All sector reconstructions and source identities are exact consequences of
  the archive plus #269/#337.
- The `z` rankings reuse the exploratory P334 20k block and are not independent
  tests.
- No affine-TL/Potts defect identity or continuum field is assigned.
- A fresh block should freeze both the balanced F3 H4 alias and the F3
  diagonal-odd null/sector before scoring.  The latter is the genuinely new
  channel; the former is primarily the canonical twist encoding of H4.

Reproduce with:

```bash
python3 scripts/score_flat_twist_projective_archive.py \
  --births results/local-20260830/P334-projective-birth-N65-smoke/n65_20k.births.csv \
  --metadata results/local-20260830/P334-projective-birth-N65-smoke/n65_20k.metadata.json \
  --json results/local-20260830/P334-projective-birth-N65-smoke/flat_twist_score.json \
  --markdown results/local-20260830/P334-projective-birth-N65-smoke/flat_twist_score.md

python3 -m unittest discover -s tests \
  -p 'test_flat_twist_projective_archive.py'
```

