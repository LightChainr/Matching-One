# Which orientation triggers explains the main next-label mode

The exact vacancy census turns the former three-mask readout into a genuine
within/between decomposition. The main new result is that **the two-bit trigger
type carries most of the positive first/completion next-label covariance**.
The coarser safe-versus-any-birth partition obscures this direction by pooling
the two different orientation triggers.

## Fixed source and estimand

This is an analysis of the existing fork block `e32a8593`, 20,000 original
prefixes and 640,000 suffixes per size, with the deterministic census
`e9dc7a1078b2c64b319f4a36ffc1c844e8426aa0`. No new random data, DP or engine replay
was run. The counter, geometry, canonical p and paired H4 normalization are
unchanged. All errors use the original twenty common batches.

For a fixed prefix Z, let m(u)=E[X|Z,u], and let c(u) record whether each
orientation changes rank at the next site. Thus c is 00, 01, 10 or 11; this
records any rank increase, including a direct jump. If d sites are vacant,
the exact class probabilities are

```
p00 = joint_safe_count / d
p01 = (first_safe_count - joint_safe_count) / d
p10 = (second_safe_count - joint_safe_count) / d
p11 = (d - first_safe_count - second_safe_count + joint_safe_count) / d.
```

For the unbiased quartet matrix Bhat, a same-class mask has expectation
`E[1{c(U)=c(V)=c} Bhat|Z] = pc^2 Bc`, where Bc=Cov(m(u)|Z,c).
Dividing this mask by the **individual prefix's exact pc** estimates pc Bc.
The empty class contributes zero. Hence

```
Bwithin  = E_Z sum_c pc Bc
Bbetween = Btotal - Bwithin
         = E_Z Cov_c(E[m|Z,c]).
```

For the two-class safe/birth partition the same formula applies. Its mixed
mask is `ps pb (Bs + Bb + (mu_s-mu_b)(mu_s-mu_b)^T)`, not just between-class
covariance. This corrects the tempting interpretation of the original roughly
78% mixed-mask share.

## Main mechanism readout

Gamma is the F1/F2 cross entry of paired next-label B. Values below are means
plus/minus shared-batch SE; canonical p is 0.59274605079.

| N | total Gamma | binary between | four-type between | four-type within |
|---|---:|---:|---:|---:|
| 325 | .001381329 ± .000044623 | .000173427 ± .000082448 | .001297665 ± .000047892 | .000083664 ± .000039240 |
| 425 | .000997369 ± .000040007 | .000055687 ± .000029454 | .000873098 ± .000037316 | .000124271 ± .000035827 |

The four-type between direction is 93.94% / 87.54% of the canonical Gamma point
estimate. For the 01+10 checkpoint pair it is .000622660 ± .000026962 /
.000401556 ± .000025355, or 96.15% / 86.32% of that group's total Gamma.
The binary between direction of this group instead has negative Gamma,
−.000285200 ± .000074051 / −.000212625 ± .000031013. A between-class covariance
matrix can have either sign in this cross entry: the negative values are not
negative variance.

The integrated four-type between Gamma is
`2.0203453e-5 ± 1.01534e-6` / `9.2340258e-6 ± 8.84995e-7`, with point shares
89.21% / 76.99%. The remaining within direction is
`2.4443888e-6 ± 1.59768e-6` / `2.7595124e-6 ± 8.22621e-7`.

This is a positive explanation of the coarse mode: **which orientation triggers
is the large informative discrete coordinate**, while a smaller response can
survive within its level sets. Four-type within-11 Gamma is structurally zero:
if either orientation was R0 its first birth is then fixed at k0+1, and otherwise
its first birth was already determined by the prefix. Paired F1 is constant
over that class.

The full A/E next-label variances tell the same coarse story. At canonical p,
four-type between A variances are .013946343 ± .000182669 and
.009509054 ± .000149303; between E variances are .008755683 ± .000104222 and
.006016663 ± .000133467. Corresponding totals are .014703299/.010052247 for A
and .009177981/.006062771 for E. The saved score also gives every requested
all/01+10, canonical/integrated A/E and Gamma component with common covariance.

## Scope and reproducibility

This partitions **next-label conditional-mean covariance**, not total suffix
variance, total production noise or the mean topology observable. It does not
assert that trigger type is a sufficient state: the within part remains an
estimable target, and the separate same-rank/same-degree contact response probes
a different, own-orientation restriction. Signed Gamma shares are not
probabilities. Finite-sample matrices, including small negative diagonal
estimates, are left unclipped. No joint inverse, multiple-model scan or new
independent evidence block is introduced.

Run:

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_exact_birth_partition_response.py \
  --census-commit e9dc7a1078b2c64b319f4a36ffc1c844e8426aa0 \
  --census-template 'results/p334-exact-next-label-safe-census/N{N}/N{N}.csv.gz'
```

Outputs: `results/p334-exact-birth-partition-response/score.json`, `REPORT.md`,
and a low-rank common covariance factor per N. Old masks, unnormalized
single-orientation blocks and prior covariance coordinates remain available.

Scientific card: the mechanism space changes from undifferentiated positive
Gamma to a measured trigger-type mode plus internal remainder. Observer:
paired F1/F2 and complete A/E; sector: all prefixes and named 01+10; dependency:
the original `e32a8593` fork block with exact census, not a new replicate. The
next physically distinct coordinate is the already measured local contact
response within fixed immediate rank and Euler degree, not more trigger classes.
