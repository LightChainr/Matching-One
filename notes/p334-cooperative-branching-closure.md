# P334/P429: cooperative branching closure from the same checkpoint

This is a fresh dynamic continuation experiment, not another static geometry
proxy. Its parent is P429 production commit `751f8b3`.

The zero-new-sample follow-up now supplies a stronger exact result:
[real checkpoints with identical age, line, H2 and b2 have different branching
probabilities](p334-real-checkpoint-scalar-nonclosure.md). The original pilot
below establishes within-checkpoint heterogeneity; the follow-up separately
establishes insufficiency of that matched scalar state.

For a rank-one checkpoint C with d vacant sites, let b1 be the number of safe
one-site insertions and b2 the number of safe unordered two-site sets. For each
safe first site v let c_v count safe second sites. Monotonicity gives the exact
incidence identity `sum_v c_v=2*b2`. Set `s1=b1/d`,
`s2=2*b2/[d(d-1)]`, and `q_v=c_v/(d-1)`.

One common uniform update followed by two independent one-site clones has
conditional success `sum_v c_v^2/[d(d-1)^2]`. Uniform-safe-successor closure
would instead give `s2^2/s1`. Their exact difference is

```
Delta(C) = [b1 * sum_v c_v^2 - (2*b2)^2] / [d*b1*(d-1)^2]
         = s1 * Var(q_v | safe v, C) >= 0.
```

Set Delta(C)=0 when b1=0. The numerator is evaluated in exact integer
arithmetic, so positivity is not floating-point subtraction. The b2 census
already requires all c_v; saving their squared sum has no extra topology cost.
The primary ensemble estimand is the average of Delta(C), conditional on
rank-one C. This Rao–Blackwellizes the sampled common-update statistic without
changing its expectation.

Every row also retains H2, b2, the selected exact q_after and q_after^2, and
Y1/Y2/Y1Y2 from the two independently tagged clones. Their calibration rows
are evaluated jointly. A base permutation (including both paired orientations)
is one checkpoint cluster; the clones are never two sample units.

Pre-data decisions are in `analysis/p334_cooperative_closure_freeze.json`:
N325/Tg and N425/XP, 20k base permutations each, nominal alpha .01 two-coordinate
Wald gate with both coordinates positive, and no automatic extension.

The narrow eliminated mechanism, if positive, is uniform safe-successor
continuation at these finite checkpoints. It is not a proof of path memory,
nonlinear H2/b2 insufficiency, a scale law, or a continuum field identification.

## Result: cooperative excess survives checkpoint conditioning

Both fresh 20k runs finished on the local Mac, in 40.75 s and 69.34 s wall time.
The pre-data host override is `analysis/p334_cooperative_closure_host_override.json`;
the runner remains the source frozen at `6712ec5`. No earlier P429 clone-gap
sample is pooled. All 35,954 at-risk checkpoint rows are archived.

| size / orientation | at-risk C | exact-positive C | Delta_coop | cluster SE | z |
|---|---:|---:|---:|---:|---:|
| N325 first | 8,997 | 8,987 | 0.0001188143 | 0.0000021576 | 55.07 |
| N325 second | 8,966 | 8,951 | 0.0001224303 | 0.0000024034 | 50.94 |
| N425 first | 8,910 | 8,895 | 0.0000769773 | 0.0000015008 | 51.29 |
| N425 second | 9,081 | 9,072 | 0.0000754424 | 0.0000014499 | 52.03 |

The two-coordinate nominal Wald statistics are 5701.49 and 5335.75, far beyond
the pre-data 9.210 gate. The estimate is an average of nonnegative exact
checkpoint variances, so observed integer-positive numerators already refute
exact uniform-safe-successor closure on those configurations. The large Wald
numbers quantify the ensemble magnitude; they are not a universal exponent test.
Both sizes pass; **stop at 20k, no expansion**.

Primary checkpoint-cluster covariance matrices, first/second order, are

```
N325: [[ 4.6552806563e-12, -6.7653722136e-14],
       [-6.7653722136e-14,  5.7764125732e-12]]
N425: [[ 2.2524429907e-12,  9.3251174202e-16],
       [ 9.3251174202e-16,  2.1022689878e-12]]
```

The full 22-by-22 joint covariance within each size includes all exact and
sampled bridge rows in `results/local-20260831/P334-cooperative-closure/score.json`.
All exact row identities pass. The largest absolute z over common-update and
clone calibration rows is 1.66; no clone is counted as an independent checkpoint.

### What the former pooled gap was mixing

Write `S1=E_C[s1]`, `S2=E_C[s2]`, `B=E_C[s2^2/s1]`, and
`Q=E_C[E_common(q_after^2|C)]`. The exact current-archive decomposition is

```
Q-S2^2 = (1/S1-1)*S2^2 + (B-S2^2/S1) + (Q-B)
           common gate   between C     within C = Delta_coop.
```

| size / orientation | pooled gap | common gate | between checkpoints | within checkpoint |
|---|---:|---:|---:|---:|
| N325 first | 0.05266772 | 0.05125541 | 0.00129349 | 0.00011881 |
| N325 second | 0.05333549 | 0.05191219 | 0.00130087 | 0.00012243 |
| N425 first | 0.04556983 | 0.04453738 | 0.00095547 | 0.00007698 |
| N425 second | 0.04523541 | 0.04420293 | 0.00095703 | 0.00007544 |

After removing the common gate, the genuinely within-checkpoint term is
7.31–8.60% of the remaining gap. Thus the earlier conditional-on-common-safe
signal was still mostly **between-checkpoint hazard composition**. Nevertheless,
the within-checkpoint continuation heterogeneity is small, ubiquitous, and
resolved; it cannot be attributed to pooling different starting checkpoints.
This is the new mechanism-changing distinction, not just a stronger pooled z.

## Scientific card

- Changed mechanism space: eliminates exact uniform-safe-successor cooperative
  closure at these finite rank-one checkpoints, directly in the branching language.
- Not proved: path memory, full H2/b2 state insufficiency, an exponent, or CFT identity.
- Observer/sector: common update then two independent one-site survival clones;
  occupied NN ambient-H1 rank-one plateau, absorption at rank two.
- Source/geometry: fresh 20k counter permutations at N325 k0=193 and N425 k0=252;
  paired Gaussian quotients (17,6)/(18,1) and (16,13)/(19,8), exact HNF backend.
- Dependency groups: `p334-cooperative-N325-20260831` and
  `p334-cooperative-N425-20260831` are independent across sizes, paired across
  orientations, exact counts and both clones within each base checkpoint.
- Next discriminant: a branching-aware state must reproduce the measured
  within-checkpoint successor-degree dispersion, rather than matching only
  single-chain survival and the common absorption gate. This result alone does
  not prescribe or identify that richer state.

Reproduce the score with the explicit research Python environment:

```bash
/Users/lc/python-envs/research-py311/bin/python scripts/score_p334_cooperative_closure.py \
  --freeze analysis/p334_cooperative_closure_freeze.json \
  --n325-prefix results/local-20260831/P334-cooperative-closure/raw/N325 \
  --n425-prefix results/local-20260831/P334-cooperative-closure/raw/N425 \
  --runner-commit 6712ec5b00bdf3dc0c6f8733ef85eda58b86cb2f \
  --output results/local-20260831/P334-cooperative-closure/score.json
```
