# Post-H4 cross-node discriminator: a bounded replacement for same-ell multi-angle tomography

2026-09-14.

Status: deterministic finite safe-transfer analysis plus mechanism triage. This note is explicitly corrective: the earlier plan to use N1105 same-ell multi-angle safe-root tomography is not executable under the current row-memory transfer and, even at zero numerical noise, fixed-ell radial dressing confounds asymptotic angular coefficients. The exact Gaussian-circle projector algebra remains valid; the production route is replaced here by a low-memory pair of nearly equal physical circumferences.

## 1. What is being corrected

The previous same-projector axis/(3,4) values at `ell=5,10` suggested

```text
p_perp4(ell)-pc ~ C7 ell^-7,
```

because multiplying the two residuals by `ell^7` nearly froze them. That was a useful hint but it was not an angular identification: the same two-angle projector has a large surviving H8 coefficient. With only two scales, an H8 contribution can masquerade as a scalar radial power.

The later independent state-cost/confounding audit also rules out buying enough same-ell orientations to solve H0/H8/H12 directly with the current transfer algorithm. We therefore need a different geometry that changes the surviving angular coefficient while leaving the physical scale nearly unchanged.

## 2. A complementary node pair with low row memory

Use

```text
u_A=(5,2), n_A=2, ell_A^2=116,
nu_B=(3,2), n_B=3, ell_B^2=117.
```

The physical circumferences differ by only about 0.43 percent.

Their exact harmonics are

```text
H4_A =  41/841  = +0.0487514863...,
H8_A = -703919/707281 = -0.9952465852...,

H4_B = -119/169 = -0.7041420118...,
H8_B = -239/28561 = -0.00836805434....
```

Thus A is close to an H4 node while B is close to an H8 node. They are almost complementary filters.

The current lifted-gain safe transfer remains tractable:

```text
(5,2),n=2:  G4 22,994 states; G8 131,677 states; row memory 5/7.
(3,2),n=3:  G4 19,018 states; G8 206,197 states; row memory 3/5.
```

This is qualitatively different from the millions-to-astronomical-state same-circle designs.

## 3. Tight-root control

The new calculation uses sparse Perron tolerance `1e-12` and root tolerances around `1e-14`, rather than the looser exploratory settings in the original oblique script.

As an implementation regression, the same code gives

```text
(3,2), n=2:
local      0.5928240591852450
repository 0.5928240591851011
```

for a difference of about `1.4e-13`.

The tight roots used here are

```text
p_A = 0.5927450615721411,
p_B = 0.5927612086097266.
```

## 4. Remove the leading H4 term without using pc

Write the leading root law as

```text
p_i = pc - A4 H4_i ell_i^-4 + residual_i.
```

Let

```text
q_i=H4_i/ell_i^4.
```

Choose weights `w_A+w_B=1` and `w_A q_A+w_B q_B=0`. Numerically,

```text
w_A=0.934200379537736,
w_B=0.0657996204622639.
```

The leading-H4 projected root is therefore

```text
p_cross^perp4
 = w_A p_A+w_B p_B
 = 0.5927461240410858.
```

The simultaneously inferred leading H4 amplitude is

```text
A4_hat=0.2932542820.
```

Because the two circumferences are not exactly equal, a same-H4 `ell^-6` dressing is not killed identically. But using the already observed `A4(ell)=Ainf+A2/ell^2` scale with `A2~0.66`, its leakage into this projector is only about `1.6e-10`, far below the `1e-7` effects below.

## 5. The key angular fact

For the projected pair, before tiny radial-mismatch corrections, the surviving angular weights are approximately

```text
H8:     -0.93031,
H4^2:   +0.034845.
```

Compare the existing exact-equal-ell axis/(3,4), ell=10 projector:

```text
H8:     +0.6864,
H4^2:   +0.8432.
```

Using the common diagnostic threshold `pc=0.5927460507921` only after root location,

```text
axis/(3,4), ell=10:   p_perp4-pc = -7.57256e-8,
cross-node, ell~10.79: p_perp4-pc = +7.32490e-8.
```

The two post-H4 residuals reverse sign exactly when the H8 projector weight reverses sign, while the scalar weight remains positive and the H4^2 weight remains positive.

This makes the earlier statement

```text
post-H4 residual is primarily a scalar V_<1,4>-like ell^-7 term
```

unsafe. A substantial H8-like angular component is required unless one invokes a finely tuned multi-channel cancellation.

## 6. Reference-free stress tests of the old scalar interpretation

The earlier axis/(3,4) projector has the same angular weights at `ell=5` and `ell=10`. Fit those two points to a pure scalar law

```text
p_perp4(ell)=pc_fit+C7 ell^-7.
```

Because both `pc_fit` and `C7` are fit from those two roots, this test does not use an external threshold. It gives

```text
pc_fit = 0.5927460516782668,
C7     = -0.7661178948.
```

Transporting that scalar law to the cross-node projector predicts

```text
0.5927460061955850,
```

whereas the observed projected root is

```text
0.5927461240410858.
```

The miss is about

```text
1.18e-7,
```

roughly three orders above the estimated same-H4 radial-dressing leakage. Thus the old two-scale scalar story fails a new angularly changed held-out geometry.

## 7. A caution: H4^2 mixing can mimic the new geometry

There is an important counterpoint. Fit the same `ell=5,10` axis/(3,4) roots to

```text
p_perp4=pc_fit+Cmix P_perp4[H4^2] ell^-6.
```

This gives

```text
pc_fit = 0.5927461295061164,
Cmix   = -0.1831589648.
```

It predicts the cross-node projected root as

```text
0.5927461255146714,
```

only `1.5e-9` from the observed value.

So the cross-node result does **not** by itself identify an independent H8 field. A mixed term involving the leading odd H4 correction and an even lattice anisotropy can generate H4^2=(H0+H8)/2 and reproduce this geometry extremely well.

However, that pure three-point mixed fit places its asymptotic intercept about `7.9e-8` above the independent diagnostic pc. With the diagnostic pc held fixed, the larger-scale H4-projected data are described materially better when an H8-like column is included than by scalar or H4^2 columns alone. Therefore the safest current statement is:

> post-H4 data contain a large H8-like angular component; whether it is an independent dual-odd H8 operator or nonlinear mixing must be decided by a source/field calculation, not by naming the residual from two root powers.

## 8. Why the common x~4 correction does not automatically settle the mixing question

The oblique magnetic-gap controls give, roughly,

```text
(ell I0-2pi*5/48) ell^2 ~ 0.29--0.34.
```

A crude constant-plus-H4 regression gives a dominant scalar part and a much smaller H4 part. Therefore the observed common `omega~2` correction should not automatically be called an even spin-4 field. If it is mainly scalar, multiplying it by odd thermal Q4 only dresses the H4 sector and is removed by an H4 projector; only its spin-4 subcomponent generates H4^2/H8.

This is a concrete reason not to promote the q=2 mixed explanation solely because it is symmetry-allowed.

## 9. Current candidate ordering

After the corrections above, I would not order candidates by a two-point radial exponent. The useful ordering is by what remains to be demonstrated:

1. **H8-like post-H4 angular response — observed finite evidence.** The angular component is now the most robust new statement.
2. **Nonlinear odd-H4 x even-spin4 mixing — live mechanism.** It can reproduce the new cross-node geometry, but the relevant even spin-4 coupling and asymptotic intercept must be established.
3. **Independent dual-odd H8 operator — live mechanism.** Linear identity-family spin8 is one candidate only if its matching-odd coupling is shown nonzero; logarithmic/interchiral alternatives remain possible.
4. **Scalar x~33/4 block — demoted from default.** It may coexist, but the pure scalar ell^-7 explanation fails the cross-node held-out geometry.
5. **H12/higher channels — retained as residual alternatives, not currently required by this control.**

This does not change the well-supported leading result: the first non-common root correction is H4/spin4 with radial dimension consistent with x=21/4.

## 10. Next high-information analysis, no large circle

The next calculation should target the origin of the H8-like component rather than another angular census.

Two routes are higher value than N1105:

1. **source-resolved mixed derivative:** in the safe Perron/Clapeyron framework, insert one physical source known to excite the common `omega~2` correction and differentiate the odd-H4 root response with respect to it. A genuine `u4^- v4^+` mechanism predicts an H4^2 second response with a locked sign/tensor structure.
2. **rank-projected Ward/generic-Q test:** compute whether the actual safe/rank functional preserves the ordinary thermal-Q4 Ward reduction or acquires logarithmic/null/seam terms. An independent H8 or scalar block should then appear as a genuinely new continuum matrix element rather than as ordinary Q4 dressing.

The machine-readable inputs and model checks are in `results/research-control-20260914/cross-node-post-h4-projector-20260914.json`.
