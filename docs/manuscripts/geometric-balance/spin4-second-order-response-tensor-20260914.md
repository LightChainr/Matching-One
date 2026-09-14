# Second-order spin-four response tensor after the leading matching-odd root correction

2026-09-14.

Status: symmetry/RG synthesis motivated by the deterministic safe-root controls. The representation-theory decomposition below is exact at the level of rotation characters; identifying the microscopic lattice couplings with particular CFT fields and computing their response coefficients remain open.

## 1. Why the old `cos(4theta)^2=(1+cos8theta)/2` shorthand is too strong

The repository's composite-field oracle correctly states that the product of two real aligned spin-four harmonics has support

```text
H4(theta)^2 = (H0+H8(theta))/2.
```

That is an exact statement about the elementary angular monomial. It does **not** imply that a general second-order observable response must have equal scalar and spin-eight coefficients.

The distinction matters now because post-H4 safe-root data resolve a q=2 H8 response, while the simplest locked `P0=P8` interpretation is unstable.

## 2. Complex spin basis

Write the two real spin-four perturbations as

```text
T4(theta) = t_+ e^{+4 i theta} O_T^+ + t_- e^{-4 i theta} O_T^-,
I4(theta) = i_+ e^{+4 i theta} O_I^+ + i_- e^{-4 i theta} O_I^-,
```

with reflection reality `t_-=conj(t_+)`, `i_-=conj(i_+)` for the square-lattice perturbation.

At second order, rotational covariance permits two inequivalent tensor contractions:

```text
spin 0:
  t_+ i_- R_0[O_T^+,O_I^-]
+ t_- i_+ R_0[O_T^-,O_I^+],

spin 8:
  t_+ i_+ R_8[O_T^+,O_I^+]
+ t_- i_- R_8[O_T^-,O_I^-].
```

The response functionals `R_0` and `R_8` are different integrated connected correlators/OPE channels. Symmetry does not set them equal.

For a reflection-even square orientation they reduce schematically to

```text
delta O^(2)(theta)
 = C0 * H0 + C8 * H8(theta),
```

with independent `C0,C8` unless an additional factorization/dynamical identity is proved.

Thus the elementary `1:1` coefficient occurs only in a special aligned-factorized model where the observable response treats the two tensor channels with the same coefficient. It is not a generic consequence of `4 tensor 4 = 0 plus 8`.

## 3. Matching parity and radial exponent

Let the leading matching-odd thermal spin-four coupling be

```text
T4: omega_T=13/4,
```

and an even spin-four lattice anisotropy be

```text
I4: omega_I=2.
```

Their mixed second-order response is matching odd and has total irrelevant power

```text
omega_T+omega_I = 21/4.
```

For the charge/root observable, relative to the leading T4 root term this adds two powers of length:

```text
leading T4 root: ell^-4,
T4 x I4 mixed root: ell^-6.
```

The spin-zero and spin-eight tensor channels share this radial exponent but need not share amplitudes.

## 4. A separate q=2 H4 dressing channel

The common `omega≈2` finite-size correction is not purely spin four. Existing oblique magnetic-gap controls indicate that its angular-scalar component is substantially larger than its H4 component.

Denote the even scalar part by `S0`. Then

```text
T4 x S0
```

is matching odd, remains in the H4 angular sector, and also gives a root `ell^-6` correction.

Therefore the minimal q=2 root structure is naturally

```text
p_root-pc
 = c4  H4 ell^-4
 + c46 H4 ell^-6          [T4 x S0 / radial dressing]
 + c80 H0 ell^-6          [T4 x I4, spin0 tensor]
 + c88 H8 ell^-6          [T4 x I4, spin8 tensor]
 + ... .
```

There is no reason to identify `c46`, `c80`, and `c88`.

## 5. What the new Gaussian controls actually determine

Two norm-2 Gaussian lineages are especially clean because child and parent have opposite H4 but the same H8.

### Lineage A

```text
(2,1) -> (3,1), n=2,3
H4: -7/25 -> +7/25
H8: -527/625 -> -527/625.
```

### Lineage B — frozen H8-node negative control

```text
(3,2) -> (5,1), n=1,2
H4: -119/169 -> +119/169
H8: -239/28561 -> -239/28561.
```

For each lineage, `pc` and every angular-scalar contribution cancel from the child-parent difference. Fitting

```text
D_n = p_child-p_parent = alpha n^-4 + beta n^-6
```

therefore isolates the odd H4 leading/dressing pieces plus the known small effect of the common H8 value through radial rescaling.

Combining the two lineages gives, within the finite q=2 model,

```text
c46 = -0.1156140546,
c88 = -0.1381781858,
```

without using an external `pc`.

The second lineage is close to an H8 node, so its `beta` is primarily a measurement of `c46`. Its frozen pre-target prediction agreed at the few-percent level. This provides a genuine negative control for the H8 decomposition.

What these differences do **not** determine is the scalar tensor coefficient `c80`, because scalar angular pieces cancel under Gaussian orientation differences.

## 6. Consequence: do not infer scalar q=3 from an H4 projector alone

A two-angle H4 projector retains both angular-scalar and H8 pieces. The earlier axis/(3,4) `ell=5,10` sequence happened to mimic `ell^-7` extremely well. The cross-node control with nearly opposite H8 projector weight falsifies the interpretation that this was dominantly one scalar `ell^-7` field.

The correct hierarchy is now:

```text
observed leading: c4 H4 ell^-4;
observed next:    c46 H4 ell^-6 + c88 H8 ell^-6;
unknown:          c80 H0 ell^-6 and later scalar/H8 blocks.
```

A genuine scalar `x≈33/4` block may still coexist, but it is no longer required by the two-scale H4-projected residual.

## 7. Continuum field-theory target

If `T4` is the thermal-family level-four descendant and `I4` is the leading even identity/lattice anisotropy, the next nontrivial CFT calculation is not another exponent fit. It is the pair of integrated second-order matrix elements

```text
R_0  ~ integral <m | O_T^(+4) O_I^(-4) | m>_conn,
R_8  ~ integral <m | O_T^(+4) O_I^(+4) | m>_conn,
```

including all counterterms/contact terms required by conformal perturbation theory and the actual safe/rank projector.

At c=0, logarithmic mixing/null descendants can modify either channel independently. The calculation should therefore be done first at generic Q if singular collisions obstruct the Q=1 tensor decomposition.

## 8. Modular fingerprint conjecture

In the simplest ordinary-module scenario:

- the leading thermal Q4 one-point has weight-four fingerprint `E4(tau)`;
- an ordinary same-chirality spin-eight mixed response naturally lies in weight eight, whose full-modular holomorphic space is generated by `E8=E4^2`;
- the spin-zero mixed channel belongs to a different nonholomorphic/modular-scalar structure and need not share its amplitude.

This predicts that the mixed spin-eight channel vanishes at the hexagonal elliptic point where `E4=0`, more strongly than a generic unrelated H8 defect sector need do.

This is a falsifiable shape prediction, not a theorem for the Matching-One observable until its torus/map dictionary is fixed.

## 9. Practical next tests

1. Use source-resolved Perron/Feynman--Hellmann derivatives to identify a physical source with nonzero projection on the common even spin-four correction, distinct from the exact thermal-reparametrization sources already rejected.
2. Measure the mixed odd/even source Hessian and project it into spin0 and spin8 angular channels.
3. In parallel, test the weight-eight/hexagonal-zero fingerprint on existing or bounded Pell/modulus assets.
4. Only after `R_0/R_8` is theoretically or empirically constrained should a scalar post-H4 block be inferred from a two-angle H4 projector.

The machine-readable finite controls are in `gaussian-h8-node-q2-control-20260914.json` and `cross-node-post-h4-projector-20260914.json`.
