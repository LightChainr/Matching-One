# A genuine microscopic odd tangent in the self-matching C4 theory

The doubled exchange from #233 collapses to an internal finite-graph operator in the self-matching C4 family. On the same configuration space,

```text
C|A> = |A^c>,
C W(t,lambda) C = W(-t,-lambda),
C^2=I.
```

This is stronger than the non-self-matching `C4 -> K4` toy: no graph-species label is needed. Differentiating at the center produces two genuinely local, single-theory odd scores,

```text
S_t      = 4 sum_x (n_x-1/2),
S_lambda = 4 sum_x parity(x)(n_x-1/2).
```

The `(3,1)`, `N=10` quotient has five sites of each parity. Its exact Fisher Gram matrix is

```text
<S_mu S_nu> = [[40,0],[0,40]].
```

Thus `S_lambda` is already a microscopic complement-odd tangent exactly orthogonal to the uniform thermal score in the UV probability metric. The possible no-go “only a doubled defect exists, with no local odd field in one theory” is false for this self-matching family.

Fisher orthogonality is not yet RG eigen-orthogonality: the two checkerboard sublattices are inequivalent, so coarse graining may mix the two odd scores. The finite exact oracle nevertheless allows a sharper zero to be built.

Let `O_local` be the radius-one complement-odd pivotal H4 readout from #219, and let

```text
epsilon_cell = [(n_even-1/2)+(n_odd-1/2)]/2
```

on one neighboring even/odd microscopic cell. Its response is exactly `(d_t,d_lambda)=(1,0)`, while

```text
O_local: (-3/64, 11/64).
```

Therefore the one-parameter local family

```text
O_alpha = O_local + alpha epsilon_cell
```

has a unique thermal zero

```text
alpha* = 3/64,
(d_t,d_lambda) O_alpha* = (0,11/64).
```

This is an exact local matrix element, not a fit. Complement oddness holds configurationwise. It supplies a cleaner large-size readout than subtracting a noisy global wrapping response: accumulate `O_local`, `epsilon_cell`, `S_t`, and `S_lambda` in the same stream and keep `alpha=3/64` frozen.

There is a dual coupling-space zero as well. The exact global response row `(15/8,5/4)` is annihilated by

```text
(delta t,delta lambda)=(2,-3),
```

but the local H4 matrix element along this direction is `-39/64`. In site probabilities this direction is `(delta p_even,delta p_odd)=(-1,5)` up to an overall infinitesimal scale. It is a microscopic odd perturbation invisible to the N=10 global thermal proxy yet visible locally.

These zeros do not prove an irrelevant RG eigenoperator or `x=21/4`. They prove the exact local observability and remove the UV thermal component in two complementary ways. They also clarify the sector bookkeeping: complement parity here is internal `Ad_C` parity and the operators are Potts-colour singlets, so #257's `[2]` selection zero is independent.
