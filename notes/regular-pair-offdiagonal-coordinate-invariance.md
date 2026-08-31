# The distinct-site Q Hessian is not a coupling-coordinate artefact

The [frozen joint experiment](../analysis/p337_regular_pair_joint_contract.json)
uses a linear tensor parameter at every vacant vertex. This fixes what
`J2=d_logQ d_epsilon^2 U` means. There is a sharper, local-coordinate
description of the same proposed mechanism test.

## 1. Keep separate local parameters before taking the uniform path

For a fixed occupation, divide the marked colour weight by the unmarked
weight and write

```
F_A(Q,lambda)=1+sum_x lambda_x beta_x(Q,A)
               +sum_{x<y} lambda_x lambda_y beta_xy(Q,A)+... .
```

Each vertex is linear in its own lambda. For canonical Kreg every
nonempty insertion vanishes at Q1. Hence the first Q jet of the log is

```
partial_logQ log F_A|Q1
  = sum_x lambda_x a_x(A)
    + sum_{x<y} lambda_x lambda_y g_xy(A)+higher distinct-site terms.
```

Products of separately closed marks start at order (Q-1)^2 and do not
enter this jet. No same-site quadratic term is present. With the
linear original-U response functional L, the first-Q off-diagonal
Hessian of U is consequently `H_xy=L[g_xy]` for x!=y. For paired
geometries the source packet includes the corresponding contraction
in each geometry; arbitrary vertex labels do not affect the sums.

On the declared path `lambda_x=epsilon/N`,

```
J2 = (1/N^2) sum_{x!=y} H_xy.
```

Thus the proposed total and nonadjacent decisions are sums of specified
off-diagonal Hessian entries, not a newly chosen nonlinear source.

## 2. What reparameterization can and cannot remove

Under independent local changes `lambda_x=f_x(eta_x)` with f_x(0)=0
and f_x'(0)=1, mixed derivatives at distinct sites obey

```
partial_eta_x partial_eta_y = partial_lambda_x partial_lambda_y,
                                                     x!=y, at zero.
```

The f_x'' terms are diagonal only. Therefore H_xy and its fixed adjacent
or nonadjacent sums are invariant under these local unit-Jacobian
coordinate changes. Adding an explicitly same-site quadratic tensor
counterterm also cannot cancel an off-diagonal entry. A thermal source
proportional to 1 or total occupation is annihilated by the original
moving-root functional and supplies no escape either.

In contrast, an unrestricted **uniform** relabeling epsilon=f(eta)
changes the full path derivative by

```
J2_path,new = J2_path,old + f''(0)*W,   f'(0)=1,
W=d_logQ d_epsilon U.
```

It adds a diagonal/single-site path acceleration. It does not remove
the distinct-site Hessian of the declared interaction. A single
nonzero scalar path derivative with no microscopic parameter convention
would not establish this distinction.

## 3. The stronger additive null and its remaining scope

Any first-Q effective occupation action of the form
`sum_x A_x(eta_x;A)`, with no dependence on other sites' parameters,
has zero off-diagonal parameter Hessian even if each A_x is nonlinear.
Thus a nonzero fixed J2, or fixed nonadjacent part, also excludes this
site-parameter-separable first-Q description of the declared tensor
family. The word separable refers to **parameter dependence**: A_x may
still depend on the entire unperturbed occupation. A description that
depends only on a common externally tied parameter is not such a
multivariate closure.

This is a pre-readout algebraic prediction, not a claim that either
computed sum is nonzero. A nonzero sum cannot locate every H_xy entry,
identify a continuum field, or exclude models with genuine multisite
terms. The experiment retains its original contract and fixed source;
no extra data, coupling fit or post-readout contact definition is needed.
