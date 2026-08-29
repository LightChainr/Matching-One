# Issue #159: norm-2 mechanism split after the Pell character result

Status: post-reveal theory split with prospective N60/N112 predictions.  The
N30/N56 means were known before this note; no child result was inspected.

## The negative result is narrower than an H4 rejection

The real primitive-sector contrast is a nontrivial `C3` character, not a
modular scalar.  This matters at the equianharmonic point.  The scalar
ordinary-Q4 one-point is proportional to `E4` and must vanish at `rho`, but
the character construction in Issue #156 was introduced precisely so that a
spin-4 response can survive in a nontrivial stabilizer representation.

Consequently, the same signs of `C30` and `C56` reject the extra assumption

```text
C(tau) proportional to the scalar E4(tau) amplitude,
```

not spin 4 in the character channel.  The failed historical `-2` score remains
a valid rejection of that scalar-E4 bridge.  It should not be promoted into a
general H4 falsification.

## Radial power plus the dimension--spin bound

For a first-order local correction of scaling dimension `x`, a dimensionless
torus probability has

```text
L^(2-x) = N^(-(x-2)/2).
```

The observed fixed-coordinate values `N*C=0.2265,0.2120` motivate, but do not
prove, an `N^-1` law.  If that law is asymptotic, it selects `x=4`.  A local
spin-8 field obeys `x>=|s|=8`, even before imposing detailed percolation
spectrum information, and therefore begins no slower than

```text
L^-6 = N^-3.
```

Thus “H8-like phase” and “leading N^-1 local field” cannot both be true.  An
`E4^2` or two-H4 composite is even under the Pell-side sign reversal, but as a
second-order response to an `N^-1` irrelevant coupling it begins at `N^-2`.
It can be a correction, not the leading explanation of a stable `N^-1` term.

The minimal live mechanisms are therefore:

| mechanism | status under `C~N^-1` | norm-2 phase | radial law |
|---|---|---:|---:|
| identity-family/rank-4 analytic anisotropy | compatible | `-` | `N^-1` |
| nonlocal topological character with an even lift | compatible but not a local CFT field | `+` | empirically `N^-1` |
| quadratic H4 / `E4^2` composite | subleading | `+` | `N^-2` |
| first-order local H8 | excluded as the leading `N^-1` term | `+` | at most `N^-3` |
| scalar-to-character Pell leakage | subleading after actual-`tau` baseline subtraction | reverses with signed shape displacement | at most `N^-2` |

The nonlocal row deliberately has no assigned scaling dimension.  Applying
the local `x>=|s|` bound to a global homology-sector functional would assume
the conclusion.  Its embedding phase must instead be measured.

## Smallest independent discriminator

Multiply both period vectors by the Gaussian integer `1+i`.  In the declared
column convention this is the exact left action

```text
M = [[1,-1],[1,1]],       P_child=M P_parent.
```

It doubles `N`, rotates the embedding by `pi/4`, leaves `tau` unchanged, and
transports the period-basis homology labels by the identity.  For
`Z=C+iQ`, the exact angular actions are

```text
H4: Z -> -Z,
H8: Z -> +Z.
```

Combining angle and radial power gives four parameter-free transfers:

```text
rank-4 H4, N^-1:       Z_child/Z_parent = -1/2
nonlocal even, N^-1:   Z_child/Z_parent = +1/2
quadratic H4, N^-2:    Z_child/Z_parent = +1/4
local H8 at x=8:       Z_child/Z_parent = +1/8.
```

The two existing parents produce the unseen children

```text
N30  [[6,3],[0,5]] -> N60  [[6,-2],[6,8]]
N56  [[8,4],[0,7]] -> N112 [[8,-3],[8,11]].
```

N60 has Smith invariants `(2,30)`, so the general-period backend is required;
this is not a cyclic-engine target.

Score the full `(C,Q)` residuals, not scalar ratios.  The leading two tests
are

```text
R4    = 2 Z_child + Z_parent,
Reven = 2 Z_child - Z_parent.
```

The archived parent and fresh child streams are independent, so their
cross-covariance is zero; retain the complete within-geometry `(C,Q)`
covariance and combine both lineages in one fixed-model score.  Only after the
four transfers are scored may a free transfer ratio be reported.

## Scientific outcomes

- A negative child selects the rank-4 lift and rehabilitates H4 in the
  nontrivial character channel without restoring the failed scalar-E4 bridge.
- A positive half-size child selects a genuinely nonlocal/even character at
  the same radial order; this would be more novel than a local H8 label.
- A much smaller positive child supports a quadratic or high-dimension
  correction, implying the apparent `N^-1` stability at N30/N56 was
  preasymptotic or mixed.
- A nonzero transformed `Q` inconsistent with the fixed two-vector action
  signals a transport/reflection problem rather than a new harmonic.

The machine-readable contract is
`predictions/p159_norm2_character_discriminator_20260829.yaml`.
