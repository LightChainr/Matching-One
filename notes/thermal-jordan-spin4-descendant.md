# Thermal Jordan cell at level 4: a conditional LCFT mechanism for the P48 `S'` drift

Status: theory/design note written **before fresh N=185/265 P48 targets are scored**. It does not change the original P48 pure-power primary or the already-frozen `q=2` secondary correction.

Frozen tertiary artifact: `predictions/p48_sprime_jordan_log_20260828.yaml`.

## 1. External LCFT input

Vasseur, Jacobsen and Saleur showed that critical percolation at `Q=1` contains logarithmic mixing between the energy operator and a field creating two propagating clusters. In two dimensions the colliding bulk scaling dimension is `x=5/4`; the resulting pair is a logarithmic/Jordan structure rather than two independently diagonal scaling fields.

Reference:

- R. Vasseur, J. L. Jacobsen, H. Saleur, *Logarithmic observables in critical percolation*, arXiv:1206.2312, J. Stat. Mech. L07001 (2012).

For the general Virasoro representation-theory framework of rank-2 non-diagonalisable `L0` modules, see:

- K. Kytola, D. Ridout, *On Staggered Indecomposable Virasoro Modules*, arXiv:0905.0108, J. Math. Phys. 50, 123503 (2009).

This literature input makes a same-power logarithm in the thermal family a concrete LCFT mechanism, not merely an empirical curve-fitting option.

## 2. Jordan structure is inherited by a non-null descendant

Let a chiral rank-2 Jordan pair obey

```text
L0 |phi> = h |phi>,
L0 |psi> = h |psi> + lambda |phi>.
```

Let `U_n` be any homogeneous Virasoro descendant operator of total level `n`, so

```text
[L0,U_n] = n U_n.
```

Then

```text
L0 U_n|phi> = (h+n) U_n|phi>,
L0 U_n|psi> = (h+n) U_n|psi> + lambda U_n|phi>.
```

Therefore the Jordan action persists at level `n` whenever the relevant descendant is non-null in the physical quotient. The same statement applies in the bulk to left/right descendants with the scaling dimension shifted by the total descendant level.

For the percolation thermal collision,

```text
h = hbar = 5/8,
x_primary = 5/4.
```

A chiral level-4 / antichiral level-0 quasiprimary and its reflected partner have

```text
spin s = +/- 4,
x = 5/4 + 4 = 21/4.
```

Thus, conditional on the level-4 state surviving the relevant percolation module quotient, the same candidate that gives the observed spin-4 radial law can itself sit in a rank-2 Jordan pair.

## 3. Finite-size consequence

A torus correction produced by a bulk field of dimension `x` scales as `L^(2-x)`. Since the Gaussian-integer site count is `N proportional to L^2`, `x=21/4` gives

```text
L^(2-21/4) = L^-13/4 = N^-13/8.
```

The Jordan partner permits an additional logarithm at the *same* power:

```text
N^-13/8 [A + B log L]
= N^-13/8 [A_tilde + (B/2) log N].
```

A constant rescaling of the microscopic definition of `L` only shifts the non-log amplitude, so the existence of the logarithmic term is invariant under that convention.

## 4. Matching parity and the derivative spectrum

This step is conditional on the same unresolved assumption recorded in #61: the lattice matching exchange must have a well-defined scaling-limit action that commutes with dilatations on the relevant generalized eigenspace.

If it does, an involution is diagonalisable and can be resolved into matching-parity sectors. A rank-2 Jordan cell can then lie inside a fixed matching parity because the nilpotent part of dilatation commutes with the involution.

For the proposed matching-odd thermal spin-4 sector (`eta=-1`), write a unified projected scaling function schematically as

```text
F_N(z) = N^-13/8 [ f(z) + log(N) g(z) + ... ],
z proportional to (p-pc) N^(3/8).
```

The matching-even/odd observable combinations select the odd/even pieces in `z`:

```text
D_N(z) = N^-13/8 [ f_even(z) + log(N) g_even(z) + ... ],
S_N(z) = N^-13/8 [ f_odd(z)  + log(N) g_odd(z)  + ... ].
```

At the intrinsic center,

```text
P4[D](0)  ~ N^-13/8 [ f(0)  + log(N) g(0)  ],
P4[S'](0) ~ N^-5/4  [ f'(0) + log(N) g'(0) ].
```

The derivative power follows because `dz/dp ~ N^(3/8)`.

This gives an important interpretation guard:

> a small or null Jordan coefficient in the *central* `D` channel does not imply a small coefficient in `S'`.

The two observables probe different Taylor coefficients, `g(0)` and `g'(0)`, of the same logarithmic scaling function. The successful pure-power Gaussian-doubling tests of central `D` therefore constrain `g(0)` but do not by themselves rule out a derivative-specific log visible in P48 `S'`.

A future stronger LCFT test should measure enough `p` dependence to reconstruct more than one Taylor coefficient of the same projected scaling function, rather than fitting unrelated logarithms channel by channel.

## 5. Retrospective P48 diagnostic and frozen log adversary

The already-known P48 retrospective values are used only for model development. With the full synchronized delete-one-batch covariance, fit

```text
Y_N = N^(5/4) P4[S'] = A + B_logN log N.
```

The five-size GLS result is

```text
A      = -2.422594685734799
B_logN =  1.016646899281392
chi2   =  1.3757723552490904 / 3.
```

On the old development split, fitting only `N=65,85,130` and predicting the already-observed `N=145,170` gives

```text
held-out chi2 = 0.9914259507591816 / 2.
```

For comparison, the original pure power gives `10.190811422597028/2`, while the already-frozen ordinary `q=2` correction gives `1.1249401555843164/2` on the same retrospective split.

The tiny difference between the two two-parameter corrections is **not** evidence for the logarithm. Over this narrow size range, correction shapes are strongly confounded. The purpose of this note is instead to show that the log has an independent LCFT mechanism and should be frozen now, not invented after the new targets.

Frozen source-fit predictions are

```text
N=185:
  N^(5/4) P4[S'] = 2.88466387698 +/- 0.23442255741
  P4[S']         = 0.004227956941 +/- 0.000343585430

N=265:
  N^(5/4) P4[S'] = 3.25002034068 +/- 0.32328609209
  P4[S']         = 0.003039686094 +/- 0.000302363719
```

These uncertainties include source-fit covariance only; target sampling covariance must be added when fresh data are scored.

## 6. Protected scoring order

The existence of a better theoretical motivation for the logarithm does not justify moving goalposts. The chronological scoring order is now:

1. original P48 four pure-power laws, unchanged;
2. zero-effect benchmark;
3. the previously frozen `q=2` ordinary-RG correction to `P4[S']`, no target refit;
4. this newly frozen rank-2 Jordan-log adversary, no target refit;
5. only afterward a free exponent or additional operator sectors.

If `q=2` wins, prefer the ordinary correction explanation for the P48 drift while retaining LCFT structure elsewhere. If the log wins, the next target should be a joint scaling-function test linking `D`, `S'`, higher derivatives and Gaussian-multiplier residuals. If both fail, neither coefficient should be retuned on N=185/265.

## 7. What this does not prove

This note does **not** establish that the lattice matching involution is the required interchiral automorphism, that the particular level-4 generalized state has nonzero lattice coupling, or that its logarithmic coupling is fixed by the primary-level normalization in arXiv:1206.2312. Those are separate operator-identification problems. It establishes only a conditional representation-theoretic route by which the known percolation energy Jordan structure can propagate to the `x=21/4, s=4` sector and generate the exact form now frozen for prospective discrimination.
