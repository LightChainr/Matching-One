# Beyond witness counting: an incremental topological action spectrum

Date: 2026-09-14

Status: correction to the earlier single-fugacity synthesis.  Existing BK bounds and Wulff/network results retain their original status.  The new proposal is to organize rare topology by constrained incremental actions, not by assuming integer powers of one universal fugacity.

## 1. Why the one-fugacity picture was too strong

The earlier shorthand

```text
epsilon_w = exp[-w kappa(p)]
```

is useful for the baseline probability of one horizontal essential connection/component.  It became too strong when it was used to suggest that every additional topology-changing effect costs exactly one more factor of `epsilon_w`.

There are at least three distinct mechanisms:

```text
merger of two earlier essential lineages,
wrap-sensitive correction to cylinder mass,
spectral tunnelling producing the #800 slow doublet.
```

Only the first currently has a direct disjoint-witness probability bound.  The second is an irreducible-decoration perturbation problem; the third is a visible operator/eigenvalue splitting.  They need not share an action.

## 2. Incremental action definition

Let `B_w` be a baseline topological event with

```text
P(B_w)=exp[-w a_B(p)+o(w)].
```

For an additional defect/event `E_w`, define when meaningful

```text
sigma_E(p)
 = liminf_(w->inf)
   -1/w log [ P(E_w intersect B_w) / P(B_w) ].
```

This is the incremental large-deviation cost conditioned on the baseline topology.

For spectral perturbations that are not literally event probabilities, use the analogous exponential rate of the matrix element/eigenvalue splitting and keep it separately typed.

## 3. What #780 really proves

A merger of two already-essential lineages implies two vertex-disjoint essential occupied witnesses in the final configuration.  The one-witness probability has exponent `kappa`; site BK gives an absolute two-witness upper probability with exponent at least `2 kappa`.

Since the baseline component intensity has exponent `kappa`, the rigorous consequence is an incremental lower bound

```text
sigma_merger >= kappa
```

in the exponential-rate sense used by the proof.

Nothing in the argument proves

```text
sigma_merger = kappa.
```

The merger can be strictly rarer because of connectivity/anchor constraints or because the optimal two-witness geometry costs more than two independent cheapest crossings.

## 4. #760 and #800 require distinct actions

### Periodic mass locality

Pure periodization preserves the zero Fourier mode, so `gamma_w-kappa` comes from an irreducible piece that sees a periodic image.  Define its leading action

```text
sigma_wrap(p).
```

Then generically one expects

```text
gamma_w-kappa
 = poly/subexp(w) * exp[-w sigma_wrap(p)]
```

if a unique leading defect exists.  `sigma_wrap=kappa` is one candidate, not the definition.

### Slow doublet

For the projected slow-pole splitting define

```text
sigma_split(p)
 = liminf -w^-1 log Delta gamma_w.
```

The observed `Delta gamma/(w nu)=O(1)` at NN `p=1/4` suggests

```text
sigma_split ~= kappa
```

but this remains a diagnostic until the tagged operator eigenvectors identify the actual tunnelling mechanism.

Opposite residues indicate an approximately antisymmetric source/readout, not by themselves a one-witness probability event.

## 5. Wulff-network actions are the natural replacement

Subcritical 2D connectivity already comes with a directional norm `tau_p`.  A constrained topology is naturally assigned a network cost by minimizing sums of `tau_p` over a graph/skeleton compatible with that topology.

Therefore propose:

```text
sigma_T(p)
 = minimum Wulff/network excess cost for defect topology T
   relative to the baseline winding object.
```

This makes #758's loop/branch variational programme part of the same language: its rate is not generally an integer multiple of `kappa`, because branches can share geometry and the optimum can change topology.

A witness-count argument gives useful inequalities on `sigma_T`; it need not give the exact minimizer.

## 6. Near-critical conjecture: universal action ratios, not integer powers

If rotational restoration holds near criticality so that

```text
tau_p(v)/kappa(p) -> |v|_2,
```

then every fixed finite network topology has a Euclidean variational constant `c_T`.  This motivates

```text
boxed:
sigma_T(p)/kappa(p) -> c_T
```

for the relevant class of defects as `p->pc` from the subcritical side.

Then with

```text
s=w kappa(p),
```

the massive tail is organized schematically by

```text
sum_T exp[-c_T s] P_T(s),
```

not necessarily by integer powers `exp[-j s]`.

This is a stronger but safer unification: the same mass sets the unit of action, while topology determines the coefficient `c_T`.

## 7. What remains true about `w kappa`

The #780 merger upper bound implies that

```text
w kappa -> infinity
```

is sufficient to suppress macro-merger errors in the fixed-subcritical fragmentation picture.

It does **not** prove that every other defect becomes small with exactly `exp[-w kappa]`, nor that finite `s=w kappa` alone determines the whole process.

At finite near-critical `s`, a split--merge/hard-core process remains a plausible continuum target, but its rates may involve several `c_T` values.

## 8. Minimal high-information tests

1. `#800`: identify the two slow eigenvectors and the matrix element coupling them; only then decide whether the tunnelling skeleton has action `kappa`.
2. `#760`: characterize the minimal wrap-sensitive irreducible piece before fitting an exponential rate.
3. `#780`: no need to simulate deeper fixed-subcritical no-merger; existing theorem already gives suppression.  A future finite-`s` calculation should measure one merger observable and compare it with a declared network action.
4. `#758`: use the variational skeleton to compute candidate `c_T` values which can be compared with spectral/process rates.

## 9. Claim boundary

Retained exact/author-level inputs:

- baseline winding/component exponential mass `kappa`;
- BK double-witness upper bound for merger;
- zero-Fourier periodization identity;
- projected slow-doublet numerical diagnostics;
- existing directional/Wulff network definitions.

New conjecture:

- each topological/spectral defect has its own incremental action `sigma_T`;
- near criticality, ratios `sigma_T/kappa` approach geometry/topology constants `c_T` when isotropic Wulff scaling applies.

Superseded simplification:

```text
one extra witness => exactly one universal factor exp[-w kappa]
```

is no longer used as a general law.