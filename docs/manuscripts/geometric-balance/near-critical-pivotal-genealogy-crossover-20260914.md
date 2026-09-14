# Near-critical pivotal genealogy as the crossover completion of the fixed-subcritical splitting process

2026-09-14.

Status: synthesis/conjectural bridge for square-site Matching One, with a rigorous positive-control framework available for triangular-lattice near-critical percolation through the Garban--Pete--Schramm scaling limit. This note does not claim square-site universality has been proved.

## 1. Why the fixed-subcritical theorem should stop where it does

The current #780 route has essentially closed, on compact strictly subcritical parameter intervals, the steps

```text
rare complete winding components -> spatial Poisson process,
common-label birth marks -> marked PPP in intrinsic intensity clock,
merger factorial density / component intensity -> 0,
pure splitting genealogy on bounded clock windows.
```

The no-merger estimate pays one additional winding action. Schematically,

```text
M_w / nu_w <= poly(w) exp[-kappa(p) w].
```

At fixed subcritical p this vanishes exponentially. But the estimate itself announces its failure variable:

```text
x = w kappa(p).
```

When `x=O(1)`, there is no exponential reason for mergers to disappear.

## 2. The near-critical window is exactly x=O(1)

For two-dimensional percolation,

```text
xi(p) ~ |p-pc|^-4/3,
kappa(p) ~ xi^-1 ~ |p-pc|^(4/3).
```

The usual near-critical scaling

```text
p = pc + lambda w^-3/4
```

therefore gives

```text
w kappa(p) ~ const * |lambda|^(4/3)=O(1).
```

So the failure of the fixed-subcritical merger suppression is not an accidental technical gap. It occurs at exactly the universal near-critical scaling window.

## 3. Candidate continuum object: pivotal-driven essential-component genealogy

Garban--Pete--Schramm construct the near-critical scaling limit by perturbing the critical configuration with a Poissonian pivotal measure. In an O(1) lambda interval, O(1) macroscopic pivotal updates occur in a bounded continuum region. These updates can alter macroscopic connectivity, hence can both create and merge essential lineages.

This suggests the correct continuation of #780 is not a corrected Poisson barrier process but a continuum process

```text
G_lambda(tau)
```

of homology-carrying clusters/components on a torus or cylinder of fixed continuum modulus `tau`, driven by the near-critical pivotal measure.

The observable process should retain at least

```text
ambient homology rank / primitive class,
essential component partition,
first-birth marks,
merger/split genealogy,
optional neutral component count.
```

The rank process alone is a lossy projection but has a clean one-time generating function.

## 4. Scaling prediction for merger density

Draft #773/#782 uses the near-critical component-intensity ansatz

```text
nu_w(pc+lambda/(a_t w^(3/4)))
 = w^-1 I(lambda)[1+o(1)].
```

If macroscopic pivotal updates have a nondegenerate continuum rate, then the per-row factorial merger intensity should have the same dimensional scale:

```text
M_w(lambda_-,lambda_+)
 = w^-1 J(lambda_-,lambda_+)[1+o(1)],
```

rather than `o(nu_w)`.

Therefore

```text
M_w/nu_w -> J/I,
```

with a generally nonzero universal crossover function.

This is the sharp qualitative distinction from the fixed-subcritical theorem.

### IR matching condition

As lambda tends deep into the subcritical side, equivalently `x=w kappa -> infinity`, the continuum crossover must reproduce the dilute theorem:

```text
J/I -> 0.
```

The existing BK/AGG proof suggests an IR tail at least qualitatively of the form

```text
J/I ~ polynomial(x) * exp(-x)
```

up to the massive periodic-sewing prefactor. The exact polynomial and amplitude are not supplied by the current proof.

## 5. Interface to the rank-source programme

Define the near-critical rank-source marginal

```text
Z_tau(lambda,s)
 = E[ exp(s(r_lambda-1)) ]
 = Pi0(lambda;tau)e^-s + Pi1(lambda;tau) + Pi2(lambda;tau)e^s.
```

This is precisely the one-time rank projection of the proposed pivotal genealogy.

Hence #782's massive-Potts programme and the near-critical probabilistic route should be viewed as complementary:

- pivotal/CLE route: constructs the process and its merger/split semantics;
- massive Potts / modified trace: may compute one-time sector weights and IR amplitudes;
- Krushkal / affine-TL source dictionary: supplies exact finite topological/neutral fugacities.

One should not wait for a full massive TBA formula before defining the continuum genealogy.

## 6. A useful two-time observable

The natural continuation of the fixed-subcritical merger statistic is

```text
J_tau(lambda1,lambda2)
 = E sum_{C at lambda2} binom(n_C(lambda1),2),
```

with the torus/cylinder normalized to continuum size one.

At fixed aspect ratio this is a finite continuum observable of the pivotal process. On the lattice, the corresponding quantity should satisfy

```text
w M_w(lambda1,lambda2) -> J_tau(lambda1,lambda2).
```

This is more informative than asking only whether a merger was observed.

A second high-value observable is the conditional transition matrix of the aggregate homology rank

```text
P(r_lambda2=j | r_lambda1=i), i,j in {0,1,2},
```

under the common-label monotone coupling. It is not Markov-complete, but it gives a direct bridge to the existing rank-source and torus-homology assets.

## 7. Strong conjecture: fixed-subcritical pure splitting is an IR boundary condition

The proposed universal picture is

```text
near-critical pivotal genealogy at finite lambda
        |
        | lambda -> -infinity / x -> infinity
        v
Poisson barrier cloud + pure splitting clock kernel.
```

In the reverse direction, the fixed-subcritical process does not need to be manually patched with rare mergers. Its correct crossover completion is the pivotal process itself.

This also predicts that the Laguerre hierarchy and exponential gap correlations from #785 are not globally valid in lambda. They should emerge only in the dilute IR after the nonlinear change from lambda to component-intensity clock.

## 8. Square-site scope

The rigorous near-critical scaling-limit construction is currently a triangular-site positive control, not a theorem for square-site Matching One. For the present project the square-site claim is a universality conjecture.

That is still scientifically useful because it gives a precise target and a control order:

1. define/compute the rank and merger observables in the rigorous triangular near-critical coupling;
2. verify the deep-subcritical IR tends to the Poisson splitting law;
3. only then compare square-site finite-width data in the same normalized lambda/x variables.

## 9. Concrete next work

No new fixed-subcritical common-label simulation is needed for this bridge. Higher-value tasks are:

1. derive the continuum measurability of torus homology rank and the merger factorial functional in the near-critical quad-crossing/pivotal process;
2. determine whether existing near-critical CLE/pivotal results imply continuity in lambda for these functionals away from exceptional pivotal times;
3. compute or bound `J_tau(lambda1,lambda2)` in a triangular positive control, analytically if possible and otherwise with a continuum/discrete convergence experiment;
4. connect the one-time marginal `Z_tau(lambda,s)` to the massive-Potts modified-trace programme in #782;
5. identify the IR matching variable and amplitude connecting `J/I` to the fixed-subcritical BK/AGG tail.

The central decision is no longer `MERGER_YES/NO`; it is the crossover function `J/I` and its massive/near-critical normalization.
