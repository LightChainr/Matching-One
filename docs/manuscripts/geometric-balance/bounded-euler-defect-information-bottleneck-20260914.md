# The bounded Euler defect as the root-information bottleneck

Date: 2026-09-14

Status: exact consequence of the finite matching Euler identity, plus research interpretation.  No new continuum field is identified.

## 1. An extensive-looking expression is exactly three-valued

For every honest square-torus site configuration,

```text
X := r4-1
   = k4-k8-K+E-F0
   in {-1,0,+1}.
```

Each term on the right can be extensive:

```text
k4, k8, K, E, F0 = O(L^2)
```

for typical configurations.  Their combination is nevertheless bounded **configuration by configuration**, not only after expectation.

This is much stronger than an average cancellation.

## 2. Why typical-object asymptotics can miss the root

The Matching-One root solves

```text
E_p X = 0.
```

Suppose a coarse/effective theory approximates separately

```text
k4/L^2,
k8/L^2,
K/L^2,
E/L^2,
F0/L^2
```

to small absolute error.  Unless those approximations preserve the exact correlated Euler cancellation, their reconstructed error in `X` can be much larger than the true `O(1)` signal.

Thus very accurate bulk laws for cluster densities, component sizes, boundary fractions or local motifs do not automatically determine the root.

The missing information is precisely the bounded topological Euler defect left after cancellation.

This is a concrete microscopic realization of the research-compass warning that “a solved typical limit can discard the quantity that fixes the root.”

## 3. Charged-sector susceptibility is the variance of the Euler defect

Because

```text
X=-1,0,+1,
```

we have exactly

```text
E X     = P2-P0 = M,
E X^2   = P0+P2 = chi,
Var X   = chi-M^2.
```

At the balance root `M=0`,

```text
boxed: chi = Var(X).
```

So the “charged-sector susceptibility” is literally the variance of the bounded Euler defect.

The topological source generating function is

```text
Z(h)=E e^{hX}=P1+P2 e^h+P0 e^-h.
```

At the root, odd source cumulants vanish and even cumulants are fixed by `chi`; for example

```text
kappa_2(X)=chi,
kappa_4(X)=chi-3chi^2.
```

No separate high-dimensional cluster model is required to know the source law once `(M,chi)` are known.

## 4. The dangerous direction is reconstruction, not direct evaluation

There are two numerically/theoretically different strategies:

### Direct topological evaluation

Compute `r` or `X` from lifted homology / rank classification.  The result is bounded and stable.

### Reconstructed Euler evaluation

Compute the large pieces `k4,k8,K,E,F0` and subtract.  This is exact only if all pieces are obtained on the **same configuration** with mutually consistent conventions.

Independent marginal approximations, separately fitted asymptotics, or source derivatives with missing normalization can destroy the cancellation.

This gives a general warning for the project:

> Whenever a proposed “explanation” of the root passes through several extensive pieces, verify the configurationwise identity before interpreting residuals as new physics.

The #802 source-normalization failure is an example of the same structural hazard: an omitted common normalization turned an exact/gauge cancellation into an apparent large response.

## 5. Local insertion form

Let a white site `v` be switched to black.  Write

```text
c_b = number of distinct black NN neighbour components touched by v,
d_b = number of occupied NN neighbours of v,
f_b = number of elementary faces completed black by inserting v,
t_w = number of nonempty white-matching components produced when v is removed from its old white component.
```

Then

```text
Delta k4 = 1-c_b,
Delta k8 = t_w-1,
Delta K  = 1,
Delta E  = d_b,
Delta F0 = f_b.
```

Hence the exact Euler source gives

```text
boxed:
Delta_v X
 = 1 - c_b - t_w + d_b - f_b.                       (5.1)
```

The left side is the rank insertion jump in `{0,1,2}`.

Equation (5.1) is a local-connectivity certificate for topological birth.  It connects the pivotal language of #768/#769 to the cluster/Euler language without naming an arm field.

In particular, a direct `0->2` birth is not characterized by a bare local degree alone; it is a mismatch between black-component merging, white-component splitting and the local edge/face Euler terms.

## 6. A revised microscopic question for the leading correction

Instead of starting from

```text
“is the L^-4 correction an 8-arm field or a spin-four field?”
```

one can ask a source-defined question:

> which finite-size/angular correction first biases the distribution of the exact Euler defect `X` away from its continuum balanced law?

The answer can then be decomposed empirically into

```text
angular irrep (H0/H4/H8/...),
radial exponent,
microscopic source / pivotal channel,
```

before assigning a continuum operator.

The existing oblique data say the first visible root bias is overwhelmingly H4-like.  #808 asks what survives after that irrep is projected out.  Equation (5.1) supplies a lattice-side target for classifying the corresponding insertion events.

## 7. Consequence for effective models

A component-Palm/Poisson/fragmentation model may correctly describe most geometry while not preserving the Euler defect.  To claim it explains Matching One, it must additionally reproduce at least one of:

```text
E X,
the sourced law Z(h),
or an equivalent sector free-energy difference.
```

Matching only component density, mean span or neutral-count law is insufficient unless a theorem transports the Euler defect through the approximation.

This criterion is stricter than “the effective model fits several observables,” but directly aligned with the root.

## 8. Claim boundary

Exact:

- configurationwise bounded Euler identity;
- `chi=Var X` at the root;
- local insertion formula (5.1).

Interpretive/programmatic:

- use Euler-defect preservation as a criterion for whether an effective model explains the root;
- classify leading angular/radial corrections through the distribution of `X` before CFT field naming.

The main message is that the root lives in an exact correlated cancellation, not in any one of its extensive ingredients.