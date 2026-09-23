# Intensity-clock separation for the common-label gap process

Date: 2026-09-14

Status: exact Poisson-process algebra plus a conditional interface for #780. This note does **not** re-claim the fixed-width record-process generator already written on the #771/#764 lineage. Its purpose is to isolate the unmarked semigroup that a future large-width common-label theorem would inherit once the SITE barrier-birth clock is identified.

## 1. Abstract Poisson clock

Let `Pi(dy dLambda)` be a unit-rate Poisson random measure on longitudinal space times cumulative barrier intensity. At clock value `Lambda`, the active barriers form a homogeneous PPP of rate `Lambda` on the line.

Fix spatial position 0 and let `X^-_Lambda, X^+_Lambda` be the distances to the nearest barriers on the two sides. Define

```text
U^-_Lambda = Lambda X^-_Lambda,
U^+_Lambda = Lambda X^+_Lambda.
```

At every fixed clock value, `U^-` and `U^+` are independent `Exp(1)`, and the fixed-location gap

```text
S_Lambda = U^-_Lambda + U^+_Lambda
```

is `Gamma(2,1)`. This is the uniform-location protocol, not the once-per-component Palm spacing law.

## 2. Exact two-time kernel

Let `Lambda_2 = r Lambda_1`, `r>=1`. On one side, conditional on `U_1=u`, new barriers form an independent rate `r-1` process in units of `Lambda_1`. If `Y~Exp(r-1)`, then

```text
U_2 = r min(U_1,Y).                                      (2.1)
```

Equivalently, there is an atom at `v=ru` of mass `exp[-(r-1)u]`, and on `0<v<ru` the transition density is

```text
k_r(u,v) = (r-1)/r * exp[-(r-1)v/r].                    (2.2)
```

Integrating against the stationary `Exp(1)` law gives the exact one-side joint Laplace transform

```text
J_r(s,t)
 = E exp[-s U_1 - t U_2]
 = (r+s) / ((1+s) [s+r(1+t)]).                         (2.3)
```

The two sides are independent, hence

```text
E exp[-s S_1 - t S_2] = J_r(s,t)^2.                    (2.4)
```

In particular

```text
Cov(S_1,S_2)=2/r,
Corr(S_1,S_2)=1/r.                                      (2.5)
```

Writing `x=log Lambda`, this becomes `Corr(S_x,S_{x+h})=exp(-|h|)`.

## 3. Log-clock generator

For one side, in an infinitesimal increment `dx`, deterministic rescaling sends `u` to `e^{dx}u`, while a new nearest record occurs at rate `u` and replaces the physical distance by a uniform fraction of its previous value. The stationary generator is therefore

```text
L f(u)
 = u f'(u)
   + u integral_0^1 [f(uv)-f(u)] dv.                   (3.1)
```

The two-sided state `(U^-,U^+)` has generator `L_-+L_+` and stationary density `exp(-u-v)`.

The event that neither side is cut between clocks `Lambda_1` and `Lambda_2` has probability

```text
r^-2 = exp[-2(log Lambda_2-log Lambda_1)].              (3.2)
```

This is the unmarked projection of the stronger length/cavity record process already present on the #771/#764 lineage.

## 4. Clock--kernel separation principle

For the actual common-label SITE model, the nontrivial model-dependent object should be the cumulative complete-barrier intensity

```text
Lambda_w(p).
```

A sufficient process theorem would prove, uniformly on bounded intensity-clock intervals,

1. marked complete-barrier births converge to a Poisson random measure in `(y,Lambda)`;
2. a final localized barrier contains two macroscopically distinct earlier essential lineages with probability `o(1)`;
3. localization/anchor motion errors are uniform in the clock window.

Under those three statements, every fixed-location multitime white-gap law follows from (2.1)--(3.2) after the deterministic time change `p -> Lambda_w(p)`. The microscopic lattice, mass normalization, and possible massive field theory enter the clock; the lineage kernel does not need to be re-derived for each model.

This suggests the following separation for #780/#782:

```text
massive / SITE insertion physics  ->  Lambda_w(p)
Poisson birth theorem             ->  space x intensity cloud
universal record geometry         ->  multitime gap genealogy.
```

## 5. Sampling boundary

A once-per-component Palm gap has one-time law `Exp(1)` after rate normalization. The fixed-location interval above has `Gamma(2,1)`. Tracking a component descendant therefore requires a declared rule (anchor descendant, uniform-position descendant, all children, or size-biased child). The semigroup in this note is the fixed-location protocol and must not be silently relabelled component-Palm.

## 6. Claim boundary

- Exact: Sections 1--3 for the abstract Poisson cloud.
- Conditional: Section 4 as an interface theorem for the actual large-width common-label SITE filtration.
- Not claimed: uniform-in-parameter AGG for the SITE anchors, absence of mergers, or a square-site near-critical clock formula.
