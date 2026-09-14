# A two-ended branch explanation for the span doublet

Date: 2026-09-14

Status: a mechanism-level conjecture motivated by the visible-pole diagnostic in `span-visible-pole-splitting-20260914.md` and by the loop--branch decomposition in #758. It is more specific than the generic metastable-doublet interpretation and yields direct tests for #762.

## 1. The algebraic signature

The stored complete-component span spectra show two slow positive visible poles

```text
rho_+ > rho_-
```

with nearly opposite scalar residues. This is exactly the signature produced by a convolution of two one-sided geometric/exponential waiting laws.

If two nonnegative branch lengths `A_-`, `A_+` have generating functions

```text
B_-(z) ~= C_-/(1-rho_- z),
B_+(z) ~= C_+/(1-rho_+ z),
```

then their sum has

```text
B_-(z) B_+(z)
 ~= const / [(1-rho_- z)(1-rho_+ z)].                    (1.1)
```

The coefficient is a difference of two exponentials with opposite residues:

```text
[z^h] B_- B_+
 proportional to
 (rho_+^(h+1)-rho_-^(h+1))/(rho_+-rho_-).                 (1.2)
```

As the two one-sided masses coalesce, `rho_+-rho_- -> 0`, (1.2) tends to

```text
const * (h+1) rho^h.                                      (1.3)
```

Thus a near-double pole naturally produces both observations in #800:

- coalescing visible roots;
- residues with ratio near `-1`.

No metastable tunnelling interpretation is required for this algebra.

## 2. Geometric interpretation for complete winding components

In the saturated large-span branch of the #758 candidate, write schematically

```text
L = L_core + A_- + A_+,
```

where `L_core` is the vertically saturated winding core and `A_-`, `A_+` are lower/upper longitudinal decorations or branches.

The rate candidate above saturation pays a linear cost only for the **total** excess branch length. This leaves a one-dimensional split degeneracy between the two ends. If the two far branch excursions are asymptotically separated by the saturated core and share the same one-sided cylinder mass, then the natural first approximation is precisely the convolution in Section 1.

This gives a direct bridge:

```text
#758 linear branch rate
        +
two-ended branch split
        ->
double-pole / Erlang-type subexponential factor.          (2.1)
```

The LDP exponent can therefore be correct even when the fixed-w scalar response never looks like a single exponential on height `O(w)`.

## 3. New subexponential prediction

In the exactly coalesced idealization,

```text
d_h ~= C h e^{-gamma h}.                                  (3.1)
```

Consequently the tail also has a linear polynomial factor,

```text
P(L>=h) ~= C' (h+c0) e^{-gamma h},                        (3.2)
```

and the finite-height hazard satisfies

```text
gamma_eff(h)
 = -log[T(h+1)/T(h)]
 = gamma - 1/h + O(h^-2)                                  (3.3)
```

in the pre-splitting regime `h |gamma_+-gamma_-| << 1`.

For `h=A w`, the `log h` prefactor contributes only `O(log w/w)` to the finite-size rate and therefore does not alter the #758 linear LDP slope. It does, however, matter for prefactor diagnostics and for any attempt to identify the ultimate fixed-w pole from moderate heights.

The current NN `p=1/4` hazards are qualitatively consistent with this correction; this note does not claim that the coefficient `1` in (3.3) has already been numerically certified for the full SITE model.

## 4. A direct #762 prediction: uniform branch split

If in the asymptotic branch regime the two one-sided excesses are independent exponentials with the same rate,

```text
A_- ~ Exp(gamma),
A_+ ~ Exp(gamma),
```

then conditional on their total

```text
R=A_-+A_+,
```

the fraction

```text
U=A_-/(A_-+A_+)
```

is exactly

```text
boxed: U | R  ~ Uniform(0,1).                              (4.1)
```

Equivalently `(A_-/R,A_+/R)` is `Dirichlet(1,1)`.

This turns one of #762's morphology options into a sharply motivated test rather than a generic candidate. In the saturated branch regime, the two-ended convolution mechanism predicts simultaneously:

1. `L_core/w` saturates near the #758 optimizer;
2. the excess `L-L_core` carries the linear tail cost;
3. `U` approaches Uniform(0,1) away from zero-excess cases;
4. the scalar span spectrum develops a double-pole/Erlang prefactor.

Failure of (4.1), especially strong endpoint condensation, would favor an asymmetric single-long-branch mechanism and would require a different explanation of the opposite residues.

## 5. Split poles as a finite-width deformation

The empirical poles are close but not exactly equal. A natural finite-width deformation is

```text
gamma_- = gamma_bar - Delta/2,
gamma_+ = gamma_bar + Delta/2.                            (5.1)
```

Then the branch-sum coefficient is the hypoexponential form (1.2). The near-equal masses may come from weak endpoint/core asymmetry or from interaction between the two branch excursions through the finite core.

This interpretation is distinct from, but not mutually exclusive with, a metastable two-sector transfer block. The decisive check is the actual slow eigenvector geometry:

- if the two modes localize on upper/lower branch or two serial excursion stages, the convolution mechanism is supported;
- if they instead distinguish a global topology bit unrelated to branch orientation, the metastable-sector mechanism is more plausible.

## 6. Relation to the loop--branch rate

Above the saturated core size `r_*(p)`, the #758 candidate is

```text
I_p(A)=kappa A + c_*,
```

so all partitions of the excess span between the two ends have the same leading exponential cost. The continuum of branch splits is exactly the kind of zero mode that produces a polynomial prefactor while leaving the large-deviation rate unchanged.

A sharpened conjecture is therefore:

```text
P_Palm(L≈A w)
 = w^beta C_p(A) exp[-w I_p(A)] [1+o(1)],                (6.1)
```

with an additional `beta=1` contribution from the two-ended split degeneracy in the strictly saturated linear branch, relative to a convention where the core location and one branch split are otherwise fixed. The total exponent beta also contains rooting/unrooting and transverse fluctuation factors, so `beta=1` is **not** claimed as the final complete-component prefactor.

The robust part is the predicted linear-in-excess split measure and Uniform(0,1) conditional fraction, not a final absolute power of w.

## 7. What to test next

### Existing operator, no new width production

Extract the two slow visible eigenvectors and inspect whether their state mass distinguishes upper/lower branch stages or a global topology bit.

### #762 morphology

At `L>=A w` for `A=1,2`, record the already planned

```text
L_core,
A_-, A_+,
U=A_-/(A_-+A_+).
```

The strongest mechanism test is not the mean of U but its conditional law at fixed excess-span bins.

### Spectrum

Fit the late sequence directly to both models:

```text
M1: c1 rho1^h + c2 rho2^h,
M2: (a+b h) rho^h,
```

on nested windows. If `M2` becomes competitive as w grows while the fitted pole split shrinks, that is the expected coalescing-convolution signature.

## 8. Claim boundary

- Exact algebra: convolution of two simple poles gives opposite residues and a double-pole/Erlang limit; equal-rate exponentials give the Uniform branch fraction.
- Data-supported hypothesis: the SITE span doublet may be realizing this two-ended branch convolution.
- Not claimed: asymptotic independence of `A_-` and `A_+`, identification of the slow eigenvectors, or a final complete-component w-prefactor.
