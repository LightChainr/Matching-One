# A bulk-invisible rank projection reverses the strong-coupling U tail

## Two fixed laws give opposite signs

Keep the already defined source `Sstar=C+F4+Bvac`. Compare it with the
existing counterfactual `Sdrop=Sstar+r`, obtained by dropping only the
local-colour representation's topological factor `m^(-r)`. There is no new
fitted source coefficient. Put lambda=exp(-t), N=L², A_N=N^(13/8)/2 and
Delta=cos4(theta_axis)-cos4(theta_companion)>0. On the axis L×L torus,
L>=5, paired with a same-area honest quotient with ell1>=L+2,

```text
U_star/A_N = -(L²-6L+6)/Delta * lambda^(2L+1) + O(lambda^(2L+3)),
U_drop/A_N = +(L-2)/Delta * lambda^(2L-2+2/L)
            + O(lambda^(2L-2+4/L)).
```

The first law is execution's [square-family theorem at762dbaf4](https://github.com/LightChainr/Matching-One/blob/762dbaf4c3afd9925f7e39b27220274312db4dc4/notes/closed-source-square-family-leading-law.md).
The second is derived here from the same classified winding layer. For each
fixed L, one law approaches zero from below and the other from above. This
excludes treating removal of the rank projection as a harmless thermal
coordinate change for the original global U, despite its vanishing pressure
density effect at fixed coupling and growing N.

## Why deleting the projection breaks the minimal-stripe cancellation

The [existing colour-law comparison](https://github.com/LightChainr/Matching-One/blob/85fd492312b597b3fa102ea913e4bcc7aeae2acf/notes/closed-source-local-colour-gas.md)
fixes the unprojected coordinate and defect action:

```text
d = [p/((1-p)m)] m^(2/N),
weight_drop = d^K lambda^eta,
eta = g-r+2K/N,   g=Bmix-2C_B+r.
```

The empty and full weights are1 and d^N. Every other configuration has
eta>=2+2/N. Both geometries are normalized separately before pooling.
The unique simple root therefore has
`d0=1+O(lambda^(2+2/N))` and `Q_d(d0)=N/2+O(lambda^(2+2/N))`.
The common thermal-coordinate Jacobian cancels in U; this is not a new
observer normalization.

Execution's boundary classification gives the complete minimal rank-one
layer on the axis: `g=2L-1`, exactly2L straight stripes for each width
w=1,...,L-1, with K=Lw. On this layer

```text
eta_w=2L-2+2w/L,
Z1_drop,min=2L sum_(w=1..L-1) d^(Lw) lambda^eta_w.
```

Under Sstar all these widths had the same source cost, allowing their
normalized reciprocal polynomial to cancel its thermal derivative at the
balanced root. Under Sdrop different widths have different eta. The unique
lowest-cost width is1; its leading normalized rank-one probability is
`2L d^L/(1+d^N) * lambda^(2L-2+2/L)`.

At d=1, differentiating *after normalization* gives

```text
[P1_d]_leading = L(L-N/2),
[E_d]_leading = L(N/2-L),    E=1-P1.
```

Divide by Q_d=N/2 and by Delta. This gives the coefficient `(L-2)/Delta>0`.
Root motion enters strictly later and cannot restore the original reciprocal
shell cancellation.

## No omitted sector or shell can precede this term

Every rank-one occupied graph contains a noncontractible simple cycle,
so K>=ell1. On the axis K>=L. Higher axis rank-one layers have
g>=2L+1, because g is odd in rank1; hence their eta is at least2L+2/L,
strictly beyond the leading stripe and the next width. The companion has
g>=2ell1-1>=2L+3 and K>=ell1>=L+2, so it also starts later.
The next minimal-axis width raises eta by2/L, while normalization and
root corrections raise the first term by at least2+2/N. This proves the
displayed remainder. These are finite polynomial/Puiseux expansions in
xi=lambda^(1/N), not fitted powers or uncontrolled differentiation of an
unspecified big-O term.

The same argument gives the complete minimal-strip band

```text
U_drop/A_N = sum_(w=1..L-1) (L-2w)/Delta
              * lambda^(2L-2+2w/L) + O(lambda^(2L+2/L)).
```

Only a safe remainder bound is used here; the first-term comparison requires
no higher-shell classification or additional coefficients.

## Explicit size and relative-response predictions

For the dilated Gaussian family (5k,0)/(4k,3k), Delta=1152/625:

| N | L | Projected leading U/A_N | Projection-deleted leading U/A_N |
|---:|---:|---|---|
|25|5|`-lambda^11/Delta`|`+3 lambda^(42/5)/Delta`|
|100|10|`-46 lambda^21/Delta`|`+8 lambda^(91/5)/Delta`|
|225|15|`-141 lambda^31/Delta`|`+13 lambda^(422/15)/Delta`|

Thus, at each fixed size,

`U_drop/U_star ~ -(L-2)/(L²-6L+6) * exp[(3-2/L)t]`.

The N25 leading support can be read directly from the locked `(K,g,q)`
histogram. N100/N225 here are combinatorial predictions, not new measured
rows. The pressure-density discrepancy between these laws is at most2t/N;
its vanishing at fixed t does not imply equality of finite topological
observers. The two size/coupling limits have not been interchanged.

This settles the strong-coupling sign of the already named counterfactual.
It does not supply a finite-t crossover location, accessible sampling budget,
thermodynamic transition or continuum H4 field. No arbitrary rank-fugacity
scan, added coupling point or rescue of the completed P154/P334/F4 decisions
is attached to this result. The next open comparison is finite-coupling
transport of these two fixed laws, with the pressure and observer criteria
kept distinct.
