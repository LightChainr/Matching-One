# Common-label marked Poisson bridge from the existing component-Poisson proof

Date: 2026-09-14

Status: **candidate author-proof / proof-audit draft**.  This note reorganizes estimates already written in `poisson-birth-windows.md` into the parameter-marked statement needed by #780.  It does not claim independent verification.  The point is to identify whether a genuinely new probability mechanism is missing, or whether the common-label bridge is already latent in the existing AGG/BK proof.

## 1. Setup

Fix a compact strictly subcritical interval

```text
I compactly contained in (0,pc(G)).
```

Let `p0 in I`.  On the infinite cylinder `C_w x Z`, use one common family of iid labels

```text
U_v ~ Uniform(0,1),
```

and declare a site occupied at parameter `p` when `U_v<=p`.

Let

```text
nu_w(p)
```

be the once-per-complete-horizontal-winding-component anchor intensity per longitudinal row, with the same lowest-row/tie-break convention as `poisson-birth-windows.md`.

The reused results are:

1. for `p` in compact subcritical intervals,

```text
-log nu_w(p)/w -> kappa(p)
```

uniformly along convergent parameter sequences;

2. with localization height `H=w^2`, the intensity of components whose vertical span exceeds a generic `r>=w` obeys

```text
rho_w^long(p;r)
 := E[# {complete winding C: anchor_y(C)=0, L(C)>r}]
 <= C_I w exp(-c_I r);
```

3. localized anchor indicators have dependency range `O(wH)` and, for two distinct full winding components in overlapping windows,

```text
E[I_i I_j] <= U_w(p)^2,
U_w(p)=poly(w,H) exp[-(w-1) kappa(p)],
```

by two disjoint occupied winding witnesses and site BK;

4. disjoint localization windows are functions of disjoint site-label sets, hence the AGG `b3` term is exactly zero.

No Ornstein--Zernike prefactor, `p`-analyticity, or affine parameter clock is assumed below.

## 2. Merger factorial intensity

For `p_-<p_+` in `I`, let `C` range over final `p_+` complete winding components.  Let

```text
n_C(p_-)
```

be the number of distinct `p_-` complete winding components contained in `C`.  Monotonicity of occupied components makes this containment unique.

Define the per-row merger factorial intensity

```text
M_w(p_-,p_+)
 = E sum_{C: anchor_y(C)=0} binom(n_C(p_-),2).
```

The intended conclusion is

```text
limsup_w (1/w) log M_w(p_-,p_+) <= -2 kappa(p_-)
```

up to the usual `o(w)` uniformity for moving parameters, and consequently

```text
M_w(p_-,p_+) / nu_w(p_-) -> 0
```

whenever `p_-,p_+ -> p0` in a compact subcritical interval.  For fixed `p_-<p_+`, the ratio conclusion follows whenever the displayed exponential comparison is used with `nu_w(p_-)=exp[-kappa(p_-)w+o(w)]`; the final-component long-span term is superexponential and does not affect it.

### 2.1 Short final components

Put `H=w^2`.  Suppose a final `p_+` component `C` has

```text
L(C)<=H
```

and contains two distinct `p_-` winding components `C1,C2`.

Since `C1,C2 subset C`, each has span at most `H`, and their anchor rows differ by at most `H`.  Therefore both are counted by the localized `p_-` anchor indicators used in the original AGG proof, and their localization windows overlap.

Each unordered merger pair is thus a member of the larger set

```text
{unordered pairs of distinct localized p_- winding components
 whose anchor windows overlap}.
```

The original BK pair estimate gives

```text
E[I_i^- I_j^-] <= U_w(p_-)^2.
```

There are `w` anchor indices per longitudinal row and at most

```text
D_w = w(2H+3)
```

dependency neighbours per index.  By stationarity / mass transport from the merger pair to (say) its lexicographically first early anchor,

```text
M_w^short
 <= (1/2) w D_w U_w(p_-)^2
 <= poly(w) exp[-2 kappa(p_-) w+o(w)].
```

The important point is conceptual: this is exactly the old AGG `b2` scale divided by observation length.  No new two-lineage probability estimate is needed.

### 2.2 Long final components

If `L(C)>H`, then trivially

```text
n_C(p_-) <= |C| <= w L(C),
```

so

```text
M_w^long
 <= (w^2/2)
    E sum_{C:anchor_y(C)=0} L(C)^2 1{L(C)>H}.
```

Use the discrete tail identity

```text
L^2 1{L>H}
 <= H^2 1{L>H}
    + sum_{r>=H} (2r+1) 1{L>r}.
```

The uniform component-span intensity bound then gives

```text
M_w^long
 <= poly(w,H) exp(-c_I H)
 = exp[-Theta(w^2)].
```

Thus the final parameter `p_+` enters only through a superexponentially negligible localization failure in this estimate; the two-witness exponential action is set by the earlier barrier parameter `p_-`.

## 3. Window-level no-merger consequence

Let an observation window have `m` longitudinal rows.  The expected number of merger pairs attached to final components whose anchors lie in the window is

```text
m M_w.
```

Hence Markov gives

```text
P(at least one final component in the window
  contains two earlier winding lineages)
 <= m M_w.
```

For the natural intensity scale

```text
m = O(1/nu_w(p0)),
```

and `p_-,p_+ -> p0`, the estimate above gives

```text
m M_w
 <= exp[-kappa(p0) w+o(w)] ->0.
```

If the observable is defined by components intersecting, rather than anchored inside, the window, enlarge by `H` rows at each endpoint.  Since

```text
H nu_w(p0) ->0,
```

and long components have superexponential intensity, the same conclusion holds.

This is the precise scale-separation statement needed by #780: on a macroscopic `1/nu0` longitudinal window, distinct earlier essential lineages almost surely do not coalesce into one final macroscopic barrier.

## 4. Local birth marks

Now take a bounded parameter window represented by monotone functions

```text
p_w(x),  x in [a,b],
p_w(x) -> p0 uniformly,
```

and write

```text
nu0 = nu_w(p0),
Lambda_w(x)=nu_w(p_w(x))/nu0.
```

Assume only that on the chosen bounded `x` interval

```text
Lambda_w(x) -> Lambda(x)
```

at the finite set of clock values used below; continuity / strict monotonicity of the limiting clock can be imposed later when passing to a continuous mark space.

Fix the top parameter

```text
p_+ = p_w(b).
```

For each localized final `p_+` winding component `C`, define its birth mark

```text
tau(C)
 = inf{p in [p_w(a),p_+]:
       C contains a p-winding ancestor}.
```

On the no-merger event, that ancestor lineage is unique throughout the bounded clock window.  With continuous iid labels, the first threshold is almost surely unique; a direct `rank0 -> rank2` birth inside one component causes no ambiguity for the scalar mark.

If the final component is localized inside its `H`-row window, all lower-`p` ancestors are subsets of that same final component and all guard sites closed at `p_+` remain closed for lower `p`.  Hence `tau(C)` is a measurable function of the same finite label window as the final anchor.

Rare merger configurations may be assigned an arbitrary deterministic tie-break; their total intensity is `o(nu0)` by Section 2, so they disappear in the limiting marked process.

## 5. Finite mark partitions inherit the same AGG proof

Partition `[a,b]` into finitely many mark bins `A_1,...,A_k`.  Define typed localized indicators

```text
I_{j,x,r}
 = 1{a final p_+ component is anchored at (j,x)
     and tau(C) belongs to A_r}.
```

For disjoint localization windows, the entire typed families are independent because they depend on disjoint `U` variables.

For overlapping windows, two distinct typed anchors still imply two distinct final `p_+` winding components.  Their occupied winding witnesses are disjoint, so the same BK estimate gives

```text
E[I_{i,r} I_{j,s}] <= U_w(p_+)^2
```

for every pair of types `r,s`.

Thus the AGG process theorem applies without a new dependency argument.  On an observation length

```text
m = O(1/nu0)
     = exp[kappa(p0)w+o(w)],
```

and because `p_+->p0`,

```text
m poly(w) U_w(p_+)^2 ->0.
```

Therefore the vector of counts in finitely many disjoint spatial arcs and mark bins is asymptotically a vector of independent Poisson variables with its actual finite-`w` means.

## 6. The mark mean measure is determined by component intensity

The remaining task is to identify those means.  For an intermediate parameter `p<=p_+`, consider all final `p_+` components and the number `n_C(p)` of `p` winding ancestors they contain.

By the translation-covariant containment map from each early component to its unique final component, the mass-transport principle gives the exact intensity identity

```text
nu_w(p)
 = E sum_{C: final anchor_y(C)=0} n_C(p).
```

Define the cumulative final-component birth intensity

```text
mu_w((-,p])
 = E sum_{C: final anchor_y(C)=0} 1{n_C(p)>=1}.
```

For every nonnegative integer `n`,

```text
0 <= n-1{n>=1} <= binom(n,2).
```

Hence

```text
0 <= nu_w(p)-mu_w((-,p])
   <= M_w(p,p_+).
```

Uniformly on a bounded clock window with `p_w(x)->p0`, Section 2 therefore gives

```text
mu_w((-,p_w(x)]) / nu0
 - Lambda_w(x)
 ->0.
```

Consequently, if `Lambda_w(x)->Lambda(x)`, then for every finite collection of continuity points

```text
mu_w((p_w(x1),p_w(x2)]) / nu0
 -> Lambda(x2)-Lambda(x1).
```

No derivative of `nu_w(p)` is required.  The intensity itself is the clock.

## 7. Spatial anchor motion disappears on the intensity scale

A localized final component and any of its ancestors lie in the same `H`-row window.  Therefore any two translation-covariant anchor choices attached to that lineage differ by at most `O(H)` rows on the good event.

The rescaled longitudinal coordinate is

```text
nu0 y.
```

Since

```text
H nu0
 = w^2 exp[-kappa(p0)w+o(w)]
 ->0,
```

first-birth anchors, final anchors, lowest-row anchors, or any other bounded-window covariant representative have the same limiting spatial coordinate.  The bad long-component event is already superexponentially negligible.

## 8. Candidate marked-PRM theorem

The preceding lemmas suggest the following theorem, conditional only on routine finite-partition/tightness bookkeeping and the audit of Sections 2 and 6.

> **Candidate theorem.**  On every bounded intensity-clock interval for fixed strictly subcritical `p0`, the marked complete-barrier birth process
>
> ```text
> Xi_w
> = sum_C delta_(nu0 y_C, Lambda_w(tau(C)))
> ```
>
> converges on bounded spatial windows to a unit-rate Poisson random measure
>
> ```text
> dy dLambda.
> ```
>
> Equivalently, for finite disjoint spatial/clock rectangles, the counts converge jointly to independent Poisson variables with means equal to rectangle areas.

A standard route from the already proved finite mark partitions to PRM convergence is to use a countable generating algebra of spatial/clock rectangles, convergence of means, and local finiteness.  No new percolation estimate appears at that stage.

## 9. Pure splitting is then a corollary, not an additional model assumption

At clock `Lambda`, retain all birth points with mark at most `Lambda`.  A Poisson random measure `dy dLambda` gives a homogeneous PPP of barriers of spatial rate `Lambda`, coupled monotonically in the clock by adding new points.

Because macroscopic mergers vanish, an old barrier lineage does not disappear into another old lineage on the scale of interest.  Therefore the complementary white intervals evolve by pure cuts in the limiting process.

All previously derived abstract consequences then become legitimate consequences of the SITE mapping:

```text
fixed-location one-sided gap semigroup,
two-sided Gamma(2,1) fixed-location gap,
Corr(S_x,S_{x+h})=exp(-|h|) in log-intensity time,
Laguerre hierarchy,
record/cut genealogy.
```

They are not independent evidence for the map; the marked-PRM theorem is the bridge that licenses them.

## 10. Audit checklist / possible failure points

Before promoting the candidate theorem, check only the following concrete points.

1. **Generic span tail.**  Confirm that the argument behind `poisson-birth-windows.md` (4.5) indeed yields the unnormalised tail-intensity estimate uniformly for every `r>=w`, not only after substituting `H=w^2`.
2. **Pair-to-final mass transport.**  Write the stationary mass transport from an unordered early merger pair to its containing final component, so the per-row inequality in Section 2.1 has no hidden anchor-offset factor.
3. **Containment semantics.**  Verify that every `p` complete winding component is contained in exactly one `p_+` component and that final localization implies ancestor localization in the same guard window.
4. **Uniformity in the clock window.**  State compact-`p` constants once and use `sup_x |p_w(x)-p0|->0` to make `M_w/nu0->0` uniform.
5. **Typed AGG statement.**  Quote the process/multitype contraction explicitly; no new bound is expected.
6. **Finite partitions to PRM.**  Supply the standard point-process tightness argument.
7. **Finite torus transfer.**  If the theorem is stated on finite exponentially long tori, reuse the original `Bad_H` coupling and endpoint enlargement rather than rebuilding the proof.

If one of these fails, preserve the exact failed arrow.  In particular, failure of a bookkeeping step is not evidence for a physical coalescent; a surviving merger mechanism would require failure of the two-witness/localization scale separation itself.

## 11. Research consequence

If the audit passes, #780 should no longer spend primary compute on estimating merger rates.  The scientific status changes from

```text
abstract pure-cut kernel + missing SITE map
```

to

```text
SITE component Poisson proof
 -> no-macro-merger by its own b2 scale
 -> local birth marks
 -> marked Poisson rain in intensity time
 -> pure splitting / record kernel.
```

This would be exactly the kind of cross-model/scale connection demanded by `research-compass-beyond-exactness-20260914.md`: a previously isolated solvable process becomes the consequence of a controlled Bernoulli-site map rather than a growing independent exact-model programme.
