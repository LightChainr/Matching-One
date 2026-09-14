# In the fixed-subcritical dilute regime, void free energy equals winding-component intensity to first order

2026-09-14.  Author-level rare-event pressure argument joining the component-Poisson and fixed-width charge-free-energy descriptions.

The important distinction is:

- at finite width, the no-winding strip free energy `I^0_w` and the complete winding-component intensity `nu_w` are not identical;
- for fixed subcritical `p` and `w->infinity`, winding components become an exponentially dilute locally dependent gas, and the pressure/activity distinction disappears at first order.

The result claimed here is

\[
\boxed{
\frac{I^0_{G,w}(p)}{\nu^G_w(p)}\longrightarrow1,
\qquad w\to\infty,
\quad 0<p<p_c(G)\text{ fixed}.}                              \tag{1}
\]

It applies separately to `G4` or `G8` at a fixed subcritical parameter for that graph.  It is **not** a near-critical statement with `p=p_w->p_c`.

## 1. Complete-component activity

On the infinite cylinder `C_w x Z`, every occupied component is vertically finite almost surely.  Count each horizontally essential component exactly once by the bottom-row/tie-break anchor from `poisson-birth-windows.md`.  Its row intensity is

\[
\nu_w=
E[\#\{\text{complete winding-component anchors in row }0\}]. \tag{1.1}
\]

The existing fixed-subcritical argument gives

\[
-\frac1w\log\nu_w\to\kappa_G(p)>0.                           \tag{1.2}
\]

## 2. Linear-height localization is enough for pressure

The birth-window proof used `H=w^2` to make localization errors superexponentially small.  For the present first-order pressure statement choose instead

\[
H=Cw,                                                         \tag{2.1}
\]

where `C` is a sufficiently large fixed constant depending on the declared subcritical `p`.

Uniform cylinder height tails give

\[
0\le\nu_w-\nu_{w,H}
\le C_1w e^{-c_1H}.                                          \tag{2.2}
\]

Choose `C` so that

\[
c_1C>\kappa_G(p)+2\eta                                      \tag{2.3}
\]

for some `eta>0`.  Then by (1.2),

\[
\boxed{\nu_w-\nu_{w,H}=o(\nu_w).}                            \tag{2.4}
\]

Thus it is enough to compute the void pressure of the finite-window anchor field.

## 3. Dependency graph and k-fold witness bound

Let `I_i`, `i=(j,x)`, be the localized anchor indicators with height cutoff `H`.  Each is measurable in `O(H)` rows and one circumference, so the dependency graph has maximum degree

\[
D_w=O(wH)=O(w^2).                                             \tag{3.1}
\]

Disjoint windows are independent.

Let `B_w` be the first-span upper bound for the event that a window union of `O(H)` rows contains a horizontal winding witness.  At fixed subcritical `p`,

\[
B_w\le\operatorname{poly}(w)e^{-\kappa_G(p)w}.               \tag{3.2}
\]

The nonmonotone anchors themselves are NOT fed into BK.  Instead, if `k` distinct anchors occur, they belong to `k` distinct complete components.  Choose one occupied horizontal-winding witness inside each component.  These occupied witness site sets are disjoint.  Hence, for every distinct tuple whose relevant windows lie in one connected dependency cluster,

\[
\{I_{i_1}=\cdots=I_{i_k}=1\}
\subset E_1\square\cdots\square E_k,                         \tag{3.3}
\]

where the `E_a` are increasing enclosing winding events.  Repeated SITE BK gives

\[
\boxed{
E\prod_{a=1}^k I_{i_a}
\le B_w^k.}                                                   \tag{3.4}
\]

This is stronger than the pair bound used for Chen--Stein.

## 4. A local-dependence pressure lemma

Consider a stationary family of zero-one variables on `m` rows, with finitely many types per row and a dependency graph of degree at most `D`.  Suppose joint moments of every distinct connected `k`-set obey

\[
E\prod_{i\in S}I_i\le B^{|S|}.                               \tag{4.1}
\]

Let

\[
Z_m=E\prod_{i\in\Lambda_m}(1-I_i)
=P(\text{no marked event in the length-}m\text{ window}).     \tag{4.2}
\]

Because joint moments factor across disconnected dependency-graph components, the logarithm of (4.2) has the standard connected-cluster expansion.  The first-order contribution is

\[
-\sum_i EI_i.                                                 \tag{4.3}
\]

Every higher term is supported on a connected index set.  A fixed root belongs to at most

\[
(eD)^{k-1}                                                     \tag{4.4}
\]

connected ordered/tree-coded `k`-clusters up to an inessential universal factor.  The usual tree-graph/cumulant bound then gives, whenever `cDB<1`,

\[
\left|-rac1m\log Z_m
-\frac1m\sum_iEI_i\right|
\le
C\,n_{type}\,D B^2
\sum_{r\ge0}(cDB)^r.                                          \tag{4.5}
\]

Here `n_type` is the number of anchor types per row (`w` in the raw site-anchor convention).  The constants are universal and irrelevant for the exponential comparison.

Equation (4.5) is the standard small-activity connected-cluster fact behind dependency-graph/LLL polymer expansions: **log void probability contains only connected event clusters.**  Scott--Sokal / Dobrushin cluster-expansion criteria provide an abstract framework; the special moment bound (4.1) makes the present estimate elementary by tree counting.

For the winding anchor field,

\[
DB_w=\operatorname{poly}(w)e^{-\kappa w}\to0,                \tag{4.6}
\]

so the expansion is absolutely convergent for all large `w`.

## 5. Apply the pressure lemma

Per row the one-anchor term is exactly

\[
\sum_{x=0}^{w-1}EI_{(0,x)}=\nu_{w,H}.                        \tag{5.1}
\]

The higher connected-cluster correction is bounded by

\[
R_w
\le\operatorname{poly}(w)e^{-2\kappa w}.                  \tag{5.2}
\]

Since

\[
\nu_w=e^{-\kappa w+o(w)},                                    \tag{5.3}
\]

we have

\[
\boxed{R_w=o(\nu_w).}                                        \tag{5.4}
\]

Therefore the localized anchor void pressure `J_{w,H}` satisfies

\[
\boxed{J_{w,H}=\nu_{w,H}[1+o(1)].}                           \tag{5.5}
\]

Together with (2.4),

\[
J_{w,H}=\nu_w[1+o(1)].                                       \tag{5.6}
\]

## 6. From anchor void pressure to the safe-strip free energy

The strip free energy

\[
I^0_{G,w}
=-\lim_{m\to\infty}\frac1m
\log P(\text{no horizontal occupied homology in an }m\text{-row strip})\tag{6.1}
\]

is insensitive to `O(H)` boundary rows.

Insert empty guard rows at the two strip ends.  On the guarded interior, every horizontal essential component is a **complete** cylinder component and therefore has one anchor; conversely every localized complete winding component creates strip horizontal homology.  Removing/adding the `O(H)` margins changes `-log` probability by at most a boundary cost independent of `m`, which vanishes after division by `m`.

The only discrepancy is a component whose vertical span exceeds `H`; its per-row activity is the tail in (2.2), already `o(nu_w)` by the choice of `C`.

Hence

\[
I^0_{G,w}
=J_{w,H}+o(\nu_w).                                           \tag{6.2}
\]

Combining with (5.6) proves (1).

## 7. Quantitative form

The argument yields the schematic estimate

\[
\boxed{
I^0_{G,w}
=\nu_w
+O(\operatorname{poly}(w)e^{-2\kappa w})
+O(w e^{-cCw}).}                                              \tag{7.1}
\]

The polynomial and constants have not been optimized.  Since the leading activity is only `e^{-kappa w+o(w)}`, both errors are relatively negligible after choosing `C` as in (2.3).

This is a pressure/activity statement, not an OZ prefactor theorem: it does not determine the polynomial prefactor of `nu_w` itself.

## 8. Consequences

### 8.1 Poisson births and void Perron roots are the same dilute gas to first order

The component-Poisson theorem uses `nu_w` as its natural clock.  The fixed-width charge theory uses `I^0_w` as the exponential no-winding cost.  Equation (1) shows that in the fixed-subcritical large-width regime

\[
\boxed{\text{activity }\nu_w
\quad=\quad
\text{pressure }I^0_w\,[1+o(1)].}                            \tag{8.1}
\]

Their finite-width difference is precisely a higher connected-cluster correction of the winding-component gas.

### 8.2 The distinction reappears near criticality

At the charge/eigenvalue root, `p_w->p_c` and the scaled width `w/xi(p_w)` is `O(1)`.  The anchor gas is no longer exponentially dilute in `w`.  There is no reason for `I^0_w/nu_w->1` in that regime.

Thus the rare-component Poisson language and the critical transfer/CFT language are not competing descriptions.  They are two regimes of the same object, separated by loss of small activity.

### 8.3 Cluster corrections have a concrete meaning

The first correction to `I^0=nu` is generated by connected pairs/triples of nearby winding components.  A future second-order calculation can therefore target a winding-component virial coefficient rather than an unexplained finite-width discrepancy.

This may be a cleaner route to finite-width corrections than fitting `I^0/nu-1` directly.

## 9. Literature boundary

Connected-cluster expansions for dependency graphs / abstract polymer gases are classical; Scott--Sokal (2003) and cluster-expansion Lovasz-local-lemma refinements give general nonvanishing/convergence frameworks.  The percolation-specific content here is the component-anchor localization and the `k`-fold disjoint winding-witness BK bound that drives the activity parameter exponentially small.

No claim of novelty is made for the generic pressure/activity expansion.

## 10. Claim boundary

The proof uses the existing author-level cylinder localization, full-component anchor construction, fixed-subcritical mass rate, and site BK inputs from #739.  The connected-cluster estimate is standard but should receive a line-by-line combinatorial audit if promoted into the final manuscript.  The result is only for fixed subcritical `p`; it is not a near-critical uniform theorem.
