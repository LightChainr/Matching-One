# Research frontier after the structural reductions — 2026-09-14

Status: proof programmes and conjectures unless a subsection explicitly says otherwise. This note is intended to consume the structural results in `structural-consequences-20260914.md`, not to create a new queue of parallel production tasks.

## 1. Homological free energy: energy minus opportunity entropy

For a period lattice `Lambda` of area/index `N` and a primitive period class `u in Lambda/±`, let `tau_p(u)` be the subcritical directional connection cost. A class-`u` winding component has a natural number of transverse opportunities of order

\[
h_u\asymp N/|u|.
\]

If its once-per-component intensity has a form

\[
\nu_u(p)\approx C(p,\hat u)|u|^{-\beta(p,\hat u)}e^{-\tau_p(u)},
\]

then the leading logarithm of the expected number of such rare components is

\[
\log h_u-\tau_p(u)+O(\log|u|).
\]

This suggests the variational barrier

\[
\boxed{
\Psi_\Lambda(p)=
\min_{u\in\Lambda_{\rm prim}/\pm}
\left\{\tau_p(u)-\log\frac{N}{|u|}\right\}.}
\]

The leading first-rank birth should occur near `Psi_Lambda(p)=0`; the component prefactor controls only the finer centre correction.

For the axial torus with periods `(w,0),(0,m)`, the cheapest class is `(w,0)`, so `tau_p(u)=w kappa(p)` and `N/|u|=m`. The condition `Psi≈0` reduces exactly to

\[
w\kappa(p)\approx\log m,
\]

which is the centre relation already used in the exponential-aspect analysis.

### Deterministic simplification in exponential elongation

If `u` is a Euclidean shortest period with `|u|=ell`, `h=N/ell`, then any nonparallel `v in Lambda` obeys `|det(u,v)|>=N` and hence `|v|>=h`. At fixed subcritical `p`, norm equivalence gives `tau_p(v)>=c_p h`, while its opportunity entropy is at most `log(N/|v|)<=log ell`. Under `log h/ell -> d>0`, nonparallel classes therefore have exponentially larger energy cost and cannot compete with `±u` in this variational problem. This reinforces the directional-separation lemma already recorded for #765.

The genuinely new regime is a lattice deliberately tuned so that two or more primitive classes have `tau_p(u)-log(N/|u|)` within `O(1)`. Nonparallel essential components cannot be independent species because intersection/topological rank constraints couple them. A useful language is a **homological hard-core gas**: parallel species can coexist; nonparallel species force intersections/mergers and tend to accelerate the rank-two birth.

Falsifiable prediction: for a family with a finite set of candidate primitive classes, the first birth is controlled by `min_u Psi_u=0`; if two nonparallel classes are nearly tied, the rank-one plateau should shorten relative to a single-minimizer geometry. If a class with significantly larger `Psi_u` dominates, the opportunity factor or microscopic normalization in the conjecture is wrong.

## 2. Small-p actual-SITE prefactor theorem via a column transfer operator

The current branch already supplies three pieces of a possible proof:

1. `dilute-winding-crossover.md` proves for the actual NN site model that, when `w p^2 -> 0`,
   \[
   \nu_w(p)=p^w I_0(2wp)(1+o(1)),
   \]
   and therefore exhibits a genuine `w^{-1/2}` closure factor in the subregime `wp -> infinity`, `wp^2 -> 0`;
2. `sewing-with-memory.md` writes the **complete-component** activity using exact three-column local weights, so distinct external boundary-site weights are not fundamentally nonlocal;
3. the tagged all-height resolvent gives an exact finite-`w` representation of the complete component and identifies mixing, one transverse mode, and Palm normalization as the real missing hypotheses.

This suggests proving the first fixed-`p` complete-component prefactor theorem on a nonempty interval `0<p<p_0`, rather than trying to solve the full subcritical phase at once.

### Candidate state and operator

Scan around the circumference column by column. A canonical state should retain:

- the occupied row sets in the two most recent columns, modulo a common vertical translation;
- the exposed connectivity partition;
- horizontal lift/winding data;
- enough root/unrooting information that a closed complete component is counted once.

The row coordinate is unbounded, so the state space is countable. For sufficiently small fixed `p`, large vertical excursions have exponentially small activity. Work on a weighted Banach space penalizing slice diameter/vertex count and aim to show the transfer `T_p` is bounded and quasi-compact with a simple leading eigenvalue.

Introduce a transverse Fourier twist `theta`. Reflection symmetry should give

\[
\log\lambda_p(\theta)
=\log\lambda_p(0)-\frac12D(p)\theta^2+O(\theta^4).
\]

Closing the component after `w` columns requires zero net transverse displacement, hence a Fourier coefficient

\[
\frac1{2\pi}\int_{-\pi}^{\pi}
\lambda_p(\theta)^w A_p(\theta)\,d\theta.
\]

A spectral gap plus a nonlattice local CLT/saddle estimate would give

\[
\boxed{
\nu_w(p)=C(p)w^{-1/2}e^{-w\kappa(p)}(1+O(1/w))
}
\]

on `0<p<p_0`, with `C(p)` automatically containing the genuine complete-component unrooting and boundary weights. This would be qualitatively stronger than fitting more widths and would directly answer part of #740 for an actual SITE model.

### Possible simultaneous payoff for #760

In the same cluster-expansion/spectral regime, the distinction between a plane connection mass and a finite-width cylinder visible mass should be generated only by wrap defects whose activity requires `Omega(w)` occupied steps. The natural target is

\[
0\le\gamma_w(p)-\kappa(p)\le C(p)e^{-c(p)w}
\]

or at least `o(1/w)` on `0<p<p_0`. Turning the rare-wrap statement into a leading-eigenvalue perturbation estimate is the real theorem; a union bound alone is not enough.

A successful small-`p` package could therefore provide, on one strict interval, all of: `beta=1/2`, a real `D(p)`, an HK/Brownian span limit, cylinder/plane locality, and `p`-analyticity. That is more informative than another sparse width ladder.

## 3. Marked-Poisson and component-Palm as the default morphology interface

The long-torus component process should not be the primitive object for every new geometric observable. The existing localized anchor proof can carry finite marks at no conceptual cost. Therefore morphology questions should be reorganized around the single complete component selected under component Palm.

Recommended hierarchy:

1. define the full component and the exact mark (`span`, vertex-biconnected winding core, occupation count, distinct boundary count, seam multiplicity, etc.);
2. prove/measure its Palm law or LDP;
3. use the marked-Poisson extension to lift it to a rare cloud on an exponentially long torus.

This avoids repeatedly rediscovering the same dependence control and makes event-conditioned versus component-conditioned sampling distinctions explicit.

## 4. The loop/branch frontier should be used as a falsifiable geometry theorem

The #758 candidate rate, if correct, necessarily has one fixed transition `r_*(p)` between a bulging-loop regime and a linear branch regime. This gives a low-cost falsification protocol before attempting the full network upper bound:

- the effective large-`A` slope must tend to `kappa(p)`;
- the transition location is set by the directional norm, not by the chosen rare-event threshold;
- small-`A` curvature must agree with the `D^{-1}=partial_yy tau` relation if the Brownian range limit is also correct.

A rigorously observed slope below `kappa` would identify a genuinely cheaper network topology and would be more valuable than another finite-width goodness-of-fit comparison.

## 5. Near-critical common window: use `(rank,K)`, not a bivariate copula

For the same parameter and same labels, exact topology reduces `(W_4,W_8)` to the state `(rank,K)`. Therefore the common-window scaling object should be formulated as

- the rank-sector weights `P_0,P_1,P_2`, and
- the positive-integer rank-one count PGF `H_x(z)`.

The one-sided constrained-Poisson laws provide boundary conditions for this object. A near-critical interpolation should match those boundaries and the critical/finite-aspect regime while preserving exact support. It should **not** be built by continuously deforming the two separated-window independent Poisson variables into a common-parameter bivariate Poisson law.

The existing speculative scaling `w nu_w(p) -> Phi(w kappa(p))` can be retained as a marginal hypothesis, but it is insufficient to determine the joint black/white crossover. The missing object is the topologically constrained family `H_x` together with the rank weights.

## 6. Execution order

To keep the programme from proliferating again:

1. integrate the persistent birth reflection, exact same-parameter PGF reduction, and dual-even/odd birth coordinates into the structural layer;
2. make the alternating-barrier span comparison publication-grade and close the reciprocal mean white-span theorem; then upgrade to the exponential/Gamma laws through the existing point-process proof;
3. record marked-Poisson as the common interface for morphology work;
4. pursue the small-`p` transfer-operator theorem as the next genuinely new proof technology;
5. keep the homological free-energy and near-critical `(rank,K)` pictures as unifying maps, not immediate production queues.

No extra width-16/20/24 prefactor campaign is justified by this note alone, and no continuum-field interpretation is implied.
