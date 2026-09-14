# Parameter-regularity audit for the inverse correlation length

2026-09-14.  Targeted literature/logic audit after the earlier correction that **direction analyticity is not p-analyticity**.  This note records what may and may not be used to remove the at-most-countable exceptional `d` set in the median-centred Gumbel theorem.

## 1. What the current #739 proof already establishes for actual SITE

For each of the NN and matching SITE models on a compact subcritical parameter interval, the component-activity argument gives

\[
F_w(p)=\frac1w\log\nu_{w,w^2}(p)\to-\kappa_G(p),              \tag{1.1}
\]

with a uniform lower curvature bound

\[
F_w''(p)\ge-C_I.                                               \tag{1.2}
\]

Hence `-kappa_G` is locally semiconvex, equivalently `kappa_G` is locally semiconcave.  In one parameter this implies:

- local Lipschitz regularity;
- existence of one-sided derivatives everywhere;
- differentiability except at at most countably many points;
- moving-point derivative convergence of `F_w'` at every differentiability point.

Together with the quantitative strict parameter comparison already on the branch, every differentiability-point slope is strictly negative.

This is enough for the regular-`d` affine Gumbel theorem and the constrained subsequential description at corners.  It does **not** prove that the corner set is empty.

## 2. Campanino--Chayes--Chayes 1991: p-analyticity, but for bond percolation

Campanino, J. Chayes and L. Chayes, *Gaussian fluctuations of connectivities in the subcritical regime of percolation*, Probab. Theory Relat. Fields 88 (1991), study the `d`-dimensional **Bernoulli bond percolation model**.  Their abstract explicitly lists, throughout the subcritical regime:

1. the `L^{-(d-1)/2}` Ornstein--Zernike prefactor along an axis;
2. real analyticity of the correlation length as a function of the bond parameter;
3. a local limit theorem for the conditioned long cluster.

Thus there is genuine classical precedent for **parameter analyticity of the Bernoulli bond mass**.  This is stronger than mere direction analyticity.

However, the printed model is bond percolation.  No direct theorem statement for independent SITE on the square or matching graph was identified in this audit.  A site-to-bond gadget representation with deterministic edges is not automatically within the hypotheses of a homogeneous bond theorem, and should not be used without checking the proof class.

## 3. Campanino--Ioffe--Velenik: analytic in direction

The CIV random-cluster fluctuation theory proves, under its stated subcritical assumptions, sharp OZ asymptotics plus analyticity/strict convexity of the inverse correlation length as a function on the **direction sphere/Wulff boundary**.

This is precisely the source of the earlier repository correction: that theorem supports directional smoothness at fixed model parameter, not by itself real analyticity of `p -> kappa_G(p)` for the square SITE models.

The paper mentions parameter/inverse-temperature regularity in surrounding discussion, but the theorem needed here must be stated for the actual parameter/model before the exceptional set is removed.

## 4. Analyticity of susceptibility/local observables is not enough

There are modern results proving analyticity of the subcritical susceptibility for transitive Bernoulli percolation and uniform analyticity of many local FK observables under mixing assumptions.  These do not automatically imply analyticity of the exponential rate

\[
\kappa(p)=-\lim_n\frac1n\log P_p(0\leftrightarrow ne_1).       \tag{4.1}
\]

The limit is a nonlocal large-distance object; exchanging analytic limits uniformly requires additional OZ/transfer control.

Therefore susceptibility/local-event analyticity is useful proof technology but is **not** a citation-level closure of the SITE mass regularity gap.

## 5. Current verdict

For the actual NN and matching independent SITE models used in #739:

\[
\boxed{
\text{Do not remove the at-most-countable exceptional p/d set yet.}}
\]

The new strict enhancement result

\[
\kappa_8(p)<\kappa_4(p)                                       \tag{5.1}
\]

is logically independent of differentiability.  It proves strict graph/centre separation but does not rule out corners of either mass.

Similarly, strict convexity/analyticity in **direction** does not rule out corners in the scalar occupation parameter.

## 6. Two credible closure routes

### Route A: locate a direct SITE p-analyticity theorem

The ideal input is a theorem explicitly covering independent site percolation on finite-range quasi-transitive graphs throughout the subcritical interval, with real-analytic dependence of the inverse correlation norm on the occupation parameter.  A bond-only statement is insufficient unless its proof is shown to cover the site gadget representation.

### Route B: prove it from the SITE renewal/operator representation

The small-`p` countable-state programme on this branch would give p-analyticity on a nonempty interval if:

- the complete/open connection is represented by a quasi-compact positive operator depending analytically on `p`;
- the relevant Perron eigenvalue is simple and isolated;
- the physical root is given by an implicit Perron equation.

The analytic implicit-function theorem then gives `p`-analyticity of `kappa` there.  Extending this throughout the whole subcritical regime is essentially an SITE OZ regularity theorem rather than a one-line corollary.

## 7. Consequence for #763

The correct current hierarchy remains:

1. intensity-clock Poisson/Gumbel coordinates: no p-differentiability needed;
2. true finite-median `1/w` affine Gumbel: proved at differentiability points;
3. exceptional `d`: at most countable, with constrained convex-log-intensity subsequential limits;
4. full elimination of exceptional `d`: open until an actual SITE p-regularity theorem/proof is supplied.

No relation involving `a(d)+b(d)=1` should be substituted for this regularity question; the strict centre result now gives `a+b>1` and the differentiability issue remains separate.

## 8. Sources checked

- M. Campanino, J. Chayes, L. Chayes, *Gaussian fluctuations of connectivities in the subcritical regime of percolation*, PTRF 88 (1991): abstract/model statement is Bernoulli **bond** percolation and explicitly states real analyticity of the correlation length in the subcritical parameter.
- M. Campanino, D. Ioffe, Y. Velenik, *Fluctuation theory of connectivities for subcritical random cluster models*, Ann. Probab. 36 (2008): OZ/random-walk structure and analyticity/strict convexity of the inverse-correlation shape/direction under its assumptions.
- A. Georgakopoulos, C. Panagiotis, *Analyticity results in Bernoulli Percolation* (2018/2020 versions): subcritical susceptibility/local analytic results, not a direct SITE inverse-mass theorem.

A negative search is not an impossibility result; this note only prevents an unsupported promotion.
