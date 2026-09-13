# From a complete-component resolvent to a closed diffusive theory

2026-09-13. Same #739 probability programme. The exact finite constructions
are in `tagged-span-resolvent.md`; this note separates an elementary
conditional stability theorem from three deliberately unproved site
conjectures. No extra width ladder is commissioned.

## 1. A complete component's range is not an additive observable

If a component is decomposed into pieces at transverse positions y_i, with
local lower/upper excursions ell_i,u_i, its inclusive span is

    L=1+max_i(y_i+u_i)-min_i(y_i+ell_i).                       (1)

It is not sum_i(u_i-ell_i+1), nor the total variation of the center path.
An alternating center sequence 0,1,0,1,... has range one and total variation
proportional to its length. This deterministic distinction invalidates the
piecewise-additive-span explanation in the returned diagnostic; it does
not invalidate the computed finite-width moments.

Here is the appropriate transfer lemma.

**Conditional bush-stability theorem.** Suppose a complete component C_w
has a center path gamma_w and transverse Hausdorff error B_w, with

    E[number of pieces whose error exceeds t] <= C w^a exp(-ct)             (2)

under the **same component-Palm/closed-loop ensemble**, for fixed C,a,c>0.
Suppose gamma_w/sqrt(Dw) has a Brownian-bridge limit, with the moment bounds
required at the order under discussion. Then C_w has the same transverse
range limit and the same corresponding normalized moments.

Proof. By a union bound,

    Pr(B_w>t)<=min(1,Cw^a exp(-ct)),

so for each fixed r, E B_w^r=O((log w)^r). Deterministically,

    |range(C_w)-range(gamma_w)| <= 2 B_w.                    (3)

After division by sqrt(w), (3) tends to zero in every fixed L^r. Slutsky
and the stated center-path uniform integrability transfer the law and
moments. In particular, if the skeleton first two moments have their
Brownian limits with error e_w, then the component CV^2 differs from the
skeleton CV^2 by O(log w/sqrt(w)), plus that pre-existing error. No sign
for this correction is asserted.

This proof requires neither independence of pieces nor absence of branches.
The load-bearing condition is (2) under the **closed component law**.
Unconditional exponential tails cannot simply be conditioned on an
exponentially rare winding event and reused unchanged.

CIV section 1.3.3, equation (1.10), states O(log n) Hausdorff closeness of
its full long connection cluster to the effective path, followed by its
Brownian-bridge Theorem C. Thus invoking CIV to say the full span must be
an extensive sum is the opposite of what that text establishes in its
actual setting. Its bond model/open-endpoint conditional law remains
different from this site's periodic component-Palm law. It is a proof
architecture and a specific missing comparison, not a finished transfer.

Finite-width CV values 0.127 (NN) or 0.116 (matching) versus the target
0.04719755 do not decide this theorem's hypotheses. Even simple positive
boundary thickness alone is not a complete explanation: adding a positive
constant to a Brownian range decreases CV^2. Finite core corrections,
random decorations and covariance can produce either correction sign.

## 2. Conjecture R: formulate the fixed-p problem using the exact resolvent

For each fixed subcritical p of the chosen graph, let alpha_w,R_w,b_w be
the tagged construction and nu_w=alpha_w(I-R_w)^(-1)b_w. Let mu_w=E L,
which is now obtainable without a height cutoff. Define

    F_w(s)= z alpha_w(I-zR_w)^(-1)b_w / nu_w,
    z=exp(-s/mu_w),  s>=0.                                  (4)

**Conjecture R (complete shape, mean normalized).** At fixed strictly
subcritical p,

    F_w(s) -> E exp[-s range(B)/E range(B)]                   (5)

locally uniformly for s>=0, with convergence of the corresponding moments.
No mass, amplitude, diffusion coefficient or finite-height offset is fitted
in (5). It is equivalent to a mean-normalized full shape law, stronger than
agreement of a single CV statistic. Fixed-w rationality does not preclude
a nonrational limit as the number of states grows. Conversely, the largest
single eigenvalue of R_w does not by itself determine the distribution in
the joint w,height scaling window; source/exit overlaps and multiple modes
remain in (4).

The exact Brownian-range CDF used to compute the right-hand side is

    F_R(x)=sqrt(2pi)*pi^2/x^3 sum_(n>=1) n^2 exp[-pi^2 n^2/(2x^2)]
          =1+2 sum_(n>=1)(1-4n^2x^2)exp(-2n^2x^2),            (6)

with E R=sqrt(pi/2), E R^2=pi^2/6. The explicit target transform at s=1,2,4
is approximately 0.3762840878410, 0.1475796893947, 0.0252164362526.
These are Brownian formula values, not measured percolation quantities.

In the new all-height computations, NN p=1/4 gives at s=2:

    w=4:0.1875595963; w=6:0.1744603845; w=8:0.1671020675.

Matching p=1/8 gives 0.1765806063,0.1690681784,0.1638984502.
These move toward the proposed value over the reported widths; they do not
prove its limit and do not establish a correction exponent. The same finite
data can be compatible with several extrapolations. No new computation is
ordered merely by making the conjecture explicit.

A proof can aim directly at the normalized resolvent (4), or prove the
conditional bush-stability hypotheses with a correctly sewn loop skeleton.
It need not first invent a unique raw seam cut. An obstruction should identify
which of mixing, one transverse mode, Palm normalization, or tight decoration
control fails.

## 3. Conjecture J: joint geometry and thermodynamic marks

The direct-activity representation defines, near a fixed physical point,

    Psi_w(z,u,v)=sum_C z^L u^K v^B,

where B is the number of DISTINCT external boundary sites. This is not a
random-cluster q derivative. It supplies the correct finite joint law and
covariance of (L,K,B), all within the complete-component ensemble.

**Conjecture J (one closed diffusive band with additive marks).** If the
physical closed chain has a reflection-symmetric mixing renewal description,
then for constants D>0, densities rho_K,rho_B and a covariance matrix Sigma,

    (L/sqrt(Dw), (K-rho_K w)/sqrt(w), (B-rho_B w)/sqrt(w))
           => (range(Brownian bridge), Gaussian_2),          (7)

with independence between the range and the Gaussian mark vector.
Possible sublinear corrections to the centering must be controlled; actual
finite means are a safer alternative where that expansion is unknown.

Why independence is plausible, rather than just convenient: in a
reflection-symmetric Markov-additive effective chain, transverse increments
are odd and occupancy/boundary increments are even. Their long-run cross
covariances vanish. A joint functional CLT would give a Gaussian transverse
path independent of the additive even marks; conditioning the path to close
produces a bridge and preserves that independence. Applying this argument
requires a proved mixing/conditioning comparison for the site loop ensemble,
not merely reflection of the one-point marginal.

The prior decorated dilute limit supplies one rigorous special-regime
motivation, not a proof at fixed p. The present exact computations give
finite squared correlations Corr(L,K)^2 of about 0.867 at NN w=4,p=1/4
and 0.838 at w=5. They are not small yet. They neither establish (7) nor
contradict an asymptotic zero correlation. Under the moment hypotheses of
(7), a concrete consequence would be Corr(L,K)->0 whenever the limiting
occupation variance is nonzero.

If f(u,v)=lim_w w^(-1) log Psi_w(1,u,v) exists with enough local regularity,
its log-coordinate gradient gives rho_K,rho_B and its Hessian gives Sigma.
Along u=p,v=1-p this would imply

    -kappa'(p)=rho_K/p-rho_B/(1-p).                          (8)

Equation (8) is conditional on exchanging the limit and derivative.
The finite score identity is exact already; analyticity of finite rational
functions does not prove p-analyticity of the infinite limit.

## 4. Conjecture U: is the closed-loop residue actually one?

This is the bolder, less-supported conjecture. Keep it separate from R and J.
For an explicitly defined finite-state cyclic renewal kernel A(z,k), a
single simple critical eigenvalue gives

    w[z^w] {-log det(I-A(z,k))} ~ R(k)^(-w).

The simple logarithmic singularity has coefficient one. Integrating the
one transverse quadratic mode gives

    L_w ~ exp(-kappa w)/sqrt(2pi D w),                      (9)

with no independent endpoint-overlap amplitude. This is a statement about
that specified cyclic object, whose w/n unrooting is exact; it is not a
statement about arbitrary site clusters.

**Conjecture U (unit sewing residue).** A canonical all-regeneration-cut
construction for the actual complete site component might identify it with
such a cyclic Gibbs object, making the previous sewing factor zeta(p)=1.
All external-boundary overlaps and multiple-cut weights must be included
before making this identification. A merely bounded, nontrivial cyclic
insertion would leave zeta!=1; several leading bands, residual marking
multiplicity, or lattice periodicity could also invalidate (9).

This conjecture has a stronger falsifiable consequence than beta=1/2. If R
and the corresponding diffusion scale also hold, then

    2 exp(kappa w) nu_w E L -> 1.                            (10)

Without U, the same quantity tends to zeta(p), not necessarily one. With
only the older heat-kernel conjecture, nu_w E L has no power prefactor but
still has the unknown zeta/2 amplitude. An accurately known mass is needed
to use (10); none is manufactured from the six existing densities.

This is not a new compute order. It states which genuinely new normalization
identity would remove an unknown amplitude, and exactly how a counterexample
would change the theory. The direct activity transfer now defines the correct
microscopic object to compare, rather than a count of seam crossings.

## 5. Near-critical continuation stays downstream

The fixed-subcritical conjectures above are not uniform claims for p->pc.
Even if D(p)kappa(p)->D0 and a scaling function for w nu_w were established,
the same-parameter black/white count constraint forbids blindly taking the
independent two-window Poisson law into a common critical window. The prior
conditional log-aspect/log-log-aspect crossover remains a conjecture requiring
uniform input; the new finite matrices do not discharge it.

Current direction: use the exact all-height component law as the microscopic
object; establish or falsify the correct closed-loop mixing and decoration
control; only then infer a fixed-p amplitude or a critical crossover. This
keeps the geometric probability paper coherent while preserving the bolder
hypotheses in an explicit, testable form.

Sources: CIV, arXiv:math/0610100v2, section 1.3.3 and equation (1.10), for the
open-connection architecture under its bond hypotheses; standard phase-type
absorption formulas are prior art (Maier, 1991, doi:10.1080/15326349108807207).
No general Brownian, phase-type, Gibbs or renewal theorem is claimed new here.
