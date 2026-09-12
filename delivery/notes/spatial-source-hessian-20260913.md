# A complete finite spatial source Hessian, not another amplitude fit

Date: 2026-09-13. Exact finite analysis. No new production, width scan, field identification, or novelty claim.

## 1. Observable, field, and root coordinate

Let a finite translation-transitive square-cell torus have N vertices. Let X=r-1, and

    M(p_1,...,p_N)=E[X].

Independent Bernoulli occupation makes M separately affine in every p_i. At a uniform
interior root p_*, let M(p)=M(p,...,p), and write M', M'' for UNIFORM thermal derivatives.
Strict rank monotonicity gives M'(p_*)>0. Introduce additive site offsets a_i and define
q(a) by M(q(a)+a_1,...,q(a)+a_N)=0 near a=0. This is a finite matching root, not p_c.
All perturbed probabilities are required to remain in (0,1).

Translation invariance, the chain rule, and separate affinity imply exactly

    M_i=M'/N,       M_ii=0,       sum_j M_ij=M''/N.

Implicit differentiation therefore gives

    q_i(0)=-1/N,
    R_ij := q_ij(0) = -M_ij/M' + M''/(N^2 M').                 (1)

All quantities on the right are evaluated at the unperturbed finite root. In particular,
R*1=0. This is also forced by the exact local identity q(a+c*1)=q(a)-c.

For EVERY zero-mean spatial pattern h, sum_i h_i=0,

    q'(0;h)=0,
    q''(0;h)=-h^T M_site'' h / M'.                             (2)

This first-order null needs only transitive translations, not a reflection making h odd.
Consequently a vanishing single-source response for a zero-mean field is NOT by itself
an identification of a nontrivial selection rule or continuum sector.

For a SINGLE addressed site v the stronger exact identity holds at any uniform p:

    M(p,...,p+epsilon at v,...,p)=M(p)+epsilon*M'(p)/N.       (3)

Equation (3) is not a truncation: no second power of that site's probability exists.
It supplies a cheap independently checkable local-defect calibration.

## 2. Conditional covariance representation; the coordinate cannot be omitted

Use log-odds z=log(p/(1-p)) and a physical product field eta:

    p_i(z,eta)=logistic(z+eta_i),
    Z_j(z,eta)=sum_{r(omega)=j} exp(z*K(omega)+sum_i eta_i n_i(omega)).

The common product normalizer cancels in the ODDS ONLY. Define

    g(z,eta)=log Z_2(z,eta)-log Z_0(z,eta).

Let Delta denote rank-2-conditioned minus rank-0-conditioned expectation. Direct
finite-sum differentiation gives

    g_z = Delta E[K],
    g_eta_i = Delta E[n_i],
    g_eta_i_eta_j = Delta Cov(n_i,n_j).                     (4)

At the translation-invariant root, Delta E[n_i]=Delta E[K]/N and g_z>0. The root z_*(eta)
satisfies z_i=-1/N. For any zero-mean h put H=sum_i h_i n_i. Then

    d^2 z_*(epsilon h)/d epsilon^2|0
      = - [Var(H|r=2)-Var(H|r=0)] / Delta E[K].              (5)

Thus the mean-preserving local logit susceptibility is a DIFFERENCE OF CONDITIONAL
STRUCTURE FACTORS. It is not a generic-Q derivative, and no overall normalizer has
been guessed. The denominator is positive and is the conditional occupation-number
separation already identified in the rank-conditioning work.

Equation (5) concerns the logit root. For the root expressed back in p, multiply by
p(1-p), because z'(0)=0 for this h.

Do not interchange additive probability and additive logit perturbations. With
N^-1 sum h_i^2=1 and sum h_i=0, exact coordinate conversion yields

    q''_logit-field (in baseline p)
      = [p(1-p)]^2 q''_additive-p-field + p(1-p)(2p-1).     (6)

The last term is the mean second-order probability change under the nonlinear logit
map. It is present even for mean-zero eta. Numerically similar source definitions are
not the same experiment. Equation (6) is verified independently against (5).

For arbitrary h, the full logit-root Hessian is

    R^z_ij = -[Delta Cov(n_i,n_j)-Delta Var(K)/N^2]/Delta E[K].

This also annihilates the uniform direction in z coordinates. The same uniform
logit field has nonzero p-root curvature p(1-p)(1-2p), because p=logistic(z).
The zero-mean table below is not a claim that the uniform p-root curvature vanishes.

## 3. Complete 4x4 result

A standalone physical lifted-graph traversal enumerates all 65,536 site configurations.
Ranks have totals 36,559 / 19,932 / 9,045. For every displacement d and every remaining
occupation count k, the mixed derivative polynomial is constructed from

    Delta_0 Delta_d X(A)
      = X(A+0+d)-X(A+0)-X(A+d)+X(A),

where A excludes the two addressed vertices. Coefficients have the unnormalized
Bernstein form sum_k c_d[k] p^k(1-p)^(14-k); the diagonal is identically zero.
The output also stores rank-conditioned origin/pair count polynomials, rather than
only fitted or rounded susceptibilities.

Translation makes the Hessian convolutional. For wavevector k=(kx,ky) in (Z/4)^2,

    lambda_k(p)=sum_d M_0d(p) cos(pi*(kx*dx+ky*dy)/2).

The Fourier diagonalization is exact and integer-valued at the coefficient level.
At p=1/2, the six D4 orbits have eigenvalues

    k=(0,0):  7429/8192,  multiplicity 1;
    k=(1,0):  2665/8192,  multiplicity 4;
    k=(2,0):  3733/8192,  multiplicity 2;
    k=(1,1): -2483/8192,  multiplicity 4;
    k=(2,1): -2567/8192,  multiplicity 4;
    k=(2,2): -5355/8192,  multiplicity 1.

Their weighted sum is zero, as required by separate affinity. The uniform eigenvalue
is M''/N, NOT zero for M itself. It becomes zero for the ROOT Hessian in (1).

The finite root has the exact isolating bracket stored in full-site-hessian.json,
with midpoint approximately 0.5906721123310283. At this root the unit-RMS, zero-mean
mode curvatures are:

| k representative | multiplicity | additive probability field q'' | logit field q'' (root reported in p) |
|---|---:|---:|---:|
| (1,0) | 4 | -1.0498968867066683 | -0.0175285451213038 |
| (2,0) | 2 | -1.9320342433611963 | -0.0690955391643194 |
| (1,1) | 4 | +0.7525801283645556 | +0.0878386301804579 |
| (2,1) | 4 | +0.7978317203732835 | +0.0904838968817806 |
| (2,2) | 1 | +2.5198959892681610 | +0.1911503942542182 |

The additive-root Hessian has inertia (9 positive, 6 negative, 1 zero).
The signs of the additive curvatures use exact rational enclosures over the entire
root bracket. The finite torus root is therefore NEITHER convex NOR concave on the
mean-preserving field subspace. Equal source variance/strength does not specify its
response: spatial correlations and wavevector matter.

For real unit-RMS representative fields use h(x,y)=cos(pi*k.x/2)+sin(pi*k.x/2).
For every nonzero orbit representative in this table its entries are +/-1, the mean
is zero, and its norm squared is 16. One translation sends h to -h, so q(epsilon h)
is exactly even. The direct two-colour polynomial is independently reconstructed
from the PR708 automaton, and finite-amplitude roots at epsilon=1/64,1/128,1/256
approach the predicted curvature with O(epsilon^2) central-difference error.

Example: striped h=(-1)^x and checkerboard h=(-1)^(x+y) have the same mean and RMS.
Their additive root curvatures are -1.93203424 and +2.51989599 respectively. These are
full finite-system responses, not a claim about the relevance of random disorder
or an infinite-volume anisotropic critical surface.

## 4. Verification and scope

- Independent automaton versus lifted physical graph: all 65,536 ranks agree.
- Independent Bernoulli likelihood-score calculation: all 15 off-diagonal entries
  agree at p=1/3,1/2,2/3 (45 exact rational comparisons).
- Thermal sum rule N*sum_d M_0d=M'' holds as an integer polynomial identity.
- Conditional covariance identity (5) and coordinate conversion (6) agree to the
  specified midpoint precision; no floating covariance rank is used.
- Direct finite-amplitude bivariate polynomials are independent of the derivative
  enumeration. Their roots are high-precision controls, not interval proofs.

Files: scripts/torus_source_hessian.py, scripts/verify_source_hessian.py and the two
full-site-hessian result JSONs. Only the optional independent verifier consumes the
unchanged PR708 certificate (Git blob 50b7297deefe7c50215aea2ed534ca5810461af3).
No new source has been inserted into #275's frozen original-U comparison.
