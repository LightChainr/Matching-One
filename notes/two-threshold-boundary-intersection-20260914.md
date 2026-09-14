# Two threshold boundaries and their shared pivotal intersection

Date: 2026-09-14. Exact finite boundary calculus. Addendum to draft #773 / #768 / #769 / #778.

## 1. The two rank births are exits from two nested monotone complexes

Use the nested simplicial complexes

    Delta_0={S:r(S)=0},
    Delta_1={S:r(S)<=1}.

A monotone insertion path exits `Delta_0` at the first rank birth and exits `Delta_1` at the second. For one addressed site v and a configuration of all other sites, opening v can produce only

    0->0, 0->1, 0->2, 1->1, 1->2, 2->2.

Define the p-dependent local boundary masses

    alpha(p)=P_p(r(v=0)=0, r(v=1)=1),
    beta(p) =P_p(r(v=0)=0, r(v=1)=2),
    gamma(p)=P_p(r(v=0)=1, r(v=1)=2).                  (1)

`beta` is exactly the shared boundary: opening v exits both nested complexes at once.

## 2. Exact Russo boundary decomposition

For a transitive N-site torus, multilinearity / Russo differentiation gives

    -P0'(p) = N [alpha(p)+beta(p)],                      (2)
    +P2'(p) = N [beta(p)+gamma(p)].                      (3)

The rank-one sector itself satisfies the exact continuity equation

    P1'(p)=N [alpha(p)-gamma(p)].                         (4)

Therefore the matching-function derivative is

    M'(p)=P2'(p)-P0'(p)
         =N [alpha+2 beta+gamma]
         =N E_p[Delta_v r].                             (5)

The factor 2 multiplying beta is not a convention: a direct `0->2` insertion crosses both rank thresholds simultaneously. Equation (4) is equally important: the local matching-odd imbalance `alpha-gamma` is exactly the slope of the rank-one probability, not merely a heuristic arm observable.

This identifies the earlier jump-two share of M' as the normalized shared-boundary mass

    omega_2(p)
      = 2 beta / [alpha+2 beta+gamma].                   (6)

## 3. Random-chain interpretation

For a random permutation, the direct-jump flux at cardinality k is the fixed-cardinality version of beta. Summing over k gives

    P(D=0)=sum_k B_k
          =N integral_0^1 beta(p) dp.                    (7)

Thus `D=0` is the probability that a uniformly random maximal chain of the Boolean lattice hits the intersection of the two threshold boundaries without spending a level in the shell `Delta_1\Delta_0`.

The ordinary persistence gap D is the number of discrete levels spent in this shell.

This gives a purely finite combinatorial meaning to simultaneous topological birth, independent of any arm-language interpretation.

## 4. Matching parity of the three boundary pieces

Under primal/matching complement, the first and second threshold boundaries exchange:

    alpha_G(p) = gamma_Ghat(1-p),
    gamma_G(p) = alpha_Ghat(1-p),
    beta_G(p)  = beta_Ghat(1-p).                         (8)

Hence the shared intersection beta is matching-even, while `alpha-gamma` is matching-odd.

In a self-matching model at its symmetric point,

    alpha=gamma,                                         (9)

and beta is the independent shared-boundary channel.

This is the local-boundary version of the process parity result that persistence gap is even and midpoint is odd.

## 5. Arm-fusion interpretation is now sharply typed

At criticality, ordinary exits from either threshold are expected to be four-arm pivotal events. The shared boundary beta requires one site to complete **two independent homology ranks at once**.

The exact attachment-spine controls show two mechanisms, including abundant three-attachment theta spines. If the theta mechanism corresponds to a six-arm event, then

    beta(p_c) ~ L^(-alpha_6),

whereas the total pivotal boundary scales like `L^-alpha_4` per site. Therefore

    omega_2(p_c) ~ L^[-(alpha_6-alpha_4)] = L^-5/3.       (10)

If the shared-boundary event instead truly requires eight arms,

    omega_2(p_c) ~ L^-4.                                 (11)

Equations (10)--(11) are the same discrimination proposed earlier, but now beta has an exact finite boundary definition.

Crucially, beta is matching-even by (8). Therefore a six-arm law for beta does not by itself generate an L^-5/3 correction in the matching-odd one-point balance function.

## 6. A new odd local diagnostic

The local matching-odd boundary imbalance is

    delta_boundary(p)=alpha(p)-gamma(p).                 (12)

Equation (4) makes this observable completely explicit:

    delta_boundary(p)=P1'(p)/N.                          (13)

Under matching,

    delta_boundary^G(p)
      = -delta_boundary^Ghat(1-p).                       (14)

At a self-matching symmetric point it vanishes exactly. For square/matching pair, its finite-size scaling at the paired critical points is a more relevant local test for a low-dimensional **odd** arm channel than beta itself.

This suggests an addition to the pivotal atlas: report alpha, beta, gamma separately, not only the weighted rank jump `E Delta r`.

If a six-arm geometry appears only in beta while `alpha-gamma` has no six-arm contribution, the 8-arm odd-root hypothesis survives that local test.

### 6.1 Exact relation to the intrinsic rank-shape slope

At a balance point write `P0=P2=a`, `P1=r=1-2a`, and use the canonical coordinates

    b=1/2 log(P0/P2),
    d=log[P1/sqrt(P0P2)].

Substituting (2)--(4) gives

    b'(p) = -N[alpha+2 beta+gamma]/(2a),
    d'(p) =  N[alpha-gamma]/(2ar).                       (15)

Therefore the parameterisation-free slope of the physical rank curve is

    dd/db
      = - [alpha-gamma]
          / [r (alpha+2 beta+gamma)].                    (16)

This is exactly the previously derived `d log c/db=-2h1/g1`, now expressed as a local boundary-flux ratio.

Equation (16) sharply separates the two fusion questions:

- `beta` measures the **even shared boundary** and controls direct simultaneous births;
- `alpha-gamma` measures the **odd imbalance of the two separate boundaries** and is the local source of asymmetry in the self-normalized rank curve.

Thus an observed six-arm scaling in beta is compatible with an 8-arm-controlled odd curve. A six-arm kill of the 8-arm odd hypothesis would need six-arm scaling in `alpha-gamma` (or another matching-odd observable), not merely in beta.

## 7. Interfaces

**#769:** from the same fixed-site before/after rank table, report p-dependent or cardinality-resolved alpha/beta/gamma. No extra homology engine is needed. The existing jump2 histogram is beta; split the rank-jump-one cases into `0->1` and `1->2`.

**#768:** use `alpha-gamma`, not beta, as the local matching-odd lower-arm diagnostic.

**#778:** beta integrates to the zero-persistence atom; alpha/gamma locate entry/exit of nonzero persistence intervals.

**#775:** the rank-sector coefficient arrays determine the total boundary marginals `alpha+beta` and `beta+gamma`; #769 beta then completes the local boundary decomposition.

## 8. Boundaries

- Arm identifications in sections 5--6 remain scaling hypotheses; equations (1)--(9),(12)--(16) are exact finite identities.
- `alpha-gamma` is a local odd diagnostic, not automatically the same observable as the global matching-function correction at p_c.
