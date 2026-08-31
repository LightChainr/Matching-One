# Can thermal and closed-source reparameterization absorb one saturation defect?

This is a finite, exact model-closure question using already published integer coefficients. It is not an independent experiment, a new observable search, or a prospective statistical confirmation. The scalar-gain and thermal-only hypotheses rejected at `bc17b81d` stay rejected; this calculation asks whether the more general endpoint source family itself can describe the infinitesimal interior profiles.

## Fixed family and necessary identity

Parent N50 has Gaussian generators `(5,5)` and `(1,7)`, homogeneous B probability p and A probability `1-epsilon*(1-p)`. The physical source remains bulk `exp(t*S)`, with `S=C+F4+Bvac=C+F4+T_NN-4K+2N`. Each geometry is normalized separately. The reference is epsilon=t=0, at the already enclosed common matching root, p0=1 minus the fixed child N25 root. Geometry order, source units and saturation chart do not change.

Let `f=(<q>_first,<E>_first,<q>_second,<E>_second)`, where q=rank−1 and E=q². These are the original readouts used to build the matching root and U. Define four-vector columns at the common p0:

```text
T = partial_p f
C = partial_t f = Cov(f,S) separately in each geometry
H = partial_epsilon f = 25*(1-p)*(<f>_one_hole - <f>_intact)
```

The candidate is first-order closure of **both complete q/E profiles in both geometries** within the same endpoint family, after arbitrary common smooth coordinate changes:

```text
f(p,epsilon,t) = f_endpoint(phi(p,epsilon,t), psi(p,epsilon,t)) + O(epsilon²)
phi(p,0,t)=p; psi(p,0,t)=t.
```

The functions may depend on p and t; there is no source-independent-gain assumption. At the reference point, the chain rule requires `H=b*T+c*C`, with the same b,c for all four rows. Thus a rank-three tangent matrix `[T,C,H]` rules out every such map, including p-dependent effective source couplings. Equality of just the scalar U is a strictly weaker question and is not excluded by this criterion alone. Separate maps for each observable or geometry are outside the candidate.

## One fixed primary decision

Use the complete 4×3 matrix in the row order above. Compute all four 3×3 minors and their sum of squares, `D3=sum_I det([T,C,H]_I)^2`. No row is chosen after seeing its value. Positive invertible changes of units are unnecessary: zero versus nonzero is invariant under any invertible fixed row transformation and nonsingular reparameterization of the two endpoint columns. D3's numerical magnitude itself is basis dependent.

- A strictly positive rational lower bound for D3 rejects thermal-plus-S endpoint profile closure at this finite root. At least one extra finite response direction is needed beyond those two tangents; this is not a count of continuum fields.
- Otherwise report that this fixed root enclosure does not exclude the candidate. Zero at one root would not prove profile closure on a neighborhood. No new source, geometry, occupancy point, defect sector, data or model is introduced to force an exclusion.

Root-comoving columns `C-T*mean(C_q)/mean(T_q)` and `H-T*mean(H_q)/mean(T_q)` are also displayed to express the same obstruction after eliminating the common matching-root motion. They are derived from the same four-vector matrix, not additional tests or evidence.

Column addition leaves every third-order minor unchanged, so the root-comoving construction gives exactly the same D3. A common density source at the endpoint adds `c*p*(1-p)*T` to C and also leaves D3 unchanged. This follows because only the B occupations fluctuate at saturation; it does not discard the density term from the fixed physical source or the previously computed mixed U jet. No area exponent or cos4 normalization is needed for this rank decision.

The old mixed-U residual is insufficient to decide this candidate. For example, a pure source reparameterization `U_int(epsilon,t)=u((1+lambda*epsilon)*t)` remains inside the endpoint source family but has `U*U_(epsilon,t)-U_epsilon*U_t=lambda*u(0)*u'(0)`, possibly nonzero, while H=0 and D3=0. The new criterion tests common profile closure; it does not reinterpret the old gain rejection or count its coefficients again as independent evidence.

## Exact calculation and limits

Inputs are the four existing complete Bernstein coefficient tables and saved root enclosure at `f5c4a74a`, copied with hashes in SOURCES.json. Endpoint parent coefficients are the complemented/reversed child coefficients; q and qS reverse sign, E and S do not. Each free-B occupancy row must contain binomial(25,k) configurations, and per-K sums already include multiplicity. There is no second enumeration, no random sampling, and no cloud work.

Use rational interval evaluation of the original Bernstein sums, source centering and first p derivative. Verify the saved root endpoints by exact opposite signs and its positive pooled slope. Outward rounding is for serialization only. The result is conditional on the published finite graph counts. Different quotient Smith classes, endpoint-only scope and the previous P154/P334/F4 stop decisions remain unchanged.

The contract, pinned inputs and implementation are committed before evaluating D3. Because the underlying exhaustive block and its earlier U result were already public, this chronology does not turn the new question into an independent confirmation.
