# Rank threshold complexes: exact Alexander duality without a matroid rank

Date: 2026-09-14. Finite combinatorial addendum to draft #773 / #775 / #636.

## 1. Two natural simplicial complexes

For a finite honest torus with vertex set V and primal graph G, define

    Delta_0^G = {S subset V : r_G(S)=0},
    Delta_1^G = {S subset V : r_G(S)<=1}.                (1)

Because homology rank cannot increase when occupied vertices are removed, both are downward-closed abstract simplicial complexes.

Their f-vectors are exactly the rank-sector coefficient data:

    f_k(Delta_0^G) = C_G[0,k],
    f_k(Delta_1^G) = C_G[0,k]+C_G[1,k].                  (2)

Thus the exact transfer target in #775 is simultaneously the pair of f-vectors of two nested complexes.

## 2. Matching complement is combinatorial Alexander duality

For a simplicial complex Delta on V, use

    Delta^*={T subset V : V\T notin Delta}.               (3)

The digital-Alexander rank identity

    r_G(S)+r_Ghat(V\S)=2                                 (4)

gives immediately

    (Delta_0^G)^* = Delta_1^Ghat,
    (Delta_1^G)^* = Delta_0^Ghat.                        (5)

Proof:

    T in (Delta_0^G)^*
    iff r_G(V\T)>0
    iff r_Ghat(T)<2
    iff T in Delta_1^Ghat,

and similarly for the second identity.

At the coefficient level this implies the reversed f-vector relations

    f_k(Delta_1^Ghat)
      = binom(N,k)-f_(N-k)(Delta_0^G),                   (6)

    f_k(Delta_0^Ghat)
      = binom(N,k)-f_(N-k)(Delta_1^G).                   (7)

These contain the familiar rank-sector complement identities but package them as an actual Alexander-dual pair of complexes.

## 3. Minimal nonfaces are the topological seeds

The minimal nonfaces of `Delta_0^G` are the inclusion-minimal occupied vertex sets which create nonzero ambient homology. The minimal nonfaces of `Delta_1^G` are the inclusion-minimal rank-2 supports.

Therefore:

- the smallest minimal nonface size of `Delta_0` is the first possible rank-1 onset;
- the smallest minimal nonface size of `Delta_1` is the first rank-2 onset;
- the number of minimum-size minimal nonfaces is the geometric onset coefficient used in low-fugacity anchoring.

Under (5), these minimal generators are dual to complements of facets of the matching threshold complex. This gives a finite algebraic dictionary for the geometric anchor terms used in conditional-odds integration.

## 4. Rank is not a matroid or polymatroid rank

It is tempting to treat `S -> r_G(S)` as a low-rank matroid/polymatroid function. This is false even at the smallest honest examples.

On the L=3 square torus take one horizontal row and split it as

    A={site 0},
    B={sites 1,2}.

Both A and B are contractible and have rank zero; their intersection is empty. Their union is the complete periodic row and has rank one. Hence

    r(A)+r(B)=0
      < 1 = r(A union B)+r(A intersection B),             (8)

violating submodularity.

The same three-site counterexample appears in the standard L=3 triangular torus.

The obstruction is the essential geometry of winding: two individually contractible pieces can close each other across a periodic seam. Therefore matroid exchange/submodular-rank algorithms cannot be assumed for the topological rank process.

## 5. What survives from simplicial-complex theory

Although the rank function is not submodular, `Delta_0` and `Delta_1` are honest simplicial complexes. Their f-vectors must obey every general simplicial f-vector constraint, in particular Kruskal--Katona shadow inequalities.

This supplies an independent exact validation layer for future rank-sector transfer output:

1. verify `C0[k]` is an f-vector;
2. verify `C0[k]+C1[k]` is an f-vector;
3. verify the matching-side arrays are their Alexander-dual reversed f-vectors via (6)--(7);
4. reconstruct minimum nonface counts and compare with direct minimal-winding enumeration at small L.

A violation is a transfer/homology bug, not a new physical effect.

## 6. A new algebraic direction

Let `I_0^G` and `I_1^G` be the Stanley--Reisner ideals of the two threshold complexes. Then (5) identifies the primal/matching pair through Alexander dual ideals. Potentially useful finite invariants include:

    minimal generator degree/counts,
    graded Betti tables,
    projective dimension / regularity,
    facet/minimal-seed incidence.

These should be treated as **finite combinatorial descriptors**, not continuum fields. Their value would be to compress and certify topological onset/closure structure, especially if the resolutions stabilise or admit a transfer description.

No claim of such stabilization is made here.

## 7. Interfaces

**#775:** add Kruskal--Katona/f-vector and Alexander-dual reversal checks to every exact C[j,k] output. They are independent of the transfer implementation.

**#636:** minimal nonfaces/generators provide a finite closure/topology dictionary which does not assume a matroid or a generic operator algebra.

**#778:** the random permutation rank process is the random vertex filtration through the nested pair `Delta_0 subset Delta_1`; the two birth indices are the exit times from these two complexes.

## 8. Boundaries

- Simplicial Alexander duality here is a finite combinatorial consequence of the digital matching-rank theorem; it is not a new continuum Alexander-duality claim.
- Stanley--Reisner invariants are not automatically physically meaningful scaling fields.
- The explicit submodularity counterexample rules out matroid/polymatroid shortcuts unless additional structure is imposed on a restricted state class.
