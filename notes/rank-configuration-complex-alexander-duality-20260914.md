# Rank-configuration complexes and simplicial Alexander duality

Date: 2026-09-14. Exact finite combinatorial-topology addendum to draft #773. No scaling or CFT input is used.

## 1. Two nested configuration complexes

For a finite N-site honest torus and primal connectivity graph G, define

    Delta_0^G = {S subset V : r_G(S)=0},
    Delta_1^G = {S subset V : r_G(S)<=1}.

Because homology rank cannot increase when occupied sites are removed, both are simplicial complexes and

    Delta_0^G subset Delta_1^G.

The rank-one sector is the shell `Delta_1\Delta_0`; it is not itself a simplicial complex.

For a simplicial complex Delta on V, use the standard Alexander dual

    Delta^* = {T subset V : V\T notin Delta}.              (1)

## 2. Exact matching Alexander duality

The configuration-level digital Alexander relation is

    r_G(S)+r_Ghat(V\S)=2.                                  (2)

Therefore

    T in (Delta_0^G)^*
    <=> r_G(V\T)>=1
    <=> r_Ghat(T)<=1,

and similarly for Delta_1. Hence

    boxed: (Delta_0^G)^* = Delta_1^Ghat,                   (3)
    boxed: (Delta_1^G)^* = Delta_0^Ghat.                   (4)

For a self-matching lattice, Delta_0 and Delta_1 are Alexander dual on the same vertex set.

This theorem contains more information than sector-count reciprocity: it pairs the entire face posets, not only the number of faces at each cardinality.

## 3. f-polynomial reciprocity follows automatically

Let

    F_Delta(z)=sum_(S in Delta) z^|S|.

From the definition of Alexander dual,

    F_(Delta^*)(z)
      = (1+z)^N - z^N F_Delta(z^-1).                       (5)

Since

    Z_0^G=F_(Delta_0^G),
    Z_0^G+Z_1^G=F_(Delta_1^G),

(3)--(5) immediately reproduce

    Z_j^G(z)=z^N Z_(2-j)^Ghat(z^-1),                       (6)

and the reciprocal bivariate rank enumerator

    Q_G(z,u)=z^N u^2 Q_Ghat(z^-1,u^-1).                    (7)

Thus the Newton-polygon reciprocity is the face-enumerator shadow of simplicial Alexander duality.

## 4. Minimal witnesses and maximal safe configurations

Alexander duality exchanges minimal nonfaces with complements of facets.

Consequently:

- minimal nonfaces of `Delta_0^G` are inclusion-minimal occupied sets with nonzero homology rank: primitive rank-one winding witnesses;
- their complements are facets of `Delta_1^Ghat`: maximal matching configurations which still have rank at most one;
- minimal nonfaces of `Delta_1^G` are inclusion-minimal rank-two witnesses;
- their complements are facets of `Delta_0^Ghat`: maximal contractible matching configurations.

This gives an exact certificate duality for first and second rank birth. It also explains the low-fugacity Newton edges: on an LxL square torus the smallest rank-one witness has L sites, while a smallest rank-two witness has `2L-1` sites (one row and one column sharing a site).

## 5. Homology of the configuration complexes

Combinatorial Alexander duality gives, over a field,

    H~_i(Delta) ~= H~^(N-i-3)(Delta^*).                   (8)

An independent F2 boundary-matrix computation at L=3 gives the following nonzero ordinary Betti numbers (beta_0 includes connectedness; the higher entries are reduced as usual).

### Square NN primal / square matching 8-neighbour dual

    Delta_0^NN:       beta_0=1, beta_3=1, beta_4=6
    Delta_1^NN:       beta_0=1, beta_4=1, beta_5=2
    Delta_0^matching: beta_0=1, beta_1=2, beta_2=1
    Delta_1^matching: beta_0=1, beta_2=6, beta_3=1

The higher Betti numbers pair exactly under `i <-> N-i-3` as required by (3)--(4).

### Triangular self-matching L=3

    Delta_0: beta_0=1, beta_2=4, beta_3=2
    Delta_1: beta_0=1, beta_3=2, beta_4=4.

Again the two complexes are exact Alexander duals.

These Betti numbers are topology of the configuration-space complexes, not physical-space cluster Betti numbers and not new continuum observables by themselves.

## 6. Threshold boundaries and direct double birth

Let `partial^+ Delta` denote Boolean-lattice edges leaving a monotone complex. The first rank birth is an exit from Delta_0; the second is an exit from Delta_1. A direct `0->2` insertion lies in the intersection of the two upward boundaries.

Thus the jump-two flux / shared pivotal mass has a clean poset interpretation:

    direct double birth = partial^+ Delta_0 intersect partial^+ Delta_1.  (9)

Under Alexander duality, upward primal boundary data map to downward boundary/facet data of the matching complexes. This may provide a more efficient certificate/classification route for #769 than identifying arms from scratch.

## 7. Process flag enumerators

The shell `Delta_1\Delta_0` contains exactly the rank-one configurations. Nested chains inside this shell determine persistence moments:

    F_m(k_1,...,k_m)
      = # {S_1 subset ... subset S_m:
           |S_a|=k_a, S_a in Delta_1\Delta_0}.            (10)

These are flag-enumerator data of an Alexander-dual pair of complexes. This is the natural combinatorial object behind #786's rank-one persistence hierarchy.

## 8. New research questions

1. Can the minimal-nonface/facet duality produce a practical exact DP for direct rank births or rank-sector polynomials?
2. Do the configuration complexes have useful shellability/Cohen--Macaulay properties on honest tori? No such property is assumed; one-variable sector polynomials already fail ordinary real-rootedness.
3. Can flag-enumerator reciprocity reduce the paired-birth computation by working on whichever of the primal/dual complexes has fewer faces?
4. Which configuration-complex Betti numbers stabilize with width, and do any have a direct transfer/closure interpretation?

These are structural/computational questions. No claim is made that configuration-space Betti growth is a physical critical exponent.

## 9. Boundaries

- Equations (3)--(7) are exact consequences of the existing configuration-level matching relation.
- The L=3 Betti tables are exact F2 controls only.
- Simplicial Alexander duality here is distinct from the physical-space digital Alexander theorem used to establish (2); the latter is the model input, the former is the configuration-complex consequence.
