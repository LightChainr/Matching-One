# Rank-one homology direction is a persistent label of the birth interval

Date: 2026-09-14. Exact finite process/topology bridge. Addendum to draft #773 / #775 / #778 / #585.

## 1. The rank-one direction cannot change before the second birth

Let `A_k` be the image of the occupied graph's first homology in ambient torus `H_1(T^2)` after k monotone site insertions. Adding occupied sites preserves every old cycle, hence

    A_k subseteq A_(k+1).                                 (1)

Suppose the rank process is one at two consecutive times. Then both `A_k` and `A_(k+1)` are one-dimensional subspaces and (1) forces equality.

Therefore, on every trajectory with nonzero persistence gap

    D=J2-J1>0,

the primitive rank-one subgroup is constant for the entire interval

    k=J1,...,J2-1.                                       (2)

Denote this unoriented primitive subgroup by

    V={a,b}, gcd(a,b)=1, {a,b}={-a,-b}.

The same argument applies to the continuous Uniform-label filtration on `[T1,T2)`.

Thus the rank process naturally carries a **marked persistence bar**

    (T1,T2,V)                                             (3)

rather than only two scalar birth times.

## 2. Direction-resolved one-time probabilities are interval coverage probabilities

Define

    P_v(p)=P(r(p)=1, V=v).                                (4)

Then pathwise

    1_{r(p)=1,V=v}=1_{T1<=p<T2,V=v}.

Integrating over p gives the exact identity

    boxed: integral_0^1 P_v(p) dp
      = E[(T2-T1) 1_{V=v}].                              (5)

Summing over v recovers

    integral_0^1 P1(p) dp=E[T2-T1].                      (6)

So a static direction-resolved rank law is not the unweighted distribution of birth directions: it is a time-slice/coverage law, and its integral is persistence-length weighted.

## 3. Fixed-cardinality version

Let

    C[v,k]
      = # {S subset V_sites: |S|=k, rank(S)=1,
                            primitive subgroup v}.        (7)

A random permutation is uniform over k-subsets at time k, hence

    P(R_k=1,V=v)=C[v,k]/binom(N,k).                       (8)

Since D is the number of insertion levels spent at rank one,

    boxed: sum_k C[v,k]/binom(N,k)
      = E[D 1_{V=v}].                                    (9)

This is the direction-resolved refinement of the existing exact formula for `E D`.

Consequently #775's optional `C[L,a,b,k]` extension would immediately deliver the persistence-length-weighted direction distribution without paired trajectory storage.

## 4. Birth-direction distribution requires entry flux, not occupancy mass

The unweighted probability that a persistent interval is born with direction v is

    P(V=v,D>0)
      = sum_k A_{k,v},                                   (10)

where `A_{k,v}` is the direct `0->1` entry flux at insertion level k, resolved by the newly created primitive subgroup.

Thus there are three different direction laws:

1. **time-slice law:** `P_v(p)/P1(p)` at fixed p;
2. **length-biased law:** proportional to the left side of (9);
3. **birth law:** proportional to the entry flux (10).

They coincide only under additional independence which is not assumed.

This distinction matters when comparing finite-lattice rank-one directions with continuum Arguin/Pinson probabilities: the latter are fixed-parameter time-slice probabilities.

## 5. Matching complement preserves the rank-one subgroup

The configuration-level digital Alexander theorem gives

    A_G(S)^perp = A_Ghat(V\S)                             (11)

with respect to the torus intersection pairing. In a two-dimensional symplectic vector space, a one-dimensional isotropic subspace is its own symplectic orthogonal:

    span(v)^perp=span(v).                                 (12)

Therefore whenever primal and matching configurations are both rank one under complement, their primitive subgroup is the same:

    boxed: V_G(S)=V_Ghat(V\S).                            (13)

Consequently

    P_v^G(p)=P_v^Ghat(1-p).                              (14)

This is a direction-resolved strengthening of aggregate rank matching.

## 6. Full marked persistence matching duality

Using the paired Uniform-label construction and the already proved birth-time reflection,

    T1^G=1-T2^Ghat,
    T2^G=1-T1^Ghat,

while (13) preserves V. Hence pathwise

    boxed: (T1^G,T2^G,V^G)
      = (1-T2^Ghat,1-T1^Ghat,V^Ghat).                    (15)

For a self-matching lattice,

    (T1,T2,V) =_law (1-T2,1-T1,V).                       (16)

Thus persistence length and direction are matching-even, while the centered interval midpoint is matching-odd even after conditioning on V.

For any function f(V), self-matching gives

    E[ midpoint_centered * f(V,D) ]=0                    (17)

when integrable. It does not imply independence of midpoint, duration and direction.

## 7. Homology-character source as a process transform

For the orientationless character

    chi_theta(V)=cos(a theta_1+b theta_2),

define the direction-resolved one-time source

    H_p(theta)=E[1_{r(p)=1} chi_theta(V)].                (18)

Then

    boxed: integral_0^1 H_p(theta) dp
      = E[(T2-T1) chi_theta(V)].                          (19)

In the random-permutation chain,

    sum_k E[1_{R_k=1} chi_theta(V_k)]
      = E[D chi_theta(V)].                               (20)

So the homology-direction source has both a modular one-time interpretation and a persistence-process interpretation.

## 8. New continuum targets

At critical fixed modulus, Arguin/Pinson determines the time-slice direction law. The marked persistence process suggests two stronger universal targets:

1. the canonical joint law of `(B1,B2,V)` where `B_i=b(T_i)`;
2. the duration-biased modular transform

       E[G_b chi_theta(V)],                              (21)

   which is independent of thermal reparameterisation.

These separate whether long rank-one intervals preferentially carry axis, diagonal or more oblique homology classes.

In elongated subcritical geometries, deterministic directional cost separation predicts concentration of V onto the cheapest period. Hence `(B1,B2,V)` is a natural object connecting critical torus homology to the rare-winding regime.

## 9. Interfaces

- #775: retain primitive direction in rank-one transfer states if cheap; equations (9) and (19) become immediate outputs.
- #769: resolve `0->1` entry flux by new subgroup to obtain the unweighted birth-direction law (10).
- #778: if paired raw trajectories include the first winding vector, report joint `(J1,J2,V)` rather than only `(J1,J2)`.
- #585: use Arguin/Pinson direction law as an explicit map/connectivity positive control.
- #765: compare the critical direction law with deterministic off-critical period-cost selection.

## 10. Boundaries

- Constancy of V on a rank-one interval is exact finite monotonic topology.
- Arguin/Pinson provides critical continuum time-slice probabilities, not the marked birth copula.
- Direction V is a subgroup label, not a local field or chirality; the unoriented convention is essential.
