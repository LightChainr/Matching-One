# Matching parity of the full two-birth persistence process

Date: 2026-09-14. Exact pathwise finite identity. Addendum to draft #773 / #768 / #778.

## 1. Couple the primal and matching-complement processes by the same labels

Let `G` be the primal NN site graph and `Ghat` its matching complement graph on the same honest torus. Give every site an iid continuous label `U_v~Uniform[0,1]`.

The primal occupied set at parameter p is

    B_p={v:U_v<=p}.

Define matching labels

    V_v=1-U_v.

Then the matching occupied set at parameter q is

    What_q={v:V_v<=q}={v:U_v>=1-q}.

Ignoring probability-zero equality at a label, `What_q` is exactly the complement of `B_(1-q)`.

The digital-Alexander matching identity therefore holds **pathwise** along the coupled processes:

    r_G(p) + r_Ghat(1-p) = 2.                            (1)

## 2. The two rank births reverse under matching

Let

    T1^G = inf{p:r_G(p)>=1},
    T2^G = inf{p:r_G(p)=2},

and define the corresponding increasing-q births for `Ghat`. Since the primal rank path is

    0  for p<T1^G,
    1  for T1^G<=p<T2^G,
    2  for p>=T2^G,

(1) gives pathwise

    T1^Ghat = 1-T2^G,
    T2^Ghat = 1-T1^G.                                   (2)

This is stronger than a distributional relation.

## 3. Gap is matching-even; midpoint is matching-odd

Define raw persistence gap and centered midpoint

    G_p = T2-T1,
    M_p = (T1+T2-1)/2.

Equation (2) yields pathwise

    G_p^Ghat = G_p^G,
    M_p^Ghat = -M_p^G.                                  (3)

Thus the **entire gap distribution is a matching-even process observable**, while the centered midpoint is matching-odd.

For a self-matching lattice such as triangular site percolation, (3) is an internal symmetry of one model.

## 4. Discrete random-permutation form

For N sites, reverse a uniformly random insertion permutation to obtain the matching-complement insertion order. If `J1,J2` are the primal rank-birth indices, then pathwise under this coupling

    J1^Ghat = N+1-J2^G,
    J2^Ghat = N+1-J1^G.                                 (4)

Hence

    D^Ghat=J2^Ghat-J1^Ghat = D^G.                       (5)

In particular

    P_G(D=0)=P_Ghat(D=0)                                (6)

exactly at every finite size.

## 5. Cardinality-resolved direct-jump flux has the same parity

Let `B_k^G` be the step-k flux of direct rank `0->2` births. The complement/reversal map gives the coefficient relation

    B_k^G = B_(N-1-k)^Ghat.                             (7)

Equivalently, for a fixed addressed site v, the number of k-subsets of the other N-1 sites with `Delta_v r_G=2` equals the number of `(N-1-k)`-subsets with `Delta_v r_Ghat=2` after complementing the other sites.

In Bernoulli form,

    P_p^G(Delta_v X=2)
      = P_(1-p)^Ghat(Delta_v X=2).                       (8)

Thus the integrated simultaneous-birth probability is automatically matching-even.

## 6. Canonical b-coordinate gives the same decomposition

The finite matching relation for sector probabilities implies

    b_Ghat(1-p) = -b_G(p),

where `b=1/2 log(P0/P2)` is defined separately for each graph. Put

    B1=b(T1),
    B2=b(T2).

Combining with (2):

    B1^Ghat = -B2^G,
    B2^Ghat = -B1^G.                                    (9)

Therefore canonical gap and midpoint satisfy pathwise

    (B1-B2)^Ghat = (B1-B2)^G,
    [(B1+B2)/2]^Ghat = -[(B1+B2)/2]^G.                  (10)

The full persistence copula naturally decomposes into matching-even gap information and matching-odd location/midpoint information.

## 7. Consequence for the 6-arm / 8-arm question

The exact L=3,4 controls on #773 show that direct rank-jump-two configurations contain many theta/T3 spines, plausible six-arm fusion geometries. Equations (3),(6)--(8) show that **this simultaneous-birth observable is matching-even by construction**.

Therefore:

> the existence, abundance, or even six-arm scaling of direct `0->2` births does NOT by itself imply a six-arm contribution to the matching-odd one-point balance function `M_G(p_c)`.

The strong 8-arm root-correction hypothesis is threatened only if a six-arm channel survives in a **matching-odd** observable/sector. #768 must therefore distinguish:

    six-arm absolute/process geometry            (allowed, matching-even),
    six-arm matching-odd torus one-point sector  (the actual kill condition). (11)

This materially narrows the theoretical obligation.

## 8. Process-level universality decomposition

For any matching pair, a future canonical two-birth scaling limit should inherit

    G*^G =_law G*^Ghat,
    M*^G =_law -M*^Ghat.                                (12)

At the common continuum percolation fixed point, primal and matching microscopic distinctions should disappear after the correct universal identification, suggesting an internally symmetric midpoint law and a common even gap law.

This process parity is stronger than comparing only `P0,P1,P2` at one time.

## 9. Interfaces

**#768:** do not use the abundant theta/T3 `D=0` class as a direct refutation of the 8-arm matching-odd hypothesis. First identify whether the corresponding six-arm field is odd or even under matching/map exchange.

**#769:** report the cardinality-resolved `B_k` sequence; (7) is an exact cross-model regression if matching data are available.

**#778:** raw and canonical gap distributions should agree exactly between a primal/matching pair under the coupled convention; midpoint distributions should reflect. This is a strong paired-data validation.

**#776:** self-matching triangular site specializes (2)--(10) to exact within-model reflection identities.

## 10. Boundaries

- The pathwise identity uses continuous labels to avoid tie conventions; discrete permutation form is exact under reversal.
- Matching-even six-arm process geometry can still influence even observables and higher source responses.
- The note does not prove that every six-arm continuum field is matching-even; it proves that the specific direct-rank-jump persistence channel is.
