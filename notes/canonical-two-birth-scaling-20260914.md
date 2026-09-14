# Canonical two-birth scaling: L^(5/4) insertion gaps and complement-symmetric persistence

Date: 2026-09-14. Process-level conjecture and exact self-matching symmetry. Addendum to draft #773 / #778.

## 1. Why the insertion-index gap has scale L^(5/4)

On an LxL torus, N=L^2 sites receive iid Uniform[0,1] labels. Let T1<=T2 be the continuous occupation parameters at which rank first reaches 1 and 2, and let J1<=J2 be the corresponding insertion indices in the random permutation.

The percolation thermal window has width

    p-p_c = O(L^-3/4).

A p-window of this size contains, in insertion-count units,

    N * L^-3/4 = L^(5/4)                               (1)

sites. Therefore a nondegenerate two-birth near-critical limit predicts

    J_i - N p_c = O_P(L^(5/4)),
    D=J2-J1 = O_P(L^(5/4)).                             (2)

The existing exact controls already show the right raw scaling direction:

    square:    E D / L^(5/4) = 0.37992... (L3), 0.38795... (L4),
    triangular:E D / L^(5/4) = 0.37992... (L3), 0.38502... (L4).

Two sizes are not an exponent measurement. Equation (1) is the actual scaling argument.

## 2. Canonical b-coordinate removes the thermal metric entirely

For each finite lattice define

    b_L(p)=1/2 log[P0_L(p)/P2_L(p)].

It is strictly decreasing from +infinity to -infinity. Put

    B1=b_L(T1),
    B2=b_L(T2),
    B1>=B2.                                               (3)

The one-time rank law is exactly the interval-coverage data of the random canonical persistence interval `[B2,B1]`:

    P0(b) = P(B1 < b),
    P1(b) = P(B2 < b <= B1),
    P2(b) = P(B2 >= b),                                  (4)

up to irrelevant endpoint conventions.

Hence the marginals of B1 and B2 are fixed by the self-normalized rank curve, while their copula is the genuinely new process information.

The canonical gap

    G = B1-B2 >=0                                        (5)

has exact mean

    E G = integral_R P1(b) db.                           (6)

No p_c or thermal metric appears in (3)--(6).

## 3. Near-critical two-birth universality conjecture

**Conjecture P1.** At fixed torus modulus, after passing to the canonical b coordinate,

    (B1,B2) => (B1*,B2*)                                 (7)

with a lattice-independent ordered-pair law in the percolation universality class.

This is stronger than one-time rank-law universality: the limiting marginals are already encoded by `P_j^*(b)`, but (7) additionally specifies the persistence copula.

Consequences include universal laws for

    G*=B1*-B2*,
    midpoint M*=(B1*+B2*)/2,
    two-time coverage H*(u,v)=P(B2*<u<v<=B1*).

The raw insertion process should satisfy the same limit after one nonuniversal linear thermal conversion, explaining (2).

## 4. Exact complement symmetry in a self-matching model

For a self-matching lattice such as triangular site percolation, complementing all labels by

    U_v -> 1-U_v

preserves their joint law and exchanges occupied with vacant topology. The finite rank-complement identity sends the two births to each other in reverse order. Therefore

    (T1,T2) =_law (1-T2, 1-T1).                           (8)

Because the self-matching canonical coordinate is exactly odd around 1/2,

    b_L(1-p)=-b_L(p),

(8) becomes the exact finite process symmetry

    (B1,B2) =_law (-B2,-B1).                             (9)

In midpoint-gap variables,

    M=(B1+B2)/2,
    G=B1-B2,

this is simply

    (M,G) =_law (-M,G).                                  (10)

Thus every odd-in-M process moment vanishes exactly, including

    E M=0,
    Cov(M,G^k)=0                                         (11)

for every integrable k.

This is a much stronger self-matching regression than one-time `P0(p)=P2(1-p)`.

## 5. Discrete permutation version of the self-matching symmetry

For N sites, reversing a uniformly random permutation is again uniform. In a self-matching model,

    (J1,J2) =_law (N+1-J2, N+1-J1).                     (12)

The exact triangular L3 control satisfies

    E J1=17/4,
    E J2=23/4,
    E J1+E J2=10=N+1.

Its full joint table can be checked cellwise against (12). The square NN process is not self-matching and does not satisfy this same-model identity; its partner relation is with the matching-complement process instead.

## 6. Square/matching pair process symmetry

For a primal graph G and its matching complement Ghat, the correct finite relation is cross-model:

    (B1^G,B2^G) =_law (-B2^Ghat,-B1^Ghat),              (13)

provided both use their own canonical b coordinates and the same torus geometry. Therefore the midpoint distributions reflect while the gap distributions coincide:

    G^G =_law G^Ghat.                                    (14)

Equation (14) is a strong future check on complete paired production or exact rank-sector transfer data. It is not visible from one model's one-time rank curve alone.

## 7. Relation to simultaneous births and arm fusion

The event `D=0` is exactly a zero-length discrete persistence interval. Under (2), any continuum atom at zero would survive in the scaled process. The current arm picture predicts instead that

    P(D=0) -> 0,

with the rate diagnosing the local fusion:

    6-arm candidate: L^-5/3,
    8-arm candidate: L^-4.                               (15)

Thus process convergence and the jump-two arm question meet at the boundary behaviour of the canonical persistence-gap law near zero.

An abundant 6-arm microscopic mechanism can still coexist with an L^-4 one-point matching-root correction if its contribution cancels in the matching-odd projection. The canonical copula is the right place to separate these statements.

## 8. A rigorous triangular-site route

Garban--Pete--Schramm construct the near-critical triangular-site configuration process in the quad-crossing topology. If periodic rank events are shown to be continuity events under an appropriate torus/quotient version of that framework, then (7) for triangular site becomes a concrete theorem target rather than an unconstrained universality slogan:

1. construct the coupled near-critical process;
2. define the first lambda at which torus rank reaches 1 and 2;
3. prove almost-sure continuity/no-ambiguous-birth at the hitting levels;
4. push through the canonical rank coordinate.

The missing periodic-domain/homology continuity lemma should be stated explicitly if no published theorem supplies it.

## 9. Interfaces

**#778:** paired raw `(J1,J2)` should be converted to both raw `D/L^(5/4)` and canonical `(B1,B2)` where the finite rank polynomials are available. For self-matching controls verify (12) cellwise.

**#776:** triangular exact/transfer continuation should report process-symmetry controls if paired permutation data can be obtained without new production; the theoretical near-critical interface is the route in section 8.

**#768/#769:** the small-gap tail / D=0 atom supplies a process-level arm-fusion diagnostic independent of the one-point root correction.

## 10. Boundaries

- L^(5/4) follows from volume times the assumed/known 3/4 thermal window; square-site use remains universality-based.
- Canonical copula universality is stronger than one-time rank-law universality and is not yet proved.
- Complement symmetry is exact only in the correctly typed self-matching or primal/matching pair setting.
