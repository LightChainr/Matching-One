# Matching duality for the full two-birth persistence process

Date: 2026-09-14. Exact pathwise coupling theorem for the rank-filtration process. Addendum to draft #773 / #778 / #776.

## 1. Continuous Uniform-label coupling

Assign iid `U_v~Uniform[0,1]` to all sites. For the primal graph G define

    B_p={v:U_v<=p},

and let `r_G(p)` be its ambient homology rank. The matching complement at the same p is

    W_p=V\B_p={v:U_v>p}.

Define reflected labels `V_v=1-U_v`. Then W_p is exactly the occupied set of the matching graph Ghat at parameter `q=1-p` under labels V.

Configuration-level duality gives pathwise

    r_G(p)+r_Ghat^V(1-p)=2.                              (1)

Let

    T1^G=inf{p:r_G(p)>=1},
    T2^G=inf{p:r_G(p)=2}.

Reading (1) at the two threshold crossings yields

    boxed: T1^G = 1-T2^Ghat,                             (2)
    boxed: T2^G = 1-T1^Ghat.                             (3)

The equalities hold pathwise in the coupled label realization, not merely in distribution.

Therefore the raw persistence gap

    G_p=T2-T1

is matching-even, while the midpoint

    C_p=(T1+T2)/2-1/2

is matching-odd under the paired coupling.

## 2. Discrete random-permutation version

Let `R_k^G` be the rank after the first k sites of a uniformly random permutation `pi` have been inserted. Use the reversed permutation for Ghat. Then

    R_k^G + R_(N-k)^Ghat = 2.                            (4)

With insertion indices

    J1=min{k:R_k>=1},
    J2=min{k:R_k=2},

one obtains the exact off-by-one relations

    boxed: J1^Ghat = N-J2^G+1,                           (5)
    boxed: J2^Ghat = N-J1^G+1.                           (6)

Hence

    D^Ghat=J2^Ghat-J1^Ghat=D^G,                          (7)

while

    J1^Ghat+J2^Ghat-(N+1)
      = -[J1^G+J2^G-(N+1)].                              (8)

So persistence length is exactly matching-even and the centered birth midpoint is matching-odd.

## 3. Self-matching consequences

For a self-matching model such as honest triangular-site percolation, the reflected-label process has the same law as the original process. Therefore

    (T1,T2) =_law (1-T2,1-T1),                           (9)

and

    (J1,J2) =_law (N-J2+1,N-J1+1).                       (10)

Consequences:

- `E[T1]+E[T2]=1` exactly;
- the raw midpoint has mean 1/2;
- `D=J2-J1` is invariant under the symmetry;
- any odd function of the centered midpoint has zero expectation;
- any joint copula estimator in #778 must obey the anti-diagonal reflection symmetry.

The symmetry does not imply midpoint and gap independence.

## 4. Canonical b-coordinate version

For a self-matching finite model, the canonical odd coordinate satisfies

    b(1-p)=-b(p).

Define

    B1=b(T1),
    B2=b(T2).

Then (9) becomes

    (B1,B2) =_law (-B2,-B1).                             (11)

Thus

    M_b=(B1+B2)/2

is odd and

    G_b=B1-B2>=0

is even. In particular

    E[M_b]=0,
    Cov(M_b,G_b)=0                                       (12)

whenever the required moments exist. Again, zero covariance is a symmetry consequence, not independence.

This is a parameterisation-free process regression for #776 and #778.

## 5. Two-time kernel symmetry

For the canonical interval-coverage kernel

    H(u,v)=P(B2<u<v<=B1),  u<v,

the self-matching symmetry gives

    H(u,v)=H(-v,-u).                                     (13)

The one-time rank-one curve is its diagonal shadow:

    P1(b)=P(B2<b<=B1).

Hence self-matching makes `P1(b)` exactly even, while the full process carries the stronger two-time anti-diagonal symmetry (13).

## 6. Interfaces

- #778: use (5)--(13) as hard checks for paired raw data and reconstructed canonical copulas.
- #776: triangular exact/transfer outputs should satisfy the discrete and continuous symmetry at every finite honest size.
- #768/#769: direct-double-birth events `D=0` are matching-even; observing a six-arm law for them does not by itself produce a matching-odd one-point correction.
- #786: shell flag enumerators should obey the induced reversal/complement reciprocity.

## 7. Boundaries

These are exact finite coupling identities. They do not identify the continuum law of the gap or midpoint and do not imply Markovianity of the rank process.
