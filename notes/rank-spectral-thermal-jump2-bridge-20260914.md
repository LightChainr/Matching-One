# Spectral Dirichlet energy, thermal influence, and direct double birth are one defect

Date: 2026-09-14. Exact finite identity linking the Fourier/spectral, pivotal-boundary and persistence-process layers of draft #773.

## 1. One-site jump algebra

For the topological rank observable

    X=r-1,

a single-site insertion has

    Delta_v X in {0,1,2}.

Let

    beta(p)=P_p(Delta_v X=2)

for a fixed site v on a transitive torus. Then pointwise

    (Delta_v X)^2 = Delta_v X + 2 1_{Delta_v X=2},

hence

    E[(Delta_v X)^2]
      = E[Delta_v X] + 2 beta(p).                         (1)

The ordinary Russo derivative and the squared-influence/Dirichlet term are therefore not independent objects.

## 2. Global thermal versus spectral identity

Transitivity gives

    M'(p)=N E[Delta_v X].                                 (2)

For the p-biased Fourier spectral sample S_X,

    Var(X) E|S_X|
      = pq N E[(Delta_v X)^2].                            (3)

Combining (1)--(3),

    Var(X) E|S_X|
      = pq [M'(p)+2N beta(p)].                            (4)

Equivalently,

    [Var(X) E|S_X|]/[pq M'(p)]
      = 1 + 2 beta/[E Delta_v X].                        (5)

The second term is exactly the direct-jump share of the rank derivative. With the earlier notation

    omega_2
      = 2 beta / [alpha+2 beta+gamma]
      = 2 beta / E[Delta_v X],

we obtain the compact exact identity

    boxed: [Var(X) E|S_X|]/[pq M'] = 1+omega_2.          (6)

Thus the difference between spectral-sample mean size and the thermal influence is precisely the shared-boundary `0->2` channel.

## 3. Balance-root form

At a balance root `P0=P2=a`,

    Var(X)=2a.

Let

    g1=E[K|r=2]-E[K|r=0].

Since `pq M'=a g1`, (6) reduces to

    boxed: 2 E|S_X| / g1 = 1+omega_2.                    (7)

No microscopic probability coordinate remains in this ratio.

Exact small controls:

| model | L | 2 E|S_X| / g1 - 1 |
|---|---:|---:|
| square NN |3|0.10396624|
| square NN |4|0.07431879|
| triangular self-matching |3|0.06896552|
| triangular self-matching |4|0.03885759|

The square values reproduce the independently computed rank-jump-two contribution to M' up to the displayed numerical precision. The decreasing trend is only a tiny-size control, not a measured arm exponent.

## 4. Integrated process identity

For a random site permutation, let `D=J2-J1`. Earlier exact chain calculus gives

    P(D=0)=N int_0^1 beta(p) dp.                          (8)

Integrating the difference between squared and ordinary total influence therefore yields

    int_0^1 [N E(Delta_v X)^2 - M'(p)] dp
      = 2 P(D=0).                                        (9)

Since `M(1)-M(0)=2`,

    int_0^1 N E[(Delta_v X)^2] dp
      = 2 + 2 P(D=0).                                    (10)

But the right side is exactly the expected discrete quadratic variation of the monotone rank path:

    E sum_k (R_{k+1}-R_k)^2
      = 2 + 2 P(D=0).                                    (11)

Hence the same defect has three representations:

    direct shared-boundary flux beta
      <-> excess Fourier/Dirichlet energy over thermal influence
      <-> excess quadratic variation of the two-birth rank process.        (12)

## 5. Scaling interpretation

Suppose the ordinary one-threshold pivotal channel is four-arm and the direct double-birth channel beta is governed by a higher-arm fusion event. Then omega_2 measures the relative cost of that fusion without involving the matching-odd one-point cancellation.

Candidates:

    six-arm beta  -> omega_2 ~ L^(-5/3),
    eight-arm beta -> omega_2 ~ L^(-4).

This test belongs to the matching-even direct-jump geometry. A six-arm law here does not by itself contradict an eight-arm leading **matching-odd** correction in M(p_c).

If beta is asymptotically negligible relative to the ordinary pivotal boundary, (6) further predicts that the leading spectral-sample mean and thermal slope carry the same pivotal amplitude after the exact `pq/Var(X)` conversion.

## 6. Interfaces

- #769: report beta and omega_2 alongside the signed/absolute pair atlas.
- #784: use (6) as the strongest low-cost cross-check between the spectral-sample code and the pivotal code.
- #778/#788: use (9)--(11) to cross-check the simultaneous-birth atom and path quadratic variation.
- #768: do not use omega_2 alone as the matching-odd six-arm kill test; beta is matching-even.

## 7. Boundaries

Equations (1)--(12) are exact finite identities. The arm assignments in section 5 are scaling hypotheses. Squared influence, signed second derivative and matching-odd one-point corrections remain distinct observables.
