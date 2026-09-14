# Rank one as a random birth interval: exact static/dynamic identities

Date: 2026-09-14. Exact consequences of the monotone rank process; additive note for draft PR #773.

Use the standard monotone site coupling: assign each site an independent uniform threshold and open it when `p` exceeds that threshold. The ambient black homology rank is a nondecreasing process taking values 0,1,2. Define its two birth times

    T1 = inf{p: r(p)>=1},
    T2 = inf{p: r(p)=2},

with `T1<=T2` almost surely.

For each fixed p the three rank probabilities satisfy the exact dictionary

    P0(p)=P(T1>p),
    P1(p)=P(T1<=p<T2),
    P2(p)=P(T2<=p).                                      (1)

Thus the rank-one state is literally the random interval between the two homology births.

## 1. The declared F is the mixture birth CDF

Let

    F1(p)=P(T1<=p)=1-P0(p),
    F2(p)=P(T2<=p)=P2(p).

Then the existing normalized matching observable is exactly

    F(p)=(1+M(p))/2
        =[F1(p)+F2(p)]/2
        =E[r(p)]/2.                                      (2)

So the balance root `M=0` is a median of the equal mixture of the first- and second-birth laws. It is **not** separately a median of T1 or T2.

The missing information in F is precisely the rank-one interval probability P1: from `(F,P1)` one reconstructs the whole rank law,

    P0 = 1-F-P1/2,
    P2 = F-P1/2.                                         (3)

This is another way to see why a theorem about the balance root or the mixture CDF does not by itself determine the two births separately.

## 2. P1 is the coverage function of the random birth interval

For every configuration,

    1_{r(p)=1} = 1_{T1<=p<T2}.                            (4)

Integrating before expectation gives the exact mean-gap identity

    E[T2-T1] = int_0^1 P1(p) dp.                          (5)

More generally, for any integrable test function phi,

    E[ int_(T1)^(T2) phi(p) dp ]
       = int_0^1 phi(p) P1(p) dp.                         (6)

Taking `phi(p)=n p^(n-1)` yields

    E[T2^n-T1^n]
       = n int_0^1 p^(n-1) P1(p) dp.                      (7)

These are exact finite-size identities requiring no asymptotics or independence between the two births.

## 3. Higher gap moments require multi-time rank-one persistence

For a fixed interval `[T1,T2)`, its k-fold Cartesian volume is `(T2-T1)^k`. Hence

    E[(T2-T1)^k]
      = int_[0,1]^k
          P( r(p_1)=...=r(p_k)=1 )
        dp_1...dp_k,                                      (8)

where all ranks use the **same monotone coupling realization**.

Since rank is monotone, the joint event depends only on the extremes of the queried times:

    {r(p_i)=1 for all i}
      = {T1<= min_i p_i,  T2> max_i p_i}.                 (9)

Equivalently, after ordering `u=min p_i`, `v=max p_i`,

    E[(T2-T1)^k]
      = k(k-1) int_(0<u<v<1)
           (v-u)^(k-2)
           P(T1<=u, T2>v)
        du dv                                             (10)

for `k>=2`.

Thus a two-time rank-one persistence surface is sufficient to recover **all** moments of the birth gap. Static one-time P1 supplies only its first moment.

This gives a clean target for future exact/conditional-sampling work if gap fluctuations, rather than only marginal birth quantiles, become scientifically relevant.

## 4. Canonical coordinates interpret the plateau pointwise

In the softmax coordinates of the companion note,

    P1(p)= e^{d(p)} / [e^{b(p)}+e^{d(p)}+e^{-b(p)}].       (11)

At the balance root b=0,

    P1(p_*)=e^d/[2+e^d].                                  (12)

Hence d at the root is exactly the log-odds of the event that the balance point lies **between the two births** versus either extreme rank:

    d = log[P1/sqrt(P0P2)]
      = log[P1/P0]  at balance.                           (13)

A large positive d means that the root is located in a configuration-wise rank-one plateau with high probability, even though rank-0 and rank-2 probabilities are exactly balanced.

This is the static categorical content behind the root-versus-law separation.

## 5. A useful no-go and a useful sufficient package

A small mean gap alone does not force both birth laws to concentrate at the balance root: one may have `T1=T2` but a broad common random threshold. Likewise root convergence alone does not control the gap.

A sufficient package for both births to concentrate at pc is, for example,

1. the mixture law `F=(F1+F2)/2` concentrates at pc; and
2. `E[T2-T1]->0`.

Condition 1 alone already implies each marginal concentrates if T1<=T2 and both have equal total mass in the mixture; in the current manuscript this is expressed through the full-law theorem. Identity (5) then quantifies the disappearance of the intermediate rank-one interval rather than replacing that theorem.

Conversely, if `int P1` stays bounded away from zero along a sequence, the two births cannot both converge to the same deterministic limit.

## 6. Interface to existing computations

Rank-sector polynomial data (#775) immediately give P1(p) and therefore the exact mean gap (5) by polynomial/Beta integration, with no monotone trajectory simulation.

If

    P1(p)=sum_k c_{1,k} p^k(1-p)^(N-k),                   (14)

then

    E[T2-T1]
      = sum_k c_{1,k} Beta(k+1,N-k+1)
      = [sum_k c_{1,k}/binom(N,k)]/(N+1).                 (15)

The coefficients are configuration counts, so (15) is an exact rational number.

This gives #775 an additional zero-cost output: exact mean first-to-second homology-birth separation for each width.

The two-time surface in (10) is not determined by one-time rank polynomials and would require a joint monotone-coupling representation if higher gap moments are desired.
