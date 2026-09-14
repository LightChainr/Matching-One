# Exact bridge from random-permutation rank gaps to continuous birth-time gaps

Date: 2026-09-14. Exact finite coupling identities for draft PR #773.

Assign N sites independent Uniform(0,1) activation thresholds. Almost surely their order is a uniform random permutation, independent of the sorted threshold values

    0<U_(1)<...<U_(N)<1.

Let `R_k` be the ambient rank after the first k sites in that permutation have been activated. Define integer birth indices

    J1 = min{k:R_k>=1},
    J2 = min{k:R_k=2},
    D  = J2-J1 >= 0.                                      (1)

The continuous birth times are

    T1=U_(J1),
    T2=U_(J2).                                             (2)

A direct rank jump 0->2 at one insertion is exactly the event `D=0`.

## 1. Conditional gap law is Beta and depends only on D

The spacings

    U_(1), U_(2)-U_(1), ..., U_(N)-U_(N-1), 1-U_(N)

have the Dirichlet(1,...,1) law and are independent of the random permutation. Conditional on `D=d>0`, `T2-T1` is the sum of d exchangeable spacings. Therefore

    T2-T1 | D=d  ~  Beta(d, N+1-d).                       (3)

For d=0, `T2-T1=0` exactly.

Hence the full continuous birth-gap distribution is the mixture

    Law(T2-T1)
      = P(D=0) delta_0
        + sum_(d=1)^N P(D=d) Beta(d,N+1-d).               (4)

The integer gap histogram is therefore a sufficient statistic for the continuous gap law under the standard monotone uniform-threshold coupling.

## 2. All continuous gap moments from D

For rising factorial `(x)^(overline k)`,

    E[(T2-T1)^k | D=d]
      = d^(overline k)/(N+1)^(overline k).                 (5)

Thus

    E[(T2-T1)^k]
      = E[D^(overline k)]/(N+1)^(overline k).              (6)

In particular

    E[T2-T1]=E[D]/(N+1),                                  (7)

which is the permutation version of the one-time integral identity `int P1(p)dp`.

The order-statistic noise around `D/(N+1)` is explicit. If `D<<N`,

    Var(T2-T1|D)
      = D(N+1-D)/[(N+1)^2(N+2)]
      ~ D/N^2.                                             (8)

Under the standard critical prediction `D~L^(5/4)` with `N=L^2`, this conditional clock noise has standard deviation `L^(-11/8)`, parametrically smaller than the `L^(-3/4)` topological birth-gap scale. Thus the random-permutation integer gap is an asymptotically sharp proxy for the continuous p-gap under that hypothesis.

## 3. One-time rank polynomials give E[D]

After k insertions the active set is a uniform k-subset. If `c_{1,k}` is the number of k-site configurations of rank one, then

    P(R_k=1)=c_{1,k}/binom(N,k).                           (9)

For each permutation, the number of k values for which `R_k=1` is exactly D. Therefore

    E[D]
      = sum_(k=0)^N c_{1,k}/binom(N,k).                   (10)

Combining with (7) reproduces

    E[T2-T1]
      = [1/(N+1)] sum_k c_{1,k}/binom(N,k).               (11)

So #775's one-time rank-sector polynomials determine the first gap moment, while a saved integer D histogram determines the full continuous gap law by (4).

## 4. The zero-gap atom equals the integrated rank-jump-two pivotal measure

Fix a site v and an outside configuration A on the other N-1 sites. Let

    Delta_v X(A)=X(A union {v})-X(A).

The event `Delta_v X=2` is exactly a direct rank jump 0->2 when v is inserted after A and before every site outside A.

If `|A|=k`, the probability in a uniform random permutation that precisely A precedes v and all remaining sites follow v is

    k!(N-1-k)!/N!
      = 1/[N binom(N-1,k)].                                (12)

Summing over v and A yields

    P(D=0)
      = sum_v int_0^1 P_p(Delta_v X=2) dp.                (13)

On a transitive torus,

    P(D=0)
      = N int_0^1 P_p(Delta_0 X=2) dp.                    (14)

This identity is exact and places the single-site rank-jump-two fusion event directly inside the birth-gap law.

Since

    M'(p)=sum_v E_p[Delta_v X]
         =sum_v [P(Delta_v X=1)+2P(Delta_v X=2)],          (15)

and `M(1)-M(0)=2`, integration gives the companion sum rule

    sum_v int_0^1 P_p(Delta_v X=1) dp
       = 2[1-P(D=0)].                                      (16)

Thus the entire integrated topological rank increase splits exactly into separate +1 insertions and direct +2 insertions.

## 5. Exact L=3,4 direct-jump atom from the committed spine atlas

The rank-jump-two classifier in PR #773 stores fixed-site outside-config counts by occupation k. Equation (14) reduces to

    P(D=0)=sum_k n_k/binom(N-1,k),                         (17)

where `n_k` is the number of outside k-site configurations with `Delta_0 X=2`.

The exact controls are

    L=3: P(D=0)=3/35
                 =0.0857142857142857...,

    L=4: P(D=0)=2809/45045
                 =0.0623598623598624....                  (18)

Two sizes do not determine an exponent.

## 6. A clean fusion-exponent diagnostic

Suppose a fixed-site direct rank-jump-two event in the critical window has a bare polychromatic j-arm probability of order `L^(-alpha_j)`. The p-window contributing to (14) has width `L^(-y_t)`, while the prefactor N contributes `L^2`. Thus

    P(D=0)
      ~ L^[2-alpha_j-y_t].                                 (19)

Using `y_t=3/4`:

- six-arm fusion (`alpha_6=35/12`) predicts

      P(D=0) ~ L^(-5/3);                                  (20)

- eight-arm fusion (`alpha_8=21/4`) predicts

      P(D=0) ~ L^(-4).                                    (21)

This diagnostic concerns the **pivotal fusion mechanism of direct jumps**, not the matching-odd one-point correction `M_L(p_c)`. A six-arm result for `P(D=0)` can coexist with an eight-arm first nonzero one-point correction if the six-arm channel cancels from that different sector.

The L=3->4 effective slope is far too small-size to decide between (20) and (21); no exponent claim is made.

## 7. Interface to threshold-rank data

Any production or exact engine that already stores `(J1,J2)` or their integer gap D can reconstruct the continuous birth-gap distribution **without resimulating uniform activation times**. Use (4)--(6).

This is particularly relevant to threshold-rank histograms: a histogram of D is more informative for the birth-window problem than separately reporting only its mean or variance. The atom at D=0 is the integrated direct-jump pivotal observable (14).

If an archive stores only one-time rank counts by k and not the paired births, only E[D] is recoverable through (10); higher D moments and the atom D=0 are not determined by the one-time marginals alone.
