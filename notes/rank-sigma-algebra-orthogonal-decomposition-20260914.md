# Exact orthogonal decomposition of all rank information at balance

Date: 2026-09-14. Exact finite theorem for draft PR #773.

Let `X=r-1 in {-1,0,1}` and suppose the finite ensemble is at a rank balance root

    P(r=0)=P(r=2)=a,
    P(r=1)=r1=1-2a.                                      (1)

Define the centred even rank variable

    Y = X^2 - 2a.                                         (2)

Then

    E X=E Y=0,
    Var X = 2a,
    Var Y = 2a(1-2a),
    Cov(X,Y)=0.                                           (3)

Because every function of a three-valued variable lies in `span{1,X,X^2}`, the two orthogonal variables X,Y span the **entire nonconstant rank sigma-algebra**.

## 1. Exact projector theorem for any observable

Let A be any square-integrable configuration observable. Write

    A_j = E[A|r=j],
    Delta_A = A_2-A_0,
    h_A = A_1 - (A_0+A_2)/2.                              (4)

The two rank covariances are

    Cov(A,X) = a Delta_A,
    Cov(A,Y) = -2a(1-2a) h_A.                             (5)

Therefore the conditional mean given rank is exactly

    E[A|r]-E A
      = Cov(A,X)/Var(X) * X
        + Cov(A,Y)/Var(Y) * Y.                            (6)

There is no residual rank-conditioned direction beyond these two columns.

Taking variances gives the exact ANOVA identity

    Var(E[A|r])
      = Cov(A,X)^2/Var(X)
        + Cov(A,Y)^2/Var(Y)                               (7)

or equivalently

    Var(E[A|r])
      = (a/2) Delta_A^2
        + 2a(1-2a) h_A^2.                                (8)

The law of total variance becomes

    Var(A)
      = E Var(A|r)
        + (a/2) Delta_A^2
        + 2a(1-2a) h_A^2.                                (9)

This is exact at every finite balanced system.

## 2. Apply the theorem to total occupation K

For A=K, use

    g1 = E[K|2]-E[K|0],
    h1 = E[K|1]-(E[K|0]+E[K|2])/2.                        (10)

Since the unconditional homogeneous Bernoulli law has

    Var(K)=N p(1-p),                                      (11)

we obtain

    N p(1-p)
      = E Var(K|r)
        + (a/2) g1^2
        + 2a(1-2a) h1^2.                                 (12)

The odd/topological explained fraction is

    R_odd^2 = a g1^2/[2 N p(1-p)] = rho_KX^2,             (13)

and the even/extreme-vs-rank1 fraction is

    R_even^2 = 2a(1-2a)h1^2/[N p(1-p)].                  (14)

Their sum is **all** variance in K explained by knowing the rank:

    Var(E[K|r])/Var(K) = R_odd^2+R_even^2.                (15)

In a self-matching system at its exact symmetric point, h1=0 and only the odd channel remains.

Under the square-site near-critical hypotheses used elsewhere, `g1~L^(3/4)` while h1 is an irrelevant matching-asymmetry correction. Hence the total rank-explained fraction of K is predicted to vanish at least as `L^-1/2`, with leading term (13).

## 3. Small-information expansion

If, in addition, the conditional laws of K given rank share the same central-limit variance and differ only by small standardized mean shifts, the mutual information between rank and K has the standard local-mixture expansion

    I(r;K)
      = Var(E[K|r])/[2 Var(K)] + lower terms             (16)

in nats. At a matching-symmetric centre where h1 is negligible,

    I(r;K) ~ rho_KX^2/2 = O(L^-1/2).                     (17)

Equation (16) is **conditional on the common-variance/local-Gaussian approximation**, not part of the exact theorem. The exact content is the variance decomposition (12)--(15).

This supports a specific negative calibration for Rao--Blackwell design: global K may be useful at small sizes, but it cannot carry an asymptotically finite fraction of rank information if the critical scaling hypothesis holds. Rich cut/connectivity geometry can still carry substantial information.

## 4. Two-source exponential family and Fisher orthogonality

Introduce two finite topology sources,

    exp[s X + t X^2].                                     (18)

In canonical coordinates `(b,d)` with `d=log c`, the source acts exactly as

    b -> b-s,
    d -> d-t.                                             (19)

Thus b and d are affine natural coordinates for the full interior rank simplex, with sufficient statistics `-X` and `1_{r=1}` (equivalently X and X^2, up to signs and constants).

At balance, the Fisher cross term vanishes exactly by (3). The odd and even topology sources are therefore statistically orthogonal there.

The two-source partition function is

    Z(s,t)
      = P0 e^(-s+t)+P1+P2 e^(s+t).                       (20)

At balance,

    Z(s,t)=1-2a+2a e^t cosh s.                            (21)

Its Hessian at the origin is precisely the rank Fisher matrix. The discriminant locus of its static complex zeros is a property of this three-state closure family, not an operator Jordan diagnostic.

## 5. Relation to B_even and source dictionaries

The even rank statistic Y (or X^2) is the canonical **rank-even closure source**. It should not be silently identified with the `B_even` measure tangent in PR #746: that work already demonstrates typed even and ambient-topology contributions which need not coincide with a function of rank alone.

The useful statement is narrower:

> Within the information contained in rank itself, X and Y are a complete orthogonal odd/even basis at balance.

Any additional source response, such as spatial structure, cut geometry or a genuine measure tangent, lies outside this two-dimensional rank sigma-algebra unless it is deterministically reducible to rank.

## 6. Calibration use for original-U and conditional integration

For a candidate score A, the exact quantities

    R_odd^2(A)  = Cov(A,X)^2/[Var(A)Var(X)],
    R_even^2(A) = Cov(A,Y)^2/[Var(A)Var(Y)]                (22)

measure the fraction of its variance explained by the two orthogonal rank channels. If their sum is one, A contains no information beyond rank at first-order score level. If smaller, the residual is a genuine finite source direction, but still not a continuum field identification.

This is a cheap front-end diagnostic for #275/#578; it does not replace their full nuisance/source contracts.
